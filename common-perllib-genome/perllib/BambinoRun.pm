package BambinoRun;
# wrapper to standard Bambino system scripts

use strict;
use Configurable;
use Exporter;

@BambinoRun::ISA = qw(Configurable Exporter);
@BambinoRun::EXPORT_OK = qw();

use TdtConfig;
use Digest::MD5 qw(md5_hex);
use MiscUtils qw(dump_die);
use FileUtils qw(newer_than);

use MethodMaker qw(
	md5_cache

command_line
passive
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->md5_cache({});
  $self->configure(%options);
  return $self;
}

sub run {
  my ($self, %options) = @_;
  my $run_type = $options{"-run-type"} || die "-run-type [high|low|high_unpaired...]";
  my $bam_type = $options{"-bam-type"} || die "-bam-type [normal|rnaseq]";
  die "-bam-type must be normal or rnaseq" unless $bam_type eq "normal" or $bam_type eq "rnaseq";

  my $genome = $options{"-genome"} || die "-genome";

  my ($bam_d, $bam_g, $bam);
  if ($run_type eq "high" or $run_type eq "low") {
    # paired
    $bam_d = $options{"-bam-d"} || die "-bam-d";
    $bam_g = $options{"-bam-g"} || die "-bam-g";
  } else {
    $bam = $options{"-bam"} || die "-bam";
  }

  my $config_genome = TdtConfig::readConfig('genome', $genome) || die;

  my $genome_fa;
  if (my $fa = $ENV{FASTA_OVERRIDE}) {
    $genome_fa = $fa;
    printf STDERR "using FASTA %s\n", $fa;
  } else {
    $genome_fa = $config_genome->{FASTA} || die "no FASTA in config";
  }
  my $dbsnp_blob = $config_genome->{DBSNP};
  unless ($dbsnp_blob) {
    $dbsnp_blob = "NOT_AVAILABLE";
    $options{"-extra"} = ($options{"-extra"} || "") . " -no-dbsnp";
  }

  my $out_dir = ".";
  my @params;

  my $script_name;
  my $file_index;

  if ($run_type eq "high") {
    # $1 = NIB directory or indexed fasta file for the reference genome
    # $2 = dbsnp blob file
    # $3 = diagnosis/tumor bam
    # $4 = germline/normal bam
    # $5 = output directory
    # $6 = root output filename
    # $7 = type ("normal" (default) or "rnaseq")
    # $8,$9 = (optional) -java "JAVA_ARGS" to pass args to java (such as -Xmx)
    # $8/10... = other arguments passed to SAMStreamingSNPFinder
    $script_name = "snv_high_20_tn.sh";
    @params = (
	       $genome_fa,
	       $dbsnp_blob,
	       $bam_d,
	       $bam_g,
	       $out_dir,
	       "outfile_placeholder",
               $bam_type,
	      );
    $file_index = 5;
  } elsif ($run_type eq "low") {
    # snv_low.sh:
    # $1 = NIB directory or indexed fasta file for the reference genome
    # $2 = diagnosis/tumor bam
    # $3 = germline/normal bam
    # $4 = output directory
    # $5 = root output filename
    # $6 = type ("normal" (default) or "rnaseq")
    # $7,$8 = (optional) -java "JAVA_ARGS" to pass args to java (such as -Xmx)
    # $7/9... = other arguments passed to SAMStreamingSNPFinder
    $script_name = "snv_low_tn.sh";
    @params = (
      $genome_fa,
      $bam_d,
      $bam_g,
      $out_dir,
      "outfile_placeholder",
      $bam_type,
	);
    $file_index = 4;
  } elsif ($run_type eq "low_unpaired") {
    # snv_low_unpaired.sh:
    # NIBORFA=$1; shift
    # DBAM=$1; shift
    # OUT_DIR=$1; shift
    # OUTFILE=$1; shift
    # TYPE=$1; shift
    # if [ "$1" == "-java" ]; then JAVA_ARGS="$2"; shift 2; fi
    # ARGS="$@"
    $script_name = "snv_low_unpaired.sh";
    @params = (
      $genome_fa,
      $bam,
      $out_dir,
      "outfile_placeholder",
      $bam_type,
	);
    $file_index = 3;
  } elsif ($run_type eq "germline") {
    $script_name = "snv_germline.sh";
    @params = (
      $genome_fa,
      $dbsnp_blob,
      $bam,
      $out_dir,
      "outfile_placeholder",
      $bam_type,
	);
    $file_index = 4;
  } elsif ($run_type eq "high_unpaired") {
    $script_name = "snv_high_20_unpaired.sh";
    @params = (
      $genome_fa,
      $dbsnp_blob,
      $bam,
      $out_dir,
      "outfile_placeholder",
      $bam_type,
	);
    $file_index = 4;
  } else {
    die "unhandled run type $run_type";
  }
  die unless $script_name;

  push @params, "-java" => $options{"-java"} if $options{"-java"};
  # e.g. "-Xmx4g"

  if (my $extra = $options{"-extra"}) {
    # extra parameters passed to bambino
    if (ref $extra) {
      push @params, @{$extra};
    } else {
      push @params, $extra;
    }
  }


  my $md5 = md5_hex($run_type, @params);
  my $cache = $self->md5_cache();
#    die "duplicate" if $cache->{$md5};
  # disable: might happen if calling multiple times
  # for e.g. swarm-style run
  $cache->{$md5} = 1;
  # sanity check

  my ($output_file, $stdout_file, $stderr_file);
  if ($output_file = $options{"-outfile"}) {
    $stdout_file = $output_file . ".out";
    $stderr_file = $output_file . ".err";
  } else {
    my $label = $options{"-label"} || $md5;
    $output_file = sprintf 'bambino.%s.tab', $label;
    $stdout_file = sprintf 'bambino.%s.out', $label;
    $stderr_file = sprintf 'bambino.%s.err', $label;
  }

  $params[$file_index] = $output_file;
  printf STDERR "outfile: %s\n", $output_file;

  my $needed = $options{"-force"};

  if (-s $stderr_file) {
    $needed = 1 unless -s $output_file;
    foreach my $f_bam ($bam_d, $bam_g, $bam) {
#      printf "%s\n", join " ", $bam, $stderr_file;
      $needed = 1 if defined($f_bam) and newer_than($f_bam, $stderr_file);
    }
  } else {
    $needed = 1;
  }

  my $cmd = sprintf '%s %s',
    $script_name,
      join(" ", @params);

  unless ($self->passive) {
    $cmd .= sprintf ' > %s 2>%s', $stdout_file, $stderr_file;
  }

  $self->command_line($cmd);

  if ($needed and not($self->passive())) {
    printf STDERR "running: %s\n", $cmd;
    system $cmd;
    die "error" if $?;
  }

  return $output_file;

}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
