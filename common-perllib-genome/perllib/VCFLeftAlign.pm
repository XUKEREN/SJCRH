package VCFLeftAlign;
# wrapper to GATK LeftAlignAndTrimVariants

use strict;

use File::Basename;

use Configurable;
use Exporter;
use TdtConfig;
use Cluster;
use GenomeUtils qw(check_reference_name_compatibility);

@VCFLeftAlign::ISA = qw(Configurable Exporter);
@VCFLeftAlign::EXPORT_OK = qw();

use constant RAM_REQUIRED => 6000;

use MethodMaker qw(
genome
jar
fasta

command
ram
outfile
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self) = @_;

  $self->ram(RAM_REQUIRED);

  my @entries = grep {/gatk\/install/} split /:/, $ENV{PATH};
  die "can't identify GATK path entry" unless @entries == 1;
  my $jar = $entries[0] . "/GenomeAnalysisTK.jar";
  die unless -s $jar;
  $self->jar($jar);

  my $genome = $self->genome || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  $self->fasta($config_genome->{FASTA});
}

sub left_align {
  my ($self, %options) = @_;
  my $infile = $options{"-vcf"} || die "-vcf";

  my $jar = $self->jar() || die;
  my $fa = $self->fasta || die;

  check_reference_name_compatibility(
				     "-fasta" => $fa,
				     "-vcf" => $infile
				    ) || die sprintf "FATAL ERROR: reference sequence name incompatibility between %s and %s, LeftAlignAndTrimVariants will drop incompatible sequences.", $infile, $fa;

  my $outfile = basename($infile) . ".la.vcf";
  $self->outfile($outfile);

  my $cmd = sprintf 'java -Xmx2g -jar %s -R %s -T LeftAlignAndTrimVariants --variant %s -o %s',
    $jar,
      $fa,
	$infile,
	  $outfile;

  $cmd .= " --trimAlleles";
  # this appears to be safe

  if ($options{"-split"}) {
    # disabled by default as it strips genotype information!
    # however might be helpful for some VCFs not compatible with
    # bcftools, e.g. GiaB.
    $cmd .= " --splitMultiallelics";
    # use "bcftools norm -m-both" instead, BEFORE this is run.
  }

  if (1) {
    # WTF: some nodes apparently don't have a compatible java installation??
    # works OK on nodecn098, /hpcf/apps/java/jdk1.7.0_05/bin/java
    # java version "1.7.0_05"
    printf STDERR "DEBUG: PATH=%s\n", $ENV{PATH};
    printf STDERR "DEBUG: java path=%s\n", `which java`;
    printf STDERR "DEBUG: java version:%s\n", `java -version`;
  }

  $self->command($cmd);
  if ($options{"-cluster"}) {
    $self->run_cluster();
  } elsif ($options{"-run"}) {
    printf STDERR "running: %s\n", $cmd;
    system $cmd;
  }

  return $cmd;
}

sub run_cluster {
  my ($self, %options) = @_;

  my $c = new Cluster(
		      "-outfile" => $self->outfile,
		      "-project" => "PCGP",
		     );

  $c->node_class("");
  $c->memory_reserve_mb($self->ram());
  $c->memory_limit_mb($self->ram());
  $c->command($self->command() || die);
  $c->run();
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
