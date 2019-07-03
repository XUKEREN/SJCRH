package TARTANUtils;
# describe me

use strict;
use Configurable;
use MiscUtils qw(get_hash_option);

use File::Basename;
use Exporter;

@TARTANUtils::ISA = qw(Configurable Exporter);
@TARTANUtils::EXPORT_OK = qw(
tartan_genome_put_helper
tartan_clone_env
);

use MethodMaker qw(
	
		  );

sub tartan_genome_put_helper {
  #
  #  "tartan put" a TARTAn output directory to genome config index location
  #  - TO DO: look up path from config variable
  #
  my (%options) = @_;
  my $tartan_out = $options{"-out"} || die "-out [tartan output dir]";
  my $genome = $options{"-genome"} || die "-genome";
  my $subdir = $options{"-subdir"} || die "-subdir";
  die "$tartan_out must be tartan output dir" unless basename($tartan_out) eq "output";
  die "where is $tartan_out" unless -d $tartan_out;

  my $genome_dir = get_genome_dir($genome);
  my $target = sprintf '%s/%s', $genome_dir, $subdir;
  my $cmd = sprintf 'tartan put %s %s', $tartan_out, $target;
  printf "command: %s\n", $cmd;
  printf "OK? [y/n]: ";
  my $ok = <STDIN>;
  chomp $ok;
  if ($ok and lc($ok) eq "y") {
    system $cmd;
    die if $?;
  } else {
    print STDERR "quitting\n";
  }
}

sub tartan_clone_env {
  #
  # import a tartan item from one environment to the other (e.g. dev to prod)
  # 
  my (%options) = @_;
  my $genome = get_hash_option(\%options, "-genome");
  my $subdir = get_hash_option(\%options, "-subdir");
  my $from = get_hash_option(\%options, "-from");
  my $to = get_hash_option(\%options, "-to");
  my $genome_dir = get_genome_dir($genome);
  my $template = sprintf '%s/%s', $genome_dir, $subdir;
  $template =~ s/(\w+)\/tartan\/index/%s\/tartan\/index/ || die;
  my $current_env = $1;
  die sprintf 'current environment is %s, must be %s', $current_env, $to unless $current_env eq $to;

  my $path_from = sprintf $template, $from;
  my $path_to = sprintf $template, $to;
  die "where is $path_from" unless -e $path_from;
  if (-d $path_from) {
    my $cmd = sprintf 'tartan import %s %s', $path_from, $path_to;
    die $cmd;
  } else {
    die "only directories implemented currently";
  }

}

sub get_genome_dir {
  my ($genome) = @_;
  my $tr = $ENV{TARTAN_ROOT} || die "no TARTAN_ROOT";
  my $genome_dir = sprintf '%s/index/reference/Homo_sapiens/%s', $tr, $genome;
  unless (-d $genome_dir) {
    # GRCh38 uses "no_alt" suffix even though main config name is GRCh38  :/
    my @h = glob($genome_dir . "*");
    ($genome_dir) = @h if @h == 1;
  }
  die "where is $genome_dir" unless -d $genome_dir;
  return $genome_dir;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
