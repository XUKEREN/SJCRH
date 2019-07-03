package VCFLeftAlign;
# wrapper to bcftools norm

use strict;

use File::Basename;

use Configurable;
use Exporter;
use TdtConfig;
use Cluster;

@VCFLeftAlign::ISA = qw(Configurable Exporter);
@VCFLeftAlign::EXPORT_OK = qw();

use constant RAM_REQUIRED => 500;

use MethodMaker qw(
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
}

sub multi_allelic_split {
  my ($self, %options) = @_;
  my $infile = $options{"-vcf"} || die "-vcf";

  die "abandoned";

  my $outfile = basename($infile) . ".split.vcf";
  $self->outfile($outfile);

  my $cmd = sprintf 'java -Xmx2g -jar %s -R %s -T LeftAlignAndTrimVariants --variant %s -o %s',
    $jar,
      $fa,
	$infile,
	  $outfile;

  $self->command($cmd);
  $self->run_cluster() if $options{"-cluster"};

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
