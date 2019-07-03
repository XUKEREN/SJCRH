package ReviewableReportable;
# is a particular gene symbol germline-reviewable or germline-reportable?
# MNE 7/2018

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use FileUtils qw(read_simple_file);
use TdtConfig;
use GeneSymbolMapper qw(new_gsm_lite);

@ReviewableReportable::ISA = qw(Configurable Exporter);
@ReviewableReportable::EXPORT_OK = qw();

use MethodMaker qw(
genome
gsm_reviewable
gsm_reportable
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
  my ($self, %options) = @_;
  my $genome = $self->genome || die "-genome";

  my $gsm_reviewable = new_gsm_lite();
  my $gsm_reportable = new_gsm_lite();

  # init reviewable genes:
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $ranges = $config_genome->{CLINCLS_CANCER_RELATED_GENES_FILE} || die;
  open(TMPR, $ranges) || die;
  my %reviewable;
  while (<TMPR>) {
    chomp;
    my @f = split /\t/, $_;
    $reviewable{$f[0]} = 1;
  }
  foreach my $g (keys %reviewable) {
    $gsm_reviewable->add_gene("-gene" => $g);
  }

  # init reportable genes:
  my $reportable = read_simple_file($config_genome->{CLINCLS_GERMLINE_REPORTABLE_GENES} || die);
  foreach my $g (@{$reportable}) {
    $gsm_reportable->add_gene("-gene" => $g);
  }

  $self->gsm_reviewable($gsm_reviewable);
  $self->gsm_reportable($gsm_reportable);

}

sub is_reviewable {
  my ($self, $gene) = @_;
  return $self->gsm_reviewable->find($gene);
}

sub is_reportable {
  my ($self, $gene) = @_;
  return $self->gsm_reportable->find($gene);
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
