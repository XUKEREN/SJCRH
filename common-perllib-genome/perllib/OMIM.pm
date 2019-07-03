package OMIM;
# OMIM annotation

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use GeneSymbolMapper qw(new_gsm_lite);
use TdtConfig;

@OMIM::ISA = qw(Configurable Exporter);
@OMIM::EXPORT_OK = qw();

use MethodMaker qw(
		    gene2omim_genemap2
		    gsm
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

  my $config = TdtConfig::readConfig('species', "Homo_sapiens") || die "can't get Homo_sapiens species config";
  my $omim_dir = $config->{OMIM_DIR} || die "no OMIM_DIR";
  my $f_gm2 = $omim_dir . "/genemap2.txt.tab";
  die "where is $f_gm2" unless -s $f_gm2;
  my $df = new DelimitedFile(
			     "-file" => $f_gm2,
			     "-headers" => 1,
			     );
  my %gene2omim;
  while (my $row = $df->get_hash()) {
    my $sym = $row->{"Approved Symbol"};

    next unless $sym and $sym =~ /\w/;
    # only handle rows with approved symbols for now as these are the
    # basis for looking (including resolving ambiguity)

    push @{$gene2omim{$sym}}, $row;
    # might be more than one, e.g. NBPF15
  }

  my $gsm = new_gsm_lite();
  foreach (keys %gene2omim) {
    $gsm->add_gene("-gene" => $_);
  }

  $self->gene2omim_genemap2(\%gene2omim);
  $self->gsm($gsm);
}

sub find {
  my ($self, %options) = @_;
  my $gene_raw = $options{"-gene"} || die "-gene";
  my $result;
  if (my $gene = $self->gsm->find($gene_raw)) {
    $result = $self->gene2omim_genemap2->{$gene};
  }
  return $result;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
