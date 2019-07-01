package OMIM_mim2gene;
# provide OMIM IDs for genes.
# This is the only redistributable/public file in OMIM.

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use TdtConfig;
use GeneSymbolMapper qw(new_gsm_lite);

@OMIM_mim2gene::ISA = qw(Configurable Exporter);
@OMIM_mim2gene::EXPORT_OK = qw();

use MethodMaker qw(
		    gsm
		    mim2gene
		    gene2mim
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  $self->gsm(new_gsm_lite());
  $self->setup();
  return $self;
}

sub setup {
  my ($self, %options) = @_;

  my $f_mim2gene = $self->mim2gene();

  unless ($f_mim2gene) {
    my $config = TdtConfig::readConfig('species', "Homo_sapiens") || die "can't get Homo_sapiens species config";
    $f_mim2gene = $config->{OMIM_MIM2GENE} || die "no OMIM_MIM2GENE config var";
  }

  my $gsm = $self->gsm();
  my %gene2mim;

  open(OMIMTMP, $f_mim2gene) || die;
  my @h;
  while (<OMIMTMP>) {
    chomp;
    if (/^\#/) {
      if (/^\# MIM Number/) {
	s/^\#\s+//;
	@h = split /\t/, $_;
	foreach (@h) {
	  s/\s*\(.*//;
	}
      }
    } else {
      die unless @h;
      my %r;
      @r{@h} = split /\t/, $_, -1;
      my $type = $r{"MIM Entry Type"} || dump_die(\%r);
      next unless $type eq "gene";

      my $sym = $r{"Approved Gene Symbol"} || next;
      # some records don't have one, e.g. 102777

      my $mim_number = $r{"MIM Number"} || die "no MIM #";

      $gsm->add_gene("-gene" => $sym);
      $gene2mim{$sym}{$mim_number} = 1;

#      die sprintf '"%s"', $r{"Ensembl Gene ID"};
#      dump_die(\%r, $sym);
    }
  }
  $self->gene2mim(\%gene2mim);
}

sub find_gene {
  my ($self, $gene) = @_;
  my $result;
  if (my $omim_sym = $self->gsm->find($gene)) {
    $result = [ sort {$a <=> $b} keys %{$self->gene2mim()->{$omim_sym}} ];
  }
  return $result;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
