package TSOncoDB;
# database of tumor suppressor/LoF and oncogene/GoF gene annotations
# MNE 4/2018

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use TdtConfig;
use GeneSymbolMapper;

@TSOncoDB::ISA = qw(Configurable Exporter);
@TSOncoDB::EXPORT_OK = qw();

my $F_SYM = "sym";
my $F_IS_TS = "is_ts";
my $F_IS_ONCO = "is_onco";

use MethodMaker qw(
		    gsm
		    f_tsdb

		    db
restrict_source
exclude_source
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
  my $gsm = $self->gsm;
  my $config_species = TdtConfig::readConfig('species', "Homo_sapiens") || die "can't find hs species config";
  unless ($gsm) {
    my $f_hgnc = $config_species->{HGNC} || die "no HGNC config";

    $gsm = new GeneSymbolMapper(
      "-hgnc_file" => $f_hgnc,
      "-hgnc_lite" => 1,
      "-enable_entrez_gene" => 0,
	);
    $self->gsm($gsm);
  }

  my $f_tsdb = $self->f_tsdb;
  unless ($f_tsdb) {
    $f_tsdb = $config_species->{TS_ONCO_DB} || die "TS_ONCO_DB";
  }

  my $df = new DelimitedFile("-file" => $f_tsdb,
			     "-headers" => 1,
			     );
  my %db;
  my $restrict = $self->restrict_source();
  my $exclude = $self->exclude_source();
  while (my $row = $df->get_hash()) {
    my $sym = $row->{$F_SYM} || dump_die($row, "no $F_SYM");
    die "duplicate $sym" if $db{$sym};
    if ($restrict) {
      foreach my $f ($F_IS_TS, $F_IS_ONCO) {
	if (grep {$_ eq $restrict} split /,/, $row->{$f}) {
	  $row->{$f} = $restrict;
	} else {
	  $row->{$f} = "";
	}
      }
    } elsif ($exclude) {
      foreach my $f ($F_IS_TS, $F_IS_ONCO) {
#	my $before = $row->{$f};
	$row->{$f} = join ",", grep {$_ ne $exclude} split /,/, $row->{$f};
#	dump_die($row, "$f $before") unless $before eq $row->{$f};
      }
    }
    $db{$sym} = $row;
    $gsm->add_gene("-gene" => $sym);
  }
  $self->db(\%db);
}

sub find_row {
  my ($self, %options) = @_;
  my $gene_raw = $options{"-gene"} || die "-gene";
  my $gene_idx = $self->gsm->find($gene_raw);
  return $gene_idx ? ($self->db->{$gene_idx} || die) : undef;
}

sub is_lof {
  # loss-of-function / truncating
  my ($self, %options) = @_;
  my $result;
  if (my $row = $self->find_row(%options)) {
    $result = $row->{$F_IS_TS};
  }
  return $result;
}

sub is_gof {
  # gain-of-function / recurrent
  my ($self, %options) = @_;
  my $result;
  if (my $row = $self->find_row(%options)) {
    $result = $row->{$F_IS_ONCO};
  }
  return $result;
}

sub add_gene_info {
  my ($self, %options) = @_;
  my $gene = $options{"-gene"} || die;
  my $is_ts = $options{"-is-ts"};
  my $is_onco = $options{"-is-onco"};
  die unless defined $is_ts and defined $is_onco;
  my %r;
  $r{$F_SYM} = $gene;
  $r{$F_IS_TS} = $is_ts;
  $r{$F_IS_ONCO} = $is_onco;

  $self->db->{$gene} = \%r;
  # user entry overrides db (TO DO: option)
  $self->gsm->add_gene("-gene" => $gene);
}

sub get_oncogenes {
  my ($self) = @_;
  return [ sort map {$_->{$F_SYM}} grep {$_->{$F_IS_ONCO}} values %{$self->{db}} ];
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
