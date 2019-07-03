package DelimitedReadAhead;
# read a delimited file with a buffer of X rows.  Used to sync
# a very large delimited file with another data source in the same order.
# e.g. annovar->VCF
# MNE 6/2016

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use DelimitedFile;

my $VERBOSE = 0;

@DelimitedReadAhead::ISA = qw(Configurable Exporter);
@DelimitedReadAhead::EXPORT_OK = qw();

use MethodMaker qw(
	size
file
        hash_sub
df
db

obsolete_db
obsolete_prune_fraction_max
obsolete_prune_fraction_limit

obsolete_prune_max
obsolete_prune_limit
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->obsolete_prune_fraction_max(0.25);
  # let obsolete records grow to this fraction of the cache
  $self->obsolete_prune_fraction_limit(0.5);
  # when pruning triggered, eliminate this fraction of obsolete records
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self, %options) = @_;
  my $infile = $self->file || die "-fie";
  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			    );
  $self->df($df);
  $self->db({});
  $self->obsolete_db({});

  $self->obsolete_prune_max(int($self->size() * $self->obsolete_prune_fraction_max()));
  $self->obsolete_prune_limit(int($self->obsolete_prune_max() * $self->obsolete_prune_fraction_limit()));

  $self->fill_cache();
}

sub obsolete {
  my ($self, %options) = @_;
  my $key = $options{"-key"} || die "-key";
  my $seq = $options{"-sequence"};
  die unless defined $seq;
  my $db_obs = $self->obsolete_db() || die;
  $db_obs->{$seq} = $key;
  printf STDERR "obsolete $seq => $key\n" if $VERBOSE;

  my $pmax = $self->obsolete_prune_max();

  my $db = $self->db();
  if (scalar keys %{$db_obs} > $pmax) {
    my $plimit = $self->obsolete_prune_limit();
    my @sorted = sort {$a <=> $b} keys %{$db_obs};

    for (my $i = 0; $i < $plimit; $i++) {
      my $sequence = $sorted[$i];
      my $key = $db_obs->{$sequence};
      die "$key not in obs" unless $db_obs->{$sequence};
      delete $db_obs->{$sequence};
#      die "$key not in main" unless $db->{$key};
      # can happen if the same key appears twice
      # e.g. ClinVar:
      # 1       53676583        rs397509431     CAG     C
      # 1       53676583        rs398123153     CAG     C
      delete $db->{$key};
      print STDERR "delete $sequence => $key\n" if $VERBOSE;
    }
#    print STDERR "stop\n";
  }
}

sub fill_cache {
  my ($self) = @_;
  my $db = $self->db() || die;
  my $hash_sub = $self->hash_sub() || die "-hash_sub";
  my $wanted_size = $self->size() || die "-size";
  my $df = $self->df();

  while (scalar keys %{$db} < $wanted_size) {
    if (my $row = $df->get_hash()) {
      my $key = &$hash_sub($row);
      dump_die($row, "WARNING: duplicate $key", 1) if $db->{$key};
      $db->{$key} = $row;
      printf STDERR "store %s\n", $key if $VERBOSE;
    } else {
      last;
    }
  }
}

sub get_cache {
  return $_[0]->{db};
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
