package UniProtIDMapping;
# parser/index for UniProt ID mapping files, e.g.
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
# MNE 3/2014

use strict;

use Configurable;
use FileUtils qw(universal_open);

@UniProtIDMapping::ISA = qw(Configurable);

use MethodMaker qw(
	restrict_fields
        index_fields

all
index
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->restrict_fields([]);
  $self->index_fields(
    [
     "GeneWiki",
     "RefSeq_NT",
    ]
      );
  $self->configure(%options);
  return $self;
}

sub add_index_field {
  my ($self, $field) = @_;
  push @{$self->index_fields}, $field;
}

sub parse {
  my ($self, %options) = @_;
  my $file = $options{"-file"} || die "-file";
  my $fh = universal_open($file) || die;

  my %restrict = map {$_, 1} @{$self->restrict_fields};
  my %index_fields = map {$_, 1} @{$self->index_fields};
  my $restrict = %restrict;

  my %all;
  my %index;

  while (<$fh>) {
    chomp;
    my @f = split /\t/, $_;
    die unless @f == 3;
    my ($uniprot_id, $key, $value) = @f;
    next if $restrict and not $restrict{$key};
    push @{$all{$uniprot_id}{$key}}, $value;
    $index{$key}{$value}{$uniprot_id} = 1 if $index_fields{$key};
  }
  close $fh;

  $self->all(\%all);
  $self->index(\%index);
}

sub find {
  my ($self, %options) = @_;
  my $field = $options{"-field"} || die "-field";
  my $value = $options{"-value"};
  die "-value" unless defined $value;

  my @hits;
  if (my $hits = $self->index()->{$field}{$value}) {
    @hits = sort keys %{$hits};
  }
  return \@hits;
}

sub get {
  my ($self, $id) = @_;
  return $self->all->{$id};
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
