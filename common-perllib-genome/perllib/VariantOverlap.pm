package VariantOverlap;
# overlap between any Variant.pm type and list of significant sites
# MNE 4/2017

use strict;
use Configurable;
use Exporter;
use GenomeUtils qw(cook_chromosome_name);

@VariantOverlap::ISA = qw(Configurable Exporter);
@VariantOverlap::EXPORT_OK = qw();

use MethodMaker qw(
		    db
		    insertions_search_before_and_after
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->db({});
  $self->insertions_search_before_and_after(1);
  $self->configure(%options);
  return $self;
}

sub add_site {
  my ($self, %options) = @_;
  my $ref = cook_chromosome_name($options{"-reference"} || die);
  my $pos = $options{"-position"} || die;
  $self->db()->{$ref}{$pos} = 1;
}

sub overlaps {
  my ($self, %options) = @_;
  my $v = $options{"-variant"} || die;
  my $start = $v->start();
  my $end = $v->end();
  $end++ if $v->is_insertion and $self->insertions_search_before_and_after();
  my $overlaps = 0;
  my $rn = $v->reference_name;
  if (my $ref_set = $self->db->{$rn}) {
    for (my $i = $start; $i <= $end; $i++) {
      if ($ref_set->{$i}) {
	$overlaps = 1;
	last;
      }
    }
  }
  return $overlaps;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
