package FuzzyVariantMatcher;
# fuzzy variant matching intended for comparing results w/possibly
# remapped coordinates, e.g. Bambino vs. GATK.
# MNE 6/2015
#
# - SNVs match SNVs/MNVs/indels
# - SNVs may match MNVs
# - however, SNV-on-SNV matches must be perfect if no other evidence involved
#   e.g.
#   - chr1.234.A.G vs. chr1.234.AG.TC: MATCH (e.g. Bambino SNV only
#                                             vs. GATK MNV)
#   - chr1.234.A.G vs. chr1.234.A.T: NO MATCH
#
# - when loading or querying MNVs, convert to SNVs
#

use strict;

use Configurable;
use Exporter;

use MiscUtils qw(dump_die);
use GenomeUtils qw(cook_chromosome_name);

@FuzzyVariantMatcher::ISA = qw(Configurable Exporter);
@FuzzyVariantMatcher::EXPORT_OK = qw();

use constant GRAVITY_SUB_VS_INDEL_NT => 0;
use constant GRAVITY_INDEL_VS_INDEL_NT => 3;

use MethodMaker qw(
		    snvs
                    ranges_ins
                    ranges_del
                    ranges_complex
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  $self->ranges_ins({});
  $self->ranges_del({});
  $self->ranges_complex({});
  $self->snvs({});
  return $self;
}

sub add_variant {
  my ($self, %options) = @_;
  my $var = $options{"-variant"} || die;

  my $chrom = cook_chromosome_name($var->reference_name);

  if ($var->is_substitution()) {
    my $hash = get_chr_hash($self->snvs(), $chrom);
    my $set = get_sub_snv_set($var);
    # convert MNVs to SNVs for matching purposes
    foreach my $v (@{$set}) {
      my $pos = $v->start;
      my $key = join ".", $pos, $v->reference_allele, $v->variant_allele;
      push @{$hash->{$key}}, $var;
      # reference is to SOURCE variant rather than split/set variant
    }
  } else {
    my $ranges;
    if ($var->is_insertion()) {
      $ranges = $self->ranges_ins();
    } elsif ($var->is_deletion()) {
      $ranges = $self->ranges_del();
    } elsif ($var->is_complex) {
      $ranges = $self->ranges_complex();
    } else {
      dump_die($var, "unhandled type " . $var->get_type());
    }
    my $hash = get_chr_hash($ranges, $chrom);
    for (my $i = $var->start; $i <= $var->end; $i++) {
      push @{$hash->{$i}}, $var;
    }
  }
}


sub get_chr_hash {
  # STATIC
  my ($parent, $chrom) = @_;
  my $hash = $parent->{$chrom};
  unless ($hash) {
    $hash = $parent->{$chrom} = {};
  }
  return $hash;
}

sub get_sub_snv_set {
  # STATIC
  my ($var) = @_;
  my @set;
  if ($var->is_mnv()) {
    my @ref = split //, $var->reference_allele;
    my @var = split //, $var->variant_allele;
    my $start = $var->start;
#    printf "%s\n", $var->get_key();
    for (my $i = 0; $i < @ref; $i++) {
      my $v = new Variant();
      $v->import_generic(
			 "-reference-name" => $var->reference_name,
			 "-base-number" => $start + $i,
			 "-reference-allele" => $ref[$i],
			 "-variant-allele" => $var[$i]
			);
#      printf "%s\n", $v->get_key();
      push @set, $v;
    }
  } else {
    # SNV
    @set = $var;
  }
  return \@set;
}

sub find_variant {
  my ($self, %options) = @_;
  my $q = $options{"-variant"} || die "-variant";

  my $chrom = cook_chromosome_name($q->reference_name);

  my @search_db;
  push @search_db, get_chr_hash($self->ranges_complex(), $chrom);
  # ALL variant types match vs. complex sites

  my $gravity;
  if ($q->is_deletion) {
    # indels match against basic type only
    $gravity = GRAVITY_INDEL_VS_INDEL_NT;
    push @search_db, get_chr_hash($self->ranges_del(), $chrom);
  } elsif ($q->is_insertion) {
    # indels match against basic type only
    $gravity = GRAVITY_INDEL_VS_INDEL_NT;
    push @search_db, get_chr_hash($self->ranges_ins(), $chrom);
  } elsif ($q->is_complex) {
    $gravity = GRAVITY_INDEL_VS_INDEL_NT;
  } elsif ($q->is_substitution) {
    $gravity = GRAVITY_SUB_VS_INDEL_NT;
  }
  die unless defined $gravity;

  #
  # check site vs. appropriate databases:
  #
  my @hits;
  my %saw;
  my $is = $q->start - $gravity;
  my $ie = $q->end + $gravity;

  foreach my $db (@search_db) {
    for (my $i = $is; $i <= $ie; $i++) {
      if (my $hs = $db->{$i}) {
	foreach my $h (@{$hs}) {
	  push @hits, $h unless $saw{$h};
	  $saw{$h} = 1;
	}
      }
    }
  }

  #
  # substitutions also get a separate, more specific check:
  #
  if ($q->is_substitution()) {
    #
    # when comparing sub vs. sub, also check exact bases.
    #
    my $db_snv = get_chr_hash($self->snvs, $chrom);
    my $set = get_sub_snv_set($q);
    foreach my $v (@{$set}) {
      my $pos = $v->start;
      my $key = join ".", $pos, $v->reference_allele, $v->variant_allele;
      if (my $hs = $db_snv->{$key}) {
	foreach my $h (@{$hs}) {
	  unless ($saw{$h}) {
	    push @hits, $h;
	    $saw{$h} = 1;
	  }
	}
      }
    }
  }

  return @hits ? \@hits : undef;
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
