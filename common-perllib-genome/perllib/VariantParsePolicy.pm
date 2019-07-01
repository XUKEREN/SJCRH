package VariantParsePolicy;
# standardized interpretation of parsing variants from flatfiles,
# including column names and insertion-handling policy
# MNE 6/2017
#
# TO DO:
# - SJ post input format

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use Variant;

use constant INSERTION_BASE_BEFORE => "before";
use constant INSERTION_BASE_AFTER => "after";

@VariantParsePolicy::ISA = qw(Configurable Exporter);
@VariantParsePolicy::EXPORT_OK = qw();

my @CLOPS = (
	     "-f-chr=s",
	     "-f-pos=s",
	     "-f-ra=s",
	     "-f-va=s",
	     "-insertion-base=s",
	     "-bambino",
	     "-sj-post",
	     "-gedi",
	    );

use MethodMaker qw(
	f_chr
f_pos
f_ra
f_va
insertion_base_after
flags
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  return $self;
}

sub method {
  my ($self, %options) = @_;
  my $option = get_hash_option(\%options, "-xxx");
}

sub get_command_line_parameters {
  return @CLOPS;
}

sub setup_check {
  my ($self) = @_;
  unless ($self->f_chr) {
    my $flags = $self->flags || die;
    my ($f_chr, $f_pos, $f_ra, $f_va);
    if ($flags->{bambino}) {
      $f_chr = "Chr";
      $f_pos = "Pos";
      $f_ra = "Chr_Allele";
      $f_va = "Alternative_Allele";
      $flags->{"insertion-base"} = "after" unless $flags->{"insertion-base"};
      # unless user has manually overridden, i.e. if using 
      # bambino column names but without insertion adjustment
    } elsif ($flags->{"sj-post"}) {
      $f_chr = "Chr";
      $f_pos = "WU_HG19_Pos";
      $f_ra = "ReferenceAllele";
      $f_va = "MutantAllele";
      # this is *usally* the "after" insertion style but not always
    } elsif ($flags->{gedi}) {
      $f_chr = "chromosome";
      $f_pos = "pos";
      $f_ra = "reference_allele";
      $f_va = "non_reference_allele";
      $flags->{"insertion-base"} = "after";
    } else {
      my $msg = <<'EOS';
Variant parsing policy not specified. For predefined formats:
   * use -bambino for Bambino format
     (Chr/Pos/Chr_Allele/Alternative_Allele; assumes -insertion-position after)
   * use -gedi for GeDI format
     (chromosome/pos/reference_allele/non_reference_allele);
      assumes -insertion-position after)
   * use -sj-post for SJ postprocessed format
     (Chr/WU_HG19_Pos/ReferenceAllele/MutantAllele)
 or, specify column names manually:
   * -f-chr CHROMOSOME
   * -f-pos POSITION
   * -f-ra REFERENCE_ALLELE
   * -f-va VARIANT_ALLELE
 plus:
   * -insertion-base [before|after]
     (for insertions, specifies whether the given position refers to the 
      base number before or after the inserted sequence)
EOS

      $f_chr = $flags->{"f-chr"} || die $msg;
      $f_pos = $flags->{"f-pos"} || die $msg;
      $f_ra = $flags->{"f-ra"} || die $msg;
      $f_va = $flags->{"f-va"} || die $msg;
    }
    my $insertion_base = $flags->{"insertion-base"} || die "-insertion-base [before|after]";
    # die $insertion_base
    die "-insertion-base must be 'before' or 'after'" unless $insertion_base eq INSERTION_BASE_BEFORE or $insertion_base eq INSERTION_BASE_AFTER;

    $self->insertion_base_after($insertion_base eq "after" ? 1 : 0);
    $self->f_chr($f_chr || die "-f-chr");
    $self->f_pos($f_pos || die "-f-pos");
    $self->f_ra($f_ra || die "-f-ra");
    $self->f_va($f_va || die "-f-va");
  }
}

sub get_variant_from_row {
  my ($self, $row) = @_;
  $self->setup_check();
  my $chr = $row->{$self->f_chr} || dump_die($row, "no chrom field");
  my $pos = $row->{$self->f_pos} || dump_die($row, "no position field");
  my $ra = $row->{$self->f_ra};
  my $va = $row->{$self->f_va};
  foreach ($ra, $va) {
    dump_die($row, "undef reference or variant allele") unless defined $_;
    $_ = "-" unless $_;
  }

  my $v = new Variant();
  $v->import_generic(
		     "-reference-name" => $chr,
		     "-base-number" => $pos,
		     "-reference-allele" => $ra,
		     "-variant-allele" => $va
		    );
  if ($v->is_insertion() and $self->insertion_base_after()) {
    # adjust Bambino-style base numbering to base before insertion
    $v->start($v->start - 1);
    $v->end($v->end - 1);
  }
  return $v;
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
