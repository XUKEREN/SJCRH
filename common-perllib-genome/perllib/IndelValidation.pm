package IndelValidation;
# basic "fuzzy" validation of indels from raw Bambino output
# MNE 9/2014
#
# - run bambino restricted to target region
# - parse results
# - search for hit: same type (ins/del), same size
# - allow wiggle room if event > 1?
# - catch ambiguous hits
#
# TO DO: integrate snpshot?

use strict;

use Configurable;
use Exporter;

use File::Basename;

use BambinoRun;
use FileUtils qw(read_simple_file);
use MiscUtils qw(dump_die);
use SampleName;

my $INDEL_SITE_WIGGLE_ROOM = 25;
# how far away from source site target is allowed
# HACK: not sure what's appropriate here

@IndelValidation::ISA = qw(Configurable Exporter);
@IndelValidation::EXPORT_OK = qw();

my %TYPE2TN = (
	       "D" => "T",
	       "R" => "T",
    );


use MethodMaker qw(
	genome
		 );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  return $self;
}

sub find {
  my ($self, %options) = @_;
  my $row = $options{"-row"};
  # configure search from flatfile row
  my ($allele_ref, $allele_var, $ref_name, $pos);
  if ($row) {
    $allele_ref = $row->{ReferenceAllele} || die;
    $allele_var = $row->{MutantAllele} || die;
    $ref_name = $row->{Chr} || die;
    $pos = $row->{WU_HG19_Pos} || die;
  }
  die unless $allele_ref and $allele_var and $ref_name and $pos;
  my $bam = $options{"-bam"} || die "-bam";

  my ($is_insertion, $is_deletion, $indel_size);
  if ($allele_ref eq "-") {
    $is_insertion = 1;
    $indel_size = length($allele_var);
  } elsif ($allele_var eq "-") {
    $is_deletion = 1;
    $indel_size = length($allele_ref);
  } else {
    die;
  }

  my %info = SampleName::parse(basename($bam));
  my $type = $info{type} || die;
  my $tn = $TYPE2TN{$type} || dump_die(\%info, "need TN translation for type $type");

  my @extra = (
	       "-chr" => $ref_name,
	       "-start" => $pos - 50,
	       "-end" => $pos + 50,
	       "-tn" => $tn,
	      );

  #
  #  run Bambino on target region using low profile
  #  (more sensitive):
  #
  my $br = new BambinoRun();
  my $genome = $self->genome || die "-genome";
  my $bambino = $br->run(
			 "-bam" => $bam,
			 "-bam-type" => "rnaseq",
			 "-genome" => $genome,
			 "-run-type" => "low_unpaired",
			 "-java" => "-Xmx9g",
			 "-extra" => \@extra
			);

  #
  #  filter results to possible matches:
  #
  my $rows_raw = read_simple_file($bambino, "-as-hash" => 1);
  my @filtered;
  foreach my $row (@{$rows_raw}) {
    my $type = $row->{Type} || die;
    my $ok = 1;
    # type check:
    $ok = 0 if $type eq "SNP";
    if ($is_insertion) {
      $ok = 0 unless $type eq "insertion";
    } elsif ($is_deletion) {
      $ok = 0 unless $type eq "deletion";
    } else {
      die;
    }

    # position check:
    # (some wiggle room allowed because postprocessing may adjust):
    my $called_pos = $row->{Pos} || die;
    
    my $distance = abs($called_pos - $pos);
    $ok = 0 if $distance > $INDEL_SITE_WIGGLE_ROOM;

    # size check:
    my $called_size = $row->{Size} || die;
    if ($called_size == $indel_size) {
      if ($indel_size == 1) {
	# if indel is size 1, verify bases match.
	# done because (a) small false positive single-base indel 
	# calls are more common and (b) for multi-base events,
	# bases themselves may shift due to mapping ambiguity.
	# Hopefully that doesn't apply to single-base events!
	if ($is_deletion) {
	  die "verify single-base deletion";
	} elsif ($is_insertion) {
#	  dump_die($row, "verify single-base insertion");
	  my $called_base = $row->{Alternative_Allele} || die;
	  if ($called_base ne $allele_var) {
	    dump_die($row, sprintf "check me: single-base insertion mismatch, expected %s got %s", $allele_var, $called_base) if $ok;
	    $ok = 0;
	  }
	} else {
	  die;
	}
      }
    } else {
      $ok = 0;
    }


    push @filtered, $row if $ok;
  }

  return \@filtered;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
