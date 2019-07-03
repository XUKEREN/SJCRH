package GermlineSomatic;
# port of germline/somatic calling logic from Xiang Chen's
# /nfs_exports/genomes/1/PCGP/BucketIntermediate/XChen/Pallas/CaptureValidation.java
#
# MNE 1/2016
#

use strict;
use Carp qw(confess);
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);

use constant CALL_GERMLINE => "germline";
use constant CALL_SOMATIC => "somatic";
use constant CALL_UNKNOWN => "unknown";
use constant CALL_WILDTYPE => "wildType";
use constant CALL_UNCOVERED => "uncovered";

#use Text::NSP::Measures::2D::Fisher2::twotailed;
use Text::NSP::Measures::2D::Fisher::twotailed;
use List::Util qw(sum);

@GermlineSomatic::ISA = qw(Configurable Exporter);
@GermlineSomatic::EXPORT_OK = qw(
CALL_WILDTYPE
CALL_SOMATIC
CALL_UNCOVERED
);

use constant DEFAULT_SOMATIC_P_VALUE_CUTOFF => 0.05;

my $UNCOVERED_MAX_READS = 20;

use MethodMaker qw(
p_value_max
p_value
germline_only
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->p_value_max(DEFAULT_SOMATIC_P_VALUE_CUTOFF);
  $self->configure(%options);
  return $self;
}

sub call {
  my ($self, %options) = @_;
  my ($m_in_t, $r_in_t, $m_in_n, $r_in_n) = get_hash_option(
							    \%options,
					    "-mutant-in-tumor",
					    "-reference-in-tumor",
					    "-mutant-in-normal",
					    "-reference-in-normal"
							   );

  my $call;

  if ($m_in_t or $m_in_n) {
    #
    #  some observation of the variant allele.
    #  Not sure how well this will work for very low raw counts,
    #  i.e. not having a formal Bambino call; always appropriate??
    #
    my $npp = sum($m_in_t, $r_in_t, $m_in_n, $r_in_n);
    my $n11 = $m_in_t;
    my $n1p = $m_in_t + $r_in_t;
    my $np1 = $m_in_t + $m_in_n;

    my $p_value = calculateStatistic(
				     n11 => $n11,
				     n1p => $n1p,
				     np1 => $np1,
				     npp => $npp
				    );
    if (my $errorCode = getErrorCode()) {
      printf STDERR "ERROR: Fisher's Exact exception: %s\n", join " ", $errorCode, getErrorMessage(), $npp, $n11, $n1p, $np1, $m_in_t, $r_in_t, $m_in_n, $r_in_n;
      $self->p_value("");
      $call = CALL_UNKNOWN;
    } else {
      $self->p_value($p_value);

      if ($m_in_n * 4 > $r_in_n) {
	# mutant allele frequency too high in normal
	# (maybe better expressed fractionally?)
	$call = CALL_GERMLINE;
      } elsif ($m_in_t * $r_in_n < 10 * $r_in_t * $m_in_n) {
	# too-high tumor/normal ratio in tumor using swapped reference count (?)
	$call = CALL_GERMLINE;
#      dump_die(\%options, "TEST ME: 2nd germline call case $p_value");
	# need example to check this logic
      } elsif ($p_value < $self->p_value_max()) {
	# high-confidence somatic
	$call = CALL_SOMATIC;
      } else {
	# uncertain
	$call = CALL_UNKNOWN;
      }
    }
  } elsif ($self->germline_only ? $r_in_n <= $UNCOVERED_MAX_READS :
	   $r_in_t <= $UNCOVERED_MAX_READS) {
    # total coverage in tumor validation BAM too low
    $call = CALL_UNCOVERED;
  } else {
    # NOT validated
    $call = CALL_WILDTYPE;
  }

  return $call;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
