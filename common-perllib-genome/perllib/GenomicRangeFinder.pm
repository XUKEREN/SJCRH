package GenomicRangeFinder;
# 

use strict;
use Configurable;
use Exporter;

use Carp qw(cluck);

use GenomeUtils qw(cook_chromosome_name);
use BucketMap;
use MiscUtils qw(dump_die);

@GenomicRangeFinder::ISA = qw(Configurable Exporter);
@GenomicRangeFinder::EXPORT_OK = qw();

use constant DEFAULT_BUCKET_SIZE => 100000;

use constant STAMP_START => "__start__";
use constant STAMP_END => "__end__";

use MethodMaker qw(
	bucket_size
db
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->bucket_size(DEFAULT_BUCKET_SIZE);
  $self->db({});
  $self->configure(%options);
  return $self;
}

sub add {
  my ($self, %options) = @_;
  my $value = $options{"-value"} || die;

  my ($ref_cooked, $start, $end);

  if (my $range = $options{"-range"} || die) {
    # to do: add manually-parsed
    my @f = split /:/, $range;
    die unless @f == 2;
    my ($chrom, $r) = @f;
    if ($ref_cooked = cook_chromosome_name($chrom, "-undef-unknown" => 1)) {
      @f = split /\-/, $r;
      die unless @f == 2;
      ($start, $end) = @f;
      die if $end < $start;
      $value->{STAMP_START()} = $start;
      $value->{STAMP_END()} = $end;
    } else {
      cluck(sprintf "WARNING: invalid chrom in range %s\n", $range);
      return;
    }
  }

  die unless $ref_cooked and $start and $end;

  my $db = $self->db || die;
  my $bm = $db->{$ref_cooked};
  $bm = $db->{$ref_cooked} = new BucketMap("-chunk" => $self->bucket_size) unless $bm;

  $bm->add_range(
    "-start" => $start,
    "-end" => $end,
    "-value" => $value
      );
}

sub find {
  my ($self, %options) = @_;
  my $ref_cooked = cook_chromosome_name($options{"-chrom"} || die "-chrom");
  my $base = $options{"-base"};
  my $end = $options{"-end"};
  
  my @hits;

  if (my $bm = $self->db->{$ref_cooked}) {
    my $hits = $bm->fuzzy_find("-start" => $base,
			       "-end" => $end || $base);

    my %saw;
    foreach my $hit (@{$hits}) {
      my $hit_start = $hit->{STAMP_START()} || die;
      my $hit_end = $hit->{STAMP_END()} || die;

      next if $hit_end < $base;
      next if $hit_start > $base;
      next if $saw{$hit};
      $saw{$hit} = 1;
      push @hits, $hit;
    }
  }
  return \@hits;
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
