package BAMInsertSize;

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die median decile);
use List::Util qw(min max);

@BAMInsertSize::ISA = qw(Configurable Exporter);
@BAMInsertSize::EXPORT_OK = qw();

use MethodMaker qw(
	max_records
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
#  $self->max_records(100000);
  $self->max_records(300000);
  $self->configure(%options);
  return $self;
}

sub get_insert_size {
  my ($self, %options) = @_;
  my $bam = get_hash_option(\%options, "-bam");
  my $cmd = sprintf 'samtools view %s|cut -f9|', $bam;
  my $max = $self->max_records;
  open(BISTMP, $cmd) || die;
  my @v;
  while (<BISTMP>) {
    chomp;
    next unless $_ > 0;
    # will be 0 for unmapped reads
    push @v, $_;
    last if @v >= $max;
  }
  die "samtools didn't return any results, executable path problem?" unless @v;

#  die median(\@v);
  return decile(\@v, 0.99);
  # prefer the longer sizes, ignoring the last 1% which may contain outliers
  # ...maybe round up to nearest 100 or 25?
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
