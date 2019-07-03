package FastVariantIteratorVCF;
# provide iterator for a VCF file returning a single standardized
# record for each variant allele.
# Simpler implementation avoiding Vcf.pm which seems very slow.
# MNE 7/2015

use strict;
use Configurable;
use Exporter;

use FileHandle;
use FileUtils qw(universal_open);
use MiscUtils qw(dump_die);
use Variant;

@FastVariantIteratorVCF::ISA = qw(Configurable Exporter);
@FastVariantIteratorVCF::EXPORT_OK = qw();

use MethodMaker qw(
	file
        fh
        queue
headers
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
  my ($self) = @_;

  my $file = $self->file || die "-file";
  my $fh = universal_open($file) || die;
  while (<$fh>) {
    if (/^\#/) {
      if (/\#CHROM/) {
	chomp;
	$self->headers([ split /\t/, $_ ]);
	last;
      }
    } else {
      die;
    }
  }
  $self->fh($fh);
}

sub get_next() {
  my ($self, %options) = @_;


  my $queue = $self->queue();
  if ($queue and @{$queue}) {
    # queued records from last call exist
#    printf STDERR "NOTE: queued record!!\n";
  } else {
    my $fh = $self->fh;
    my $line = <$fh>;
    chomp $line;
    my %all;
    @all{@{$self->headers}} = split /\t/, $line;

    $queue = [];
    $self->queue($queue);
    foreach my $alt (split /,/, $all{ALT}) {
      my %r;
      $r{"#CHROM"} = $all{"#CHROM"};
      $r{POS} = $all{POS};
      $r{REF} = $all{REF};
      $r{ALT} = $alt;
      push @{$queue}, \%r;
    }
  }

  my $result;
  if ($queue and @{$queue}) {
    $result = shift @{$queue};
    my $v = new Variant();
    $v->exception_warn(1);
    $v->import_vcf_row("-row" => $result);
    $result = $v;
  }

  return $result;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
