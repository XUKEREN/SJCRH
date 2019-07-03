package SJBAMUtils;
# utilities for St. Jude BAM files

use strict;
use Configurable;
use Exporter;

use MiscUtils qw(dump_die);
use File::Basename;

@SJBAMUtils::ISA = qw(Configurable Exporter);
@SJBAMUtils::EXPORT_OK = qw(
bucket_bams
);

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  return $self;
}

sub bucket_bams {
  my (%options) = @_;
  my $files = $options{"-files"} || die "-files";
  my %bucket;
  foreach my $bam (@{$files}) {
    my %info = SampleName::parse(basename($bam), "WARN");
    next unless %info;
    # parsing problem, error message already printed
      
    my $key;
    if ($options{"-subject"}) {
      # bucket by BAM subject
      $key = $info{subject} || die;
    } elsif ($options{"-subject-and-type"}) {
      $key = join "_", @info{qw(subject type)};
    } else {
      # other types?
      die "-subject or ...?";
    }
    push @{$bucket{$key}}, $bam;
  }
  return \%bucket;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
