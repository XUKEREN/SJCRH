package MNVRemap;
# remap MNVs
#
# some MNVs may contain spurious indels, e.g. 13.28592635.ATGAT.GGACA:
#
# http://bamviewer-rt:8080/BAMViewer/aceview/splash?tumorname=/rgs01/resgen/prod/tartan/runs/RNA_mapping/QP5rSpP9/output/SJAML040523_D1/finish/SJAML040523_D1.bam&ref=hg19&region=13&center=28592635&fullPath=true
#

#
# STUB, TBD
#

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option);
use Variant;

@MNVRemap::ISA = qw(Configurable Exporter);
@MNVRemap::EXPORT_OK = qw();

use MethodMaker qw(
	
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

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
