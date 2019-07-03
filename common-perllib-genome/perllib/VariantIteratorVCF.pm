package VariantIteratorVCF;
# provide iterator for a VCF file returning a single record
# for each variant allele, optionally promoting to standardized record.
# MNE 5/2015

use strict;
use Configurable;
use Exporter;

use FileHandle;
use Vcf;
use VCFUtils;
use MiscUtils qw(dump_die);
use Variant;

@VariantIteratorVCF::ISA = qw(Configurable Exporter);
@VariantIteratorVCF::EXPORT_OK = qw();

use MethodMaker qw(
	file
vcf
vcfu
        queue
standardize
is_secondary
row_main
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

  my $vcf;
  if ($vcf = $self->vcf()) {
    # already initialized, e.g. via TabixFile.pm
  } elsif (my $file = $self->file) {
    if ($file =~ /\.bz2/i) {
      my $fh = new FileHandle();
      $fh->open("bzip2 -dc $file|") || die;
      $vcf = Vcf->new(fh => $fh);
    } elsif ($file eq "-") {
      # STDIN
      $vcf = Vcf->new(fh => *main::STDIN);
    } else {
      my $standard;
      $standard = 1 if $file =~ /\.vcf$/i;
      $standard = 1 if $file =~ /\.gz$/i;

      if ($standard) {
	# Vcf.pm apparently depends on file extensions to detect
	# compression, FAIL on dnanexus
	$vcf = Vcf->new(file => $file);
      } else {
	my $info = `file $file`;
	# bleh: perl equivalent?

	if ($info =~ /gzip/) {
	  my $fh = new FileHandle();
	  $fh->open("gzip -dc $file|") || die;
	  $vcf = Vcf->new(fh => $fh);
	} else {
	  # assume standard file even if not named .vcf (e.g. dnanexus)
	  $vcf = Vcf->new(file => $file);
	}

      }
    }
    $vcf->parse_header();
    $self->vcf($vcf);
  }
  $self->vcfu(new VCFUtils("-vcf" => $vcf));
}

sub get_next() {
  my ($self, %options) = @_;

  my $vcf = $self->vcf() || die;
  my $vcfu = $self->vcfu();

  my $queue = $self->queue();
  $self->is_secondary(0);
  if ($queue and @{$queue}) {
    # queued records from last call exist
#    printf STDERR "NOTE: queued record!!\n";
    $self->is_secondary(1);
  } else {
    my $row_main = $vcf->next_data_hash();
    $self->row_main($row_main);
    if ($row_main) {
      $queue = $vcfu->get_alt_rows("-hash" => $row_main);
      $self->queue($queue);
      # split row into one row for each alternate allele
      # with appropriately parsed-out alt info
    }
  }

  my $result;
  if ($queue and @{$queue}) {
    $result = shift @{$queue};

    if ($self->standardize()) {
      my $v = new Variant();
      $v->exception_warn(1);
      $v->import_vcf_row("-row" => $result);
      if ($self->standardize() == 2) {
	$v->{raw_row} = $result;
      }

      $result = $v;
    }

  }
  return $result;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
