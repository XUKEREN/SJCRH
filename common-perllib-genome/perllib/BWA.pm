package BWA;
# bwa wrapper, particularly for ad-hoc alignents

use strict;
use File::Basename;

use Configurable;
use Exporter;

@BWA::ISA = qw(Configurable Exporter);
@BWA::EXPORT_OK = qw();

use MethodMaker qw(
        force
	bwasw
        aln
        verbose
	bam
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->force(1);
  $self->configure(%options);
  return $self;
}

sub bwa {
  my ($self, %options) = @_;
  my $reference_file = $options{"-reference"} || die;
  my $fastq = $options{"-fastq"} || die;
  die unless $self->bwasw or $self->bwa;

  my $header_file = $self->generate_header($reference_file);
  $self->bwa_index($reference_file);
  # index reference

  my $sam_file = $fastq . ".sam";
  my $bam_raw = $fastq . ".raw.bam";
  my $bam_sorted = $fastq . ".sorted.bam";

  my $need_bwa = -s $sam_file ? 0 : 1;
  $need_bwa = 1 if $self->force();

  if ($need_bwa) {
    #
    # run bwa
    #
    my $cmd;
    if ($self->bwasw) {
      $cmd = sprintf 'bwa bwasw %s %s > %s',
	$reference_file, $fastq, $sam_file;
    } else {
      die "implement me";
    }
    $cmd .= " 2>/dev/null" unless $self->verbose;
    system $cmd;
    die "bwa error $?" if $?;
  }

  my $need_bam_raw = -s $bam_raw ? 0 : 1;
  $need_bam_raw = 1 if $self->force();
  if ($need_bam_raw) {
    #
    # merge SAM output and header into raw BAM file
    #
    my $cmd = sprintf 'samtools view -t %s -b %s > %s',
    $header_file, $sam_file, $bam_raw;
    system $cmd;
    die "samtools error $?" if $?;
  }

  my $need_bam_sort = -s $bam_sorted ? 0 : 1;
  $need_bam_sort = 1 if $self->force();
  if ($need_bam_sort) {
    #
    # sort results
    #
    my $sname = $bam_sorted;
    $sname =~ s/\.bam$//;
    # samtools appends .bam to input name  :/

    my $cmd = sprintf 'samtools sort %s %s', $bam_raw, $sname;
    system $cmd;
    die "samtools error $?" if $?;
  }

  $self->bam_index($bam_sorted);

  die unless -s $bam_sorted;

  return $self->bam($bam_sorted);
}

sub bam_index {
  my ($self, $bam) = @_;
  my $bai = $bam . ".bai";
  my $needed = -s $bai ? 0 : 1;
  $needed = 1 if $self->force;
  if ($needed) {
    system "samtools index $bam";
    die if $?;
  }
}

sub bwa_index {
  # index a reference file if not already
  my ($self, $reference_file) = @_;
  my $idx_file = $reference_file . ".bwt";
  my $needed = -s $idx_file ? 0 : 1;
  $needed = 1 if $self->force;
  if ($needed) {
    my $cmd = "bwa index $reference_file";
    $cmd .= " >/dev/null 2>&1" unless $self->verbose;
    system $cmd;
    die "error running $cmd: $!" if $?;
  }
}


sub generate_header {
  # generate BAM header file for reference FASTA file
  my ($self, $reference_file) = @_;
  my $header_file = $reference_file . ".header.tab";
  my $needed = -s $header_file ? 0 : 1;
  $needed =1 if $self->force();
  if ($needed) {
    open(TMPIN, $reference_file) || die;
    my $id;
    my %length;
    while (<TMPIN>) {
      chomp;
      if (/^>(\S+)/) {
	$id = $1;
      } else {
	$length{$id} += length($_);
      }
    }

    open(TMPOUT, ">" . $header_file) || die;
    foreach (sort keys %length) {
      printf TMPOUT "%s\n", join "\t", $_, $length{$_};
    }
    close TMPOUT;
  }
  die unless -s $header_file;
  return $header_file;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
