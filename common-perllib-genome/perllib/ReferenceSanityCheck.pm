package ReferenceSanityCheck;
# sanity-check a given reference bases vs. reference sequence
# MNE 8/2013
# TO DO: updated version to seek() to correct position based on line length

use strict;
use Carp qw(confess);

use Configurable;
use Exporter;

@ReferenceSanityCheck::ISA = qw(Configurable Exporter);
@ReferenceSanityCheck::EXPORT_OK = qw();

use MethodMaker qw(
	fasta_dir
last_reference_name
last_reference_sequence

last_reference_nt
undef_missing
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  return $self;
}

sub find_fasta_file {
  my ($self, $ref_name) = @_;
  my $local_bn = $ref_name;
  $local_bn =~ s/chr//;
  my @try = ($local_bn, "chr" . $local_bn);
  my $found;
  foreach my $try (@try) {
    my $local_fn = sprintf '%s/%s.fa',
    ($self->fasta_dir || die "no fasta dir"),
    $try;
    if (not(-s $local_fn) and $ref_name eq "M") {
      $local_fn = sprintf '%s/MT.fa',
      ($self->fasta_dir || die "no fasta dir");
    }
    if (-s $local_fn) {
      $found = $local_fn;
      last;
    }
  }
  return $found;
}

sub load_reference_sequence {
  my ($self, %options) = @_;
  my $ref_name = $options{"-ref-name"} || confess "-ref-name";
  if (($self->last_reference_name || "") ne $ref_name) {

    my $found = $self->find_fasta_file($ref_name);

    if ($found) {
      my $sequence = "";
      printf STDERR "RSC: loading %s\n", $found;
      open(RSFA, $found) || die;
      while (<RSFA>) {
	chomp;
	next if /^>/;
	$sequence .= $_;
      }
      $self->last_reference_name($ref_name);
      $self->last_reference_sequence(\$sequence);
    } elsif ($self->undef_missing) {
      my $sequence = "";
      $self->last_reference_name($ref_name);
      $self->last_reference_sequence(\$sequence);
    } else {
      confess "can't find FASTA for $ref_name in " . $self->fasta_dir;
    }

  }
}

sub get_chunk {
  my ($self, %options) = @_;
  my $start = $options{"-start"} || die "-start";
  my $length = $options{"-length"} || die "-length";
  $self->load_reference_sequence(%options);
  my $seq_ref = $self->last_reference_sequence() || die;
  return substr($$seq_ref, $start - 1, $length);
}

sub get_reference_base {
  my ($self, %options) = @_;
  $self->load_reference_sequence(%options);
  my $ref_pos = $options{"-base-number"} || die "-base-number";
  my $seq_ref = $self->last_reference_sequence() || die;
  return substr($$seq_ref, $ref_pos - 1, 1);
}

sub check {
  my ($self, %options) = @_;
  my $user_ref_base = uc($options{"-ref-base"} || die "-ref-base");
  my $local_ref_base = uc($self->get_reference_base(%options));
  $self->last_reference_nt($local_ref_base);
  my $status;
  if ($user_ref_base eq $local_ref_base) {
    $status = 1;
  } elsif ($local_ref_base eq "N") {
    # local reference is masked
    $status = -1;
  } else {
    $status = 0;
  }

  return $status;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
