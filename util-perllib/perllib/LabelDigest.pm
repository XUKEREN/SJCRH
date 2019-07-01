package LabelDigest;
# create unique short versions of CamelCase or underscore_delimited labels
# use for e.g. VCF INFO tags
# MNE 8/2017

use strict;
use Configurable;
use Exporter;

@LabelDigest::ISA = qw(Configurable Exporter);
@LabelDigest::EXPORT_OK = qw();

use MethodMaker qw(
	max_length
unique_counter
unique_out
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->max_length(2);
  $self->unique_counter({});
  $self->unique_out({});
  $self->configure(%options);
  return $self;
}

sub get_brief_label {
  my ($self, %options) = @_;
  my $label = $options{"-label"} || die "-label";
  my $max_length = $self->max_length || die;
  my $unique = $self->unique_counter();
  my $unique_out = $self->unique_out();
  my @things;
  if ($label =~ /_/) {
    # underscore_delimited
    @things = map {substr($_, 0, 1)} split /_+/, $label;
  } elsif ($label =~ /[A-Z][a-z]+[A-Z]/) {
    # CamelCase
    $label =~ s/[a-z]+//g;
    @things = split //, $label;
  } else {
    my $len = length($label);
    for (my $i = 0; $i < $max_length and $i < $len; $i++) {
      push @things, substr($label, $i, 1);
    }
  }

  my @out;
  for (my $i = 0; $i < $max_length and $i < @things; $i++) {
    push @out, $things[$i];
  }
  my $key = join "", @out;
  my $count = ($unique->{$key} || 0);
  $unique->{$key} = ++$count;
  $key .= $count if $count > 1;
#  printf STDERR "%s => %s, %d\n", $label, $key, $unique->{$key};

  if ($unique_out->{$key}) {
    die sprintf "duplicate %s for %s, last was %s", $key, $label, $unique_out->{$key};
  }

  $unique_out->{$key} = $label;
  

  return $key;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
