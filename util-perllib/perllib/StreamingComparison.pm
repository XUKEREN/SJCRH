package StreamingComparison;
# compare two sorted or semi-sorted data streams, optionally pruning 
# duplicates as we go

use strict;
use Exporter;
use FileHandle;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);

@StreamingComparison::ISA = qw(Configurable Exporter);
@StreamingComparison::EXPORT_OK = qw();

use MethodMaker qw(
		    sets
		    save_shared
		    prune_interval
		    count
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->save_shared(0);
  $self->prune_interval(1000);
  $self->configure(%options);
  $self->reset();
  return $self;
}

sub reset {
  my ($self) = @_;
  $self->count(0);
  $self->sets({});
}

sub add_entry {
  my ($self, %options) = @_;
  my $set = get_hash_option(\%options, "-set");
  my $key = get_hash_option(\%options, "-key");
  my $value = $options{"-value"};
  my $sets = $self->sets();
  if ($value) {
    $sets->{$set}{$key} = $value;
  } else {
    $sets->{$set}{$key}++;
  }
  my $count = $self->count + 1;
  $self->prune() if $count % $self->prune_interval == 0;
  $self->count($count);
}

sub prune {
  my ($self) = @_;
  unless ($self->save_shared()) {
    my $sets = $self->sets();
    my @sets = keys %{$sets};
    die "need exactly 2 sets" unless @sets == 2;
    my ($k1, $k2) = @sets;
    my $s1 = $sets->{$k1} || die;
    my $s2 = $sets->{$k2} || die;

    my @shared;
    foreach my $k (keys %{$s1}) {
      push @shared, $k if exists $s2->{$k};
    }

    foreach my $k (@shared) {
      delete $s1->{$k};
      delete $s2->{$k};
    }
    printf STDERR "pruned %d, %s=%d, %s=%d\n",
      scalar(@shared),
	$k1, scalar(keys %{$s1}),
	  $k2, scalar(keys %{$s2}) if @shared;
  }
}

sub finish {
  my ($self) = @_;
  $self->prune();
}

sub write_exclusive {
  my ($self, %options) = @_;
  my $set_key = get_hash_option(\%options, "-set");
  my $file = get_hash_option(\%options, "-file");
  my $df = $options{"-df"};
  my $fh = new FileHandle();
  $fh->open(">" . $file) || die "can't write to $file: $!";
  printf $fh "%s\n", $df->header_line() if $df;
  my $set = $self->sets->{$set_key} || die;
  foreach my $key (sort keys %{$set}) {
    printf $fh "%s", $set->{$key};
  }
  $fh->close();
}



1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
