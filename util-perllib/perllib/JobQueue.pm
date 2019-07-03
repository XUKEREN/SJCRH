package JobQueue;
# simple wrapper to perform work in sequential stages
# MNE 7/2013

use strict;
use Configurable;

@JobQueue::ISA = qw(Configurable Exporter);
@JobQueue::EXPORT_OK = qw();

use Scalar::Util qw(blessed reftype);

use MethodMaker qw(
	levels
sleep_time
max_run_time
verbose
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->sleep_time(60);
  $self->max_run_time(60 * 60 * 24);
  $self->levels([]);
  $self->verbose(1);
  $self->configure(%options);
  return $self;
}

sub add_level {
  my ($self, %options) = @_;
  push @{$self->levels}, \%options;
}

sub run {
  my ($self) = @_;
  my $levels = $self->levels();
  my $sleep_time = $self->sleep_time();
  my $start_time = time;
  my $max_run_time = $self->max_run_time();

  my $level_number;
  foreach my $level (@{$levels}) {
    # foreach level in process:
    #
    #
    #  1. call some code to do some work:
    #
    $level_number++;
    my $callback = $level->{"-callback"};
    die "specify -callback" unless defined $callback;

    if ($callback) {
      # only if specified: some use cases might have started work already
      if (reftype($callback) eq "ARRAY") {
	if (blessed($callback->[0])) {
	  my ($object, $subname, @args) = @{$callback};
	  $object->$subname(@args);
	} else {
	  die "procedural callback";
	}
      } else {
	die "unhandled callback code";
      }
    }

    #
    #  2. wait until work is done:
    #
    while (1) {
      last if defined $level->{"-wait"} and $level->{"-wait"} eq "0";
      my $count_left = 0;

      my $progress;
      if (my $files = $level->{"-outfiles"}) {
	# a list of final outfiles: if they exist
	my $count_done = 0;
	foreach my $outfile (@{$files}) {
	  if (-s $outfile) {
	    $count_done++;
	  } else {
	    $count_left++;
	  }
	}
	$progress = $count_done / @{$files};
      } else {
	die "don't know how to check status";
      }
      
      if ($count_left == 0) {
	last;
      } else {
	my $elapsed = time - $start_time;
	if ($elapsed > $max_run_time) {
	  die "job is taking too long to finish, exiting!";
	}

	if ($self->verbose) {
	  printf STDERR "%s: processing level %d of %d...",
	  scalar(localtime), $level_number, scalar @{$levels};
	  printf STDERR "(%d left, %d%%)", $count_left, $progress * 100 if $progress;
	  print STDERR "\n";
	}

	sleep $sleep_time;
      }
    }



  }
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
