package OutfileStack;

use strict;
use File::Basename;

use MiscUtils qw(log_message);
use FileUtils qw(newer_than);
use Configurable;
use Exporter;

@OutfileStack::ISA = qw(Configurable Exporter);
@OutfileStack::EXPORT_OK = qw();

use MethodMaker qw(
        start_file
	use_basename

        root
        suffixes
delete_invalid
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->delete_invalid(1);
  $self->suffixes([]);
  $self->configure(%options);

  my $use_basename = $self->use_basename;
  die "-use_basename" unless defined $use_basename;
  my $sf = $self->start_file() || die "-start_file";

  my $root = $use_basename ? basename($sf) : $sf;
  $self->root($root);

  return $self;
}

sub add_level {
  my ($self, %options) = @_;
  my $suffix = $options{"-suffix"} || die;
  $suffix = "." . $suffix unless $suffix =~ /^\./;

  my $suffixes = $self->suffixes;
  my $current = $self->get_current_file;

  unless ($self->is_first_level()) {
    die "current stack file $current doesn't exist" unless -e $current;
  }

  push @{$suffixes}, $suffix;
  $self->current_valid();
  # if desired, autodelete stale output
  return $self->get_current_file();
}

sub get_current_file {
  my ($self) = @_;
  return sprintf "%s%s", $self->root, join "", @{$self->suffixes};
}

sub get_previous_file {
  my ($self) = @_;
  my $suffixes = $self->suffixes || die;
  return sprintf "%s%s", $self->root, join "", @{$suffixes}[0 .. $#$suffixes - 1];
}

sub clone {
  my ($self) = @_;
  my $new = new OutfileStack("-use_basename" => 1, "-start_file" => "bogus");
  %{$new} = %{$self};
  $new->suffixes([ @{$self->suffixes} ]);
  return $new;
}

sub current_valid {
  my ($self, %options) = @_;
  my $f_current = $self->get_current_file;
  my $f_prev = $self->get_previous_file;
  my $valid = -e $f_current;
  if (newer_than($f_prev, $f_current)) {
    $valid = 0;
    if ($self->delete_invalid) {
      unlink $f_current;
      die "can't unlink invalid file $f_current" if -e $f_current;
    }
  }

  return $valid;
}

sub add_level_and_run {
  my ($self, %options) = @_;
  $self->add_level(%options);
  my $outfile = $self->get_current_file();
  my $infile = $options{"-infile"};
  unless (defined $infile) {
    if ($self->is_first_level()) {
      $infile = $self->start_file() || die();
    } else {
      $infile = $self->get_previous_file();
    }
  }
  die unless $infile;

#  unless (-s $outfile) {
  unless ($self->current_valid()) {
    if (my $template = $options{"-template"}) {
      my $cmd = sprintf $template, $infile;
      log_message("running $cmd");
      system $cmd;
      die "$cmd exited with $?" if $?;
    } elsif (my $callback = $options{"-callback"}) {
      if (ref $callback eq "CODE") {
	# simple subroutine callback
	&$callback($self);
      } else {
	die "implement object method callback";
      }
    } else {
      die "specify -template or -callback";
    }
  }
#  die "where is $outfile" unless -s $outfile;
  die unless $self->current_valid();
}

sub is_first_level {
  my ($self) = @_;
  my $suffixes = $self->suffixes;
  return @{$suffixes} <= 1;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
