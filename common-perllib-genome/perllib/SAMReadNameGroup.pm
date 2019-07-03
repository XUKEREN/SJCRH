package SAMReadNameGroup;
# gather SAM records by read name

use strict;
use Carp qw(confess);
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use FileUtils qw(universal_open);

@SAMReadNameGroup::ISA = qw(Configurable Exporter);
@SAMReadNameGroup::EXPORT_OK = qw();

use MethodMaker qw(
file
fh
regexp
sam_headers
queue
last_id
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
  my ($self, %options) = @_;
  my $file = $self->file() || confess "-file";
  $self->fh(universal_open($file));
}

sub get_set {
  my ($self) = @_;
  my $regexp = $self->regexp || die "-regexp";
  my $fh = $self->fh || die;
  my @headers;
  my $last_id;
  my $this_id;

  my $queue;
  if ($queue = $self->queue()) {
    $last_id = $self->last_id;
  } else {
    $queue = [];
  }

  my $results;

  while (my $line = <$fh>) {
    chomp $line;
    if ($line =~ /^\@/) {
      push @headers, $line;
    } else { 
      my @row = split /\t/, $line;
      if ($row[0] =~ /$regexp/) {
	$this_id = $1 || die 'no $1 match';
	if (not(defined $last_id) or $last_id eq $this_id) {
	  # start of file or continuation of set
	  push @{$queue}, \@row;
	} else {
	  $results = $queue;
	  $self->queue([ \@row ]);
	  $self->last_id($this_id);
	  last;
	}
	$last_id = $this_id;
      } else {
	die "$line doesn't match $regexp";
      }
    }

  }

  $self->sam_headers(\@headers) if @headers;
  return $results;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
