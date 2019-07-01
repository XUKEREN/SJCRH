package DelimitedBinarySearch;
# binary search of delimited file data.
# uses configurable comparators to enable multi-field searches.
# MNE 2/2014

use strict;
use Configurable;
use Exporter;

use POSIX qw(SEEK_SET);
use FileHandle;
use Carp qw(confess);

use MiscUtils qw(dump_die);

use constant DELIMITERS => ("\t", ",");

@DelimitedBinarySearch::ISA = qw(Configurable Exporter);
@DelimitedBinarySearch::EXPORT_OK = qw();

use MethodMaker qw(
	file
fh
is_headered
column_number
delimiter
border_start
border_end
max_line_length
blankify
headers
hashify_headers
skip_comment_lines
line_parser_callback

final_seek_start
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->max_line_length(1000);
  $self->configure(%options);
  return $self;
}

sub find {
  my ($self, %options) = @_;
  my $file = $self->file() || die;
  my $comparators = $options{"-comparators"} || die "-comparators";
  # fugly; maybe break out into skeleton/query style and populate
  # with user variables??

  my @ccode;
  my $is_hash;
  foreach my $comparator (@{$comparators}) {
    my $type = $comparator->{type} || die;
    my $code;
    my $cval = $comparator->{value};
    die unless defined $cval;

    if (exists $comparator->{column_number}) {
      #
      #  comparison based on column indices
      #
      my $cnum = $comparator->{column_number} || die;
      $cnum--;
      # convert to index
      if ($type eq "string") {
	$code = sprintf '$row->[%d] cmp "%s"', $cnum, $cval;
      } elsif ($type eq "number") {
	$code = sprintf '$row->[%d] <=> %s', $cnum, $cval;
      } else {
	die;
      }
    } else {
      #
      #   comparison based on column labels
      #
      my $cname = $comparator->{column_name} || die "need column number or name";
      $is_hash = 1;
      if ($type eq "string") {
	$code = sprintf '$row->{"%s"} cmp "%s"', $cname, $cval;
      } elsif ($type eq "number") {
	$code = sprintf '$row->{"%s"} <=> %s', $cname, $cval;
      } else {
	die;
      }
    }
    push @ccode, $code || die;
  }
  my $comparator_code = sprintf '$comparison = %s;', join ' || ', @ccode;

  my @cidx;
  my @cnames;
  if ($is_hash) {
    @cnames = map {$_->{column_name}} @{$comparators};
  } else {
    @cidx = map {$_->{column_number} - 1} @{$comparators};
  }

  my $fh = $self->fh();
  unless ($fh) {
    $fh = new FileHandle();
    $fh->open($file) || die;
    $self->fh($fh);
  }

  my $border_start = $self->border_start();
  my $is_headered = $self->is_headered();
  if ($is_headered) {
    my $first = $self->parse_line(scalar <$fh>);
    $border_start = $fh->tell() unless defined $border_start;
  }
  $border_start = 0 unless defined $border_start;

  my $delimiter = $self->delimiter();


  my $border_end = $self->border_end() || -s $file;

  my $first_data_line_pos = $border_start;
  # hack

  my $start_pos = $border_start;
  my $end_pos = $border_end;

  my $final_seek_start;

  my $seek_no = 1;

  my $stop_gap = $self->max_line_length * 3;

  my $DEBUG = $options{"-verbose"};

  # binary search: fugly, should be abstracted
  my $last_pos;
  my $comparison;

  my $wanted_value = join ".", map {$_->{value}} @{$comparators};

  while (1) {
    my $chunk = $end_pos - $start_pos + 1;
    my $pos = $start_pos + int($chunk / 2);
    if ($last_pos and $pos == $last_pos) {
      # if resolves to same seek position, just use range start.
      # maybe this can happens if configured line length is lower
      # than actual line length?  i.e. we seek to the middle of
      # a longer-than-expected line, skip it because it's the first
      # line, and then wind up in the same place again.
      printf STDERR "  same seek start, bailing out!\n" if $DEBUG;
      # prevent infinite loop!
      $final_seek_start = $start_pos;
      last;
    }

    my ($row, $line_start_pos) = $self->seek_line($pos);
    eval $comparator_code;

    my $entry_value = get_entry_value($row, $is_hash, \@cnames, \@cidx);
    printf STDERR "DEBUG: %s %s\n", $comparator_code, $entry_value if $DEBUG;

    my $gap = $end_pos - $start_pos;

    printf STDERR "count:%d wanted:%s line:%s comparison:%d seek=%d start=%d end=%d line_pos=%d, gap=%d, stop_gap=%d\n",
    $seek_no, 
    $wanted_value,
    $entry_value,
    $comparison,
    $pos,
    $start_pos,
    $end_pos,
    $line_start_pos,
    $gap,
    $stop_gap if $DEBUG;

    $seek_no++;

    die "too many seeks" if $seek_no > 50;

    if ($line_start_pos == $first_data_line_pos) {
      print STDERR "  dbg: 1st data line\n" if $DEBUG;
      $final_seek_start = $pos;
      last;
    } elsif ($comparison > 0) {
      print STDERR "  dbg: pos after wanted\n" if $DEBUG;
      if ($gap < $stop_gap) {
	# site may not exist, e.g. 10.89692912.T.G
	print STDERR "dbg: small gap, stopping with buffer\n" if $DEBUG;
	$final_seek_start = $pos - $stop_gap;
	last;
      } else {
	$end_pos = $line_start_pos;
      }
    } elsif ($comparison == 0) {
      # entry matches, but may not be first line
      print STDERR "  dbg: match but maybe not 1st\n" if $DEBUG;
#      $end_pos = $line_start_pos;
      if ($gap < $stop_gap) {
	print STDERR "  dbg: small gap, stopping with buffer\n" if $DEBUG;
	$final_seek_start = $pos - $stop_gap;
	# buffer for safety
	last;
      } else {
	$end_pos = $pos;
      }
#      die "line matches";
    } elsif ($comparison < 0 and $gap < $stop_gap) {
      # this line is before the desired entry, but seek gap is close
      # enough to just stop and parse
      print STDERR "  dbg: close enough, stopping\n" if $DEBUG;
      $final_seek_start = $pos;
      last;
    } else {
      # too far before desired entry
      print STDERR "  dbg: pos before wanted\n" if $DEBUG;
      $start_pos = $line_start_pos;
    }

    $last_pos = $pos;
  }

  $final_seek_start = 0 if $final_seek_start < 0;
  $self->final_seek_start($final_seek_start);
  $fh->seek($final_seek_start, SEEK_SET) || die "seek failed";
  my $ignore = <$fh>;
  my @hits;
  my $blankify = $self->blankify();
  die "blankify unimplemented" if $blankify;
  my $hashify_headers = $self->hashify_headers();

  while (my $line = <$fh>) {
    next if $self->skip_comment_lines and $line =~ /^\#/;

    my %r;
    my $row = $self->parse_line($line);
    my $row_value = get_entry_value($row, $is_hash, \@cnames, \@cidx);

    eval $comparator_code;
    printf STDERR "final line read: value=%s comparison=%d\n", $row_value, $comparison if $DEBUG;

    if ($comparison == 0) {
      if ($hashify_headers) {
#	die sprintf "count mismatch: row:%d headers:%d\n%s\n",
#	scalar(@{$parsed}), scalar(@{$hashify_headers}), join "\n", @{$parsed}
#	unless @{$parsed} == @{$hashify_headers};
	# split issue / null trailing?
	my %r;
	@r{@{$hashify_headers}} = @{$row};
	push @hits, \%r;
      } else {
	push @hits, $row;
      }
    } elsif ($comparison == 1) {
      # no data, or passed site we want
      last;
    }
  }

  return \@hits;
}

sub seek_line {
  my ($self, $pos) = @_;
  my $fh = $self->fh();
  $fh->seek($pos, SEEK_SET) || confess "seek to $pos failed";
  my $ignore = <$fh>;
  # might have seeked to the middle of a row, so discard
  
  my $line_start_pos;
  my $line;
  while (1) {
    $line_start_pos = $fh->tell();
    $line = <$fh>;
    if ($self->skip_comment_lines) {
      last unless $line =~ /^\#/;
      # e.g. NHLBI
    } else {
      last;
    }
  }
  confess sprintf "line length=%d, configured max=%d, line=%s", length($line), $self->max_line_length(), $line if length($line) > $self->max_line_length();

  my $parsed = $self->parse_line($line);

  return ($parsed, $line_start_pos);
}

sub parse_line {
  my ($self, $line) = @_;
  chomp $line;
  $line =~ s/\r$//;

  if (my $c = $self->line_parser_callback) {
    # custom parser
    return &$c($line);
  } else {
    my $delim = $self->delimiter();
    unless ($delim) {
      $delim = $self->delimiter($self->detect_delimiter($line));
    }
    my @f = split /$delim/, $line;
    if ($self->is_headered) {
      if ($self->headers()) {
	my %r;
	@r{@{$self->headers}} = @f;
	return \%r;
      } else {
	$self->headers(\@f);
	return \@f;
	# for first line, just return headers
      }
    }
    return \@f;
  }
}


sub detect_delimiter {
  my ($self, $line) = @_;
  my %counts;
  foreach my $d (DELIMITERS) {
    my @f = split /$d/, $line;
    $counts{$d} = scalar @f;
  }
  my ($best) = sort {$counts{$b} <=> $counts{$a}} keys %counts;
  return $best;

}

sub get_entry_value {
  # STATIC
  my ($row, $is_hash, $cnames, $cidx) = @_;
  my $entry_value;
  if ($is_hash) {
    $entry_value = join ".", map {$row->{$_}} @{$cnames};
  } else {
    $entry_value = join ".", @{$row}[@{$cidx}];
  }
  return $entry_value;
}

sub get_fh_nearby {
  # return filehandle seeked nearby upstream of result
  my ($self) = @_;
  my $fh = $self->fh();
  my $final_seek_start = $self->final_seek_start();
  $fh->seek($final_seek_start, SEEK_SET) || die "seek failed";
  return $fh;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
