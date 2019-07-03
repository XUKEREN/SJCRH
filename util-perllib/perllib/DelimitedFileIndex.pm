package DelimitedFileIndex;
# describe me

use strict;
use Configurable;
use Exporter;

use FileHandle;
use POSIX qw(SEEK_SET);

@DelimitedFileIndex::ISA = qw(Configurable Exporter);
@DelimitedFileIndex::EXPORT_OK = qw();

use constant DELIMITERS => ("\t", ",");

use MiscUtils qw(log_message);
use FileUtils qw(newer_than);
use Reporter;
use DelimitedFile;

use MethodMaker qw(
file
	delimiter
column_number
index

ping
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->ping(250000);
  $self->configure(%options);
  return $self;
}

sub get_index {
  my ($self, %options) = @_;

  return $self->index() if $self->index();

  my $file = $self->file || die "-file";
  my $column_number = $self->column_number() || die "-column_number";
  # to do: option for label
  # 1st version is dbSNP (headerless)
  my $delimiter = $self->delimiter();
  die unless $column_number;
  my $column_index = $column_number - 1;
  my $file_pos;

  my $outfile = $self->get_outfile();

  my $needed;
  if (-s $outfile) {
    $needed = newer_than($file, $outfile);
  } else {
    $needed = 1;
  }

  if ($needed) {
    #
    # index build is required
    #
    my $current_v;
    log_message(sprintf "indexing %s, column %d.  This may take a while...", $file, $column_number);
    open(DFITMP, $file) || die "can't open $file";
    my %saw;
    my %pos;
    my $last_pos;
    my @ordered;
    my $count;
    my $ping = $self->ping();
    while (1) {
      $file_pos = tell DFITMP;
      my $line = <DFITMP>;
      last unless $line;

      chomp $line;

      $count++;
      if ($ping and $count % $ping == 0) {
	log_message(sprintf "indexing: records=%d, line=%s", $count, $line);
      }

      $delimiter = $self->detect_delimiter($line) unless $delimiter;
      my @f = split /$delimiter/, $line;
      my $v = $f[$column_index];

      if (not(defined $current_v) or $v ne $current_v) {
	if ($current_v) {
	  $pos{$current_v}{last_line_index} = $last_pos;
	  $pos{$current_v}{next_record_start} = $file_pos;
	}
	$pos{$v}{first_line_index} = $file_pos;
	$current_v = $v;
	die "out-of-order appearance of $v at $file_pos" if $saw{$v};
	$saw{$v} = 1;
	log_message(sprintf "indexing: found %s...", $v);

	push @ordered, $v;
      }
      $last_pos = $file_pos;
    }
    $pos{$current_v}{last_line_index} = $last_pos;
    $pos{$current_v}{next_record_start} = $file_pos;

    my $rpt = new Reporter(
      "-file" => $outfile,
      "-delimiter" => "\t",
      "-labels" => [
	qw(
   				            value
                                            first_line_index
                                            last_line_index
                                            next_record_start
					)
      ]
	);

    foreach my $v (@ordered) {
      my %r = %{$pos{$v}};
      $r{value} = $v;
      $rpt->end_row(\%r);
    }
    $rpt->finish();
  }

  #
  #  parse built index:
  #
  my $df = new DelimitedFile("-file" => $outfile,
			     "-headers" => 1,
			     );
  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  my %index;
  while (my $row = $df->get_hash()) {
    if (0) {
      foreach (sort keys %{$row}) {
	printf "%s: %s\n", $_, $row->{$_};
      }
    }
    $index{$row->{value}} = $row;
  }

  return $self->index(\%index);
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

sub get_outfile {
  my ($self) = @_;
  my $outfile = sprintf '%s.tab_index_%d',
  ($self->file || die),
  ($self->column_number || die);
  return $outfile;
}

sub test_index {
  my ($self, %options) = @_;
  my $index = $self->get_index();

  my $fh = new FileHandle();
  $fh->open($self->file || die);
  
  foreach my $f (sort keys %{$index}) {
    printf STDERR "%s:\n", $f;
    foreach my $key (qw(
first_line_index
last_line_index
)) {
      my $i = $index->{$f}{$key};
      $fh->seek($i, SEEK_SET) || die "seek failed";
      my $line = <$fh>;
      chomp $line;
      printf STDERR "  %s: %s\n", $key, $line;
    }
  }
}

sub find {
  my ($self, $value) = @_;
  my $index = $self->get_index;
  return $index->{$value};
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
