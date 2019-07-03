#!/bin/env perl
use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;
use Cluster;
use CommandLineRebuilder;
use DelimitedFile;
use ChrBucketMap;
use Variant;
use Bambino2Annovar qw(genome2av get_annovar_db_dir);
use FileUtils qw(read_simple_file);
use TdtConfig;
use GeneSymbolMapper;

use LiftOverBatch (qw(
		       F_LO_IN_CHR
		       F_LO_IN_START
		       F_LO_IN_END
		       F_LO_OUT_CHR
		       F_LO_OUT_START
		       F_LO_OUT_END
		       F_LO_OUT_MINMATCH
		    ));

use GenomeUtils qw(reverse_complement);
use MiscUtils qw(dump_die build_argv_list);

my $QUEUE_SIZE = 1000;
my $ENABLE_RETRY = 1;
my $MINMATCH = 0.95;
# default
my $MINMATCH_FLOOR;
# if set, don't retry if minMatch falls below this threshold

my $PING = 0;

my %FLAGS;
my @clopts = (
	      "-file=s",
	      "-files=s",

	      "-genome-from=s",
	      "-genome-to=s",

	      "-fq",

	      "-queue=i" => \$QUEUE_SIZE,
	      "-retry=i" => \$ENABLE_RETRY,
	      "-ping=i" => \$PING,
	      "-minmatch=f" => \$MINMATCH,
	      "-minmatch-floor=f" => \$MINMATCH_FLOOR,

	      "-skip-invalid-positions",

	      "-out=s",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $files = build_argv_list("-flags" => \%FLAGS,
			    "-single" => "file",
			    "-set" => "files");

my $F_CHR = "chrom";
my $F_START = "loc.start";
my $F_END = "loc.end";

my $F_OUT_MINMATCH = "liftover_minmatch";

my $genome_from = $FLAGS{"genome-from"} || die "-genome-from";
my $genome_to = $FLAGS{"genome-to"} || die "-genome-to";

my ($total_lifted, $total_dropped);

foreach my $infile (@{$files}) {
  printf STDERR "processing %s...\n", $infile;
  $total_lifted = 0;
  $total_dropped = 0;
  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			    );

  my $outfile = $FLAGS{out};
  unless ($outfile) {
    $outfile = ($FLAGS{fq} ? $infile : basename($infile)) . ".liftover.tab";
  }

  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			     );

  my $headers = $rpt->labels();

  my @h_mod = ($F_CHR, $F_START, $F_END);
  my %h_mod = map {$_, 1} @h_mod;

  foreach my $h (@{$headers}) {
    $h = $h . ".orig" if $h_mod{$h};
  }
  push @{$headers}, @h_mod;
  # - rename columns where data might change, but preserve the same column number
  # - append new columns using the ordinary names for use w/annovar + medals
  push @{$headers}, $F_OUT_MINMATCH;

  my @queue;

  my $blo = new LiftOverBatch(
			      "-genome_from" => $genome_from,
			      "-genome_to" => $genome_to,
			      "-enable_retry" => $ENABLE_RETRY,
			      "-enable_strand_check" => 0,
			     );

  my $variant_count = 0;

  while (my $row = $df->get_hash()) {
    my ($chr, $start, $end) = @{$row}{$F_CHR, $F_START, $F_END};
    $variant_count++;
    if ($PING and $variant_count % $PING == 0) {
      printf STDERR "%s.%d: %d variants", $chr, $start, $variant_count;
    }

    foreach my $h (@h_mod) {
      $row->{$h . ".orig"} = $row->{$h};
    }

    die unless $chr and $start and $end;

    my $start_raw = $row->{$F_START} || die;
    my $end_raw = $row->{$F_END} || die;

    my $usable = 1;
    foreach ($start_raw, $end_raw) {
      unless (/^\d+$/) {
	my $msg = "non-integer position field detected";
	if ($FLAGS{"skip-invalid-positions"}) {
	  dump_die($row, $msg, 1);
	  $usable = 0;
	} else {
	  dump_die($row, $msg);
	}
      }
    }

    push @queue, $row if $usable;

    if (@queue >= $QUEUE_SIZE) {
      flush_queue("-queue" => \@queue, "-rpt" => $rpt, "-blo" => $blo);
      @queue = ();
    }

  }
  flush_queue("-queue" => \@queue, "-rpt" => $rpt, "-blo" => $blo);
  $rpt->finish();
}

sub flush_queue {
  my %options = @_;
  my $queue = $options{"-queue"} || die;
  my $rpt = $options{"-rpt"} || die;
  my $blo = $options{"-blo"} || die;

  #
  # assign interbase coordinates to each variant:
  #
  foreach my $row (@{$queue}) {
    my $start_raw = $row->{$F_START} || die;
    my $end_raw = $row->{$F_END} || die;

    my $i_start = $start_raw - 1;
    my $i_end = $end_raw;

    my $chr = $row->{$F_CHR} || die;
    $chr = "chr" . $chr unless $chr =~ /^chr/;

    $row->{F_LO_IN_CHR()} = $chr;
    $row->{F_LO_IN_START()} = $i_start;
    $row->{F_LO_IN_END()} = $i_end;
  }

  #
  # run batch liftover:
  #
  $blo->query_batch(
		    "-rows" => $queue,
		    "-min-match" => $MINMATCH,
		    "-min-match-floor" => $MINMATCH_FLOOR,
		   );

  #
  # translate results back into Bambino fields/coordinates:
  #
  foreach my $row (@{$queue}) {
    my $out_chr = $row->{F_LO_OUT_CHR()};
    my $out_start = $row->{F_LO_OUT_START()};
    my $out_end = $row->{F_LO_OUT_END()};
    my $out_minmatch = $row->{F_LO_OUT_MINMATCH()};
    $out_minmatch = sprintf '%.2f', $out_minmatch if $out_minmatch;

    if ($out_chr and $out_start and $out_end) {
      $row->{$F_CHR} = $out_chr;
      $row->{$F_START} = $out_start + 1;
      $row->{$F_END} = $out_end;
      # convert back from interbase to in-base
      $row->{$F_OUT_MINMATCH} = $out_minmatch;
      $rpt->end_row($row);
      $total_lifted++;
    } else {
      $total_dropped++;
    }
  }

  printf STDERR "%s: dropped:%d lifted:%d\n",
    scalar(localtime),
      $total_dropped, $total_lifted;
}
