#!/bin/env perl
# lift an intervals file used by clinical classifier

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list get_hash_option);
use DelimitedFile;
use Reporter;
use LiftOver;

my %FLAGS;
my @clopts = (
	      "-file=s",
	      "-genome-from=s",
	      "-genome-to=s",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $infile = $FLAGS{file} || die;
my $outfile = $infile . ".liftover.tab";

my $genome_from = $FLAGS{"genome-from"} || die;
my $genome_to = $FLAGS{"genome-to"} || die;

my $lo = new LiftOver();
$lo->enable_retry(0) if $FLAGS{"no-retry"};

open(IN, $infile) || die;
my $wf = new WorkingFile($outfile);
my $fh = $wf->output_filehandle();

while (my $line = <IN>) {
  chomp $line;
  my @f = split /\t/, $line;
  die unless @f == 2;
  my ($gene, $interval) = @f;

  my @f2 = split /:/, $interval;
  die unless @f2 == 2;
  my ($chrom, $range) = @f2;

  my @f3 = split /\-/, $range;
  die unless @f3 == 2;
  my ($start, $end) = @f3;

  $lo->translate_base(
		      "-from" => $genome_from,
		      "-to" => $genome_to,
		      "-chr" => $chrom,
		      "-start-interbase" => $start - 1,
		      # HACK for deletions??
		      "-end-interbase" => $end
		     );

  if ($lo->translated_ok()) {
    my $new_chr = $lo->translated_chr();
    $new_chr =~ s/^chr//;
    my $new_start = $lo->translated_start_interbase() + 1;
    my $new_end = $lo->translated_end_interbase();
    printf $fh "%s\t%s:%d-%d\n", $gene, $new_chr, $new_start, $new_end;
  } else {
    die "fail translating $chrom $start $end";
  }

}

$wf->finish();
