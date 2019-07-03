#!/bin/env perl
# sort medal ceremony output by medal (GSBClass)
# this is done by default, however parallelized output files (e.g. mux.pl)
# will lose this sorting

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list get_hash_option);
use FileUtils qw(count_file_lines);
use DelimitedFile;
use Reporter;

use constant CLASS_UNKNOWN => "Unknown";
use constant CLASS_GOLD => "Gold";
use constant CLASS_SILVER => "Silver";
use constant CLASS_BRONZE => "Bronze";

my %FLAGS;
my @clopts = (
	      "-in=s",
	      "-out=s",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $infile = get_hash_option(\%FLAGS, "in");
my $outfile = $FLAGS{out} || $infile . ".sorted.tab";

my %MEDAL_SORT_WEIGHT = (
			 CLASS_GOLD() => 300,
			 CLASS_SILVER() => 200,
			 CLASS_BRONZE() => 100,
			 CLASS_UNKNOWN() => 0.01,  # for hash presence check
			);

my $df = new DelimitedFile("-file" => $infile,
			   "-headers" => 1,
			  );
my $rpt = $df->get_reporter(
			    "-file" => $outfile,
			    "-auto_qc" => 1,
			   );

my @ordered = sort {$MEDAL_SORT_WEIGHT{$b} <=> $MEDAL_SORT_WEIGHT{$a}} keys %MEDAL_SORT_WEIGHT;

# while (my $row = $df->next("-ref" => 1)) {  # headerless
my %bucket;
while (my $row = $df->get_hash()) {
  my $medal = $row->{GSBClass} || dump_die($row, "no GSBClass");
  push @{$bucket{$medal}}, $row;
}

foreach my $medal (@ordered) {
  if (my $rows = $bucket{$medal}) {
    foreach my $row (@{$rows}) {
      $rpt->end_row($row);
    }
  }
}
$rpt->finish();

die "line count mismatch" unless count_file_lines($infile) == count_file_lines($outfile);


