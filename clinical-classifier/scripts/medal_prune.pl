#!/bin/env perl
# prune a medal results file to only those receiving a medal

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use File::Basename;
use Cluster;
use CommandLineRebuilder;

use MiscUtils qw(dump_die build_argv_list);
# use DelimitedFile;
# use Reporter;
use DelimitedFile;

my %FLAGS;
my @clopts = (
	   "-file=s",
	   "-files=s",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $files = build_argv_list("-flags" => \%FLAGS,
			    "-single" => "file",
			    "-set" => "files");

foreach my $infile (@{$files}) {
  my $outfile = $infile . ".pruned.tab";

  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			     );

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    my @medals = map {$row->{$_}} qw(Somatic_GSBClass GSBClass);
    my $usable = 0;
    foreach my $medal (@medals) {
      next unless $medal;
      next if lc($medal) eq "unknown";
      $usable = 1;
    }

    $rpt->end_row($row) if $usable;
  }
  $rpt->finish();

}
