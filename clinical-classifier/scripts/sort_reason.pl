#!/bin/env perl
# sort Reason field contents for easier comparison

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
  my $outfile = $infile . ".reason_sort.tab";

  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
  			      "-auto_qc" => 1,
			     );

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    foreach my $f (qw(Reason Somatic_Reason)) {
      if ($row->{$f}) {
	$row->{$f} = join ";", sort split /;/, $row->{$f};
      }
    }
    if (my $e = $row->{Evidence}) {
      my @things = split /:/, $e;
      my @out;
      foreach my $thing (@things) {
	my @f = split /=/, $thing;
	die unless @f == 2;
	if ($f[0] eq "COSMIC") {
	  my @c = split /,/, $f[1];
	  push @out, sprintf 'COSMIC=%s', join ",", sort @c;
#	  die join "\n", $thing, $out[$#out];
	} else {
	  push @out, $thing;
	}
      }
      $row->{Evidence} = join ":", @out;
    }

    $rpt->end_row($row);
  }
  $rpt->finish();

}


