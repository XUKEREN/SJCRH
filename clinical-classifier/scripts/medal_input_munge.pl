#!/bin/env perl
# extract columns required by medal ceremony input from other fields
# (for "close but no cigar" input files)

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

	      "-index2alleles",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $infile = $FLAGS{file} || die "-file";

my $df = new DelimitedFile("-file" => $infile,
			   "-headers" => 1,
			  );
my $outfile = basename($infile) . ".prep.tab";

my @extra;

if ($FLAGS{index2alleles}) {
  push @extra, qw(ReferenceAllele MutantAllele);
}

my $rpt = $df->get_reporter(
			    "-file" => $outfile,
			    "-extra" => \@extra,
			    "-auto_qc" => 1,
			   );

while (my $row = $df->get_hash()) {
  if ($FLAGS{index2alleles}) {
    my $idx = $row->{Index} || die "Index";
    my @f = split /\./, $idx;
    die unless @f == 5;
    $row->{ReferenceAllele} = $f[-2];
    $row->{MutantAllele} = $f[-1];
  }
  $rpt->end_row($row);
}

$rpt->finish();

