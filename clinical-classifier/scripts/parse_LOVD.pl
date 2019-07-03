#!/bin/env perl
# parse LOVD database HTML into a single, simpler tab-delimited file
# MNE 2/2014

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use File::Basename;

use Reporter;
use LOVDParser;

my %FLAGS;
GetOptions(\%FLAGS,
	   "-dir=s",
	   # directory of .htm files

	   "-limit=i",

	   "-out=s",
	  );

my $lp = new LOVDParser(
			"-field_reference" => "Variant_reference",
			"-field_aachange" => "Protein",
		       );

$lp->parse(
	   "-directory" => ($FLAGS{dir} || die "-dir"),
	   "-limit" => $FLAGS{limit}
	  );

my $outfile = $lp->export_to_file("-outfile" => $FLAGS{out});
