#!/bin/env perl
# just build tabix index for a bgzip-compressed file

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list get_hash_option);
use TabixFile;

my %FLAGS;
my @clopts = (
	      "-file=s",
	      "-f-chr=s",
	      "-f-start=s",
	      "-f-end=s",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $file = $FLAGS{file} || die "-file";
my $f_chr = $FLAGS{"f-chr"} || die "-f-chr";
my $f_start = $FLAGS{"f-start"} || die "-f-start";
my $f_end = $FLAGS{"f-end"};

my $tf = new TabixFile(
		       "-file" => $file,
		       "-f_chr" => $f_chr,
		       "-f_start" => $f_start,
		       "-f_end" => $f_end,
		       "-index" => 1,
		      );

