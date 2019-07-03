#!/bin/env perl
# configuration repair script: do not use
# Michael Edmonson

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
	      "-dev=s",
	      "-prod=s",
	      "-variable=s",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $files = build_argv_list("-flags" => \%FLAGS,
			    "-single" => "file",
			    "-set" => "files");

my $infile = get_hash_option(\%FLAGS, "goodbad");
