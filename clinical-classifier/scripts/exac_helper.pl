#!/bin/env perl
# exac utils (tartan only for now)

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;
use Cluster;
use CommandLineRebuilder;

use TARTANUtils qw(tartan_genome_put_helper tartan_clone_env);

use MiscUtils qw(dump_die build_argv_list);
# use DelimitedFile;
# use Reporter;
my $CLUSTER_RAM = 1000;

my %FLAGS;
GetOptions(
	   \%FLAGS,
	   "-tartan-index=s",
	   # output dir
	   "-genome=s",

	   "-tartan-clone",
	   "-from=s",
	   "-to=s",
	  );

my $genome = $FLAGS{genome} || die "-genome";
my $SUBDIR = "SUPPORT/ExAC";

if (my $out_dir = $FLAGS{"tartan-index"}) {
  tartan_genome_put_helper(
			   "-genome" => $genome,
			   "-subdir" => $SUBDIR,
			   "-out" => $out_dir,
			  );
} elsif ($FLAGS{"tartan-clone"}) {
  tartan_clone_env(
		   "-genome" => $genome,
		   "-subdir" => $SUBDIR,
		   "-from" => ($FLAGS{from} || die "-from"),
		   "-to" => ($FLAGS{to} || die "-to"),
		  );
} else {
  die "?";
}
