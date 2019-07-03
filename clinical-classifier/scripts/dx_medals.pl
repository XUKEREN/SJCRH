#!/bin/env perl
# run medal ceremony dnanexus
# MNE 2/2017

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
use DNANexusWrapper;

my $JOB_CHUNK_SIZE = 20000;
# this is the same value as used with mux.pl in the dnanexus bash run script.
# wish there was a way to put this value in one shared location.
my $DEFAULT_MIN_CORES = 4;
# since work is pretty intensive and process startup time is high,
# 2 cores seems too low
my $DEFAULT_MAX_CORES = 8;
# 16 and 32 seem significantly slower to start, and have lots of disk
# space we don't need

my %FLAGS;
my @clopts = (
	      "-file=s",

	      "-debug-no-wait",
	      "-no-auto-instance",
	      # if specified, default instance from dxapp.json
	      "-max-cores=i",
	      "-min-cores=i",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $infile = get_hash_option(\%FLAGS, "file");

my $dx = new DNANexusWrapper("-app" => 1);
# run as an app instead of an applet

unless ($FLAGS{"no-auto-instance"}) {
  # auto-select instance type (core count) based on # of records to process
  $dx->auto_instance($infile);
  $dx->job_chunk_size($JOB_CHUNK_SIZE);
  $dx->min_cores($FLAGS{"min-cores"} || $DEFAULT_MIN_CORES);
  $dx->max_cores($FLAGS{"max-cores"} || $DEFAULT_MAX_CORES);
}

$dx->wait_for_job(0) if $FLAGS{"debug-no-wait"};

$dx->json_setup("-fh" => \*main::DATA);

$dx->set_input(
	       "-parameter" => "infile",
	       "-value" => $infile,
	       "-clean" => 1
	      );

$dx->set_output_policy(
		       "-parameter" => "output_file",
		       "-policy" => {
				     "suffix" => "medals.tab.medals.tab"
				    }
		      );
$dx->run();

# JSON copied from application.  would be nice to get this from
# the filesystem but doubt these will be formally deployed anywhere,
# and casual users won't have a svn working copy.
__DATA__
{
  "name": "stjude_medal_ceremony",
  "title": "stjude_medal_ceremony",
  "summary": "stjude_medal_ceremony",
  "dxapi": "1.0.0",
  "version": "0.0.3",
  "inputSpec": [
    {
      "name": "infile",
      "label": "input file",
      "class": "file",
      "optional": false,
      "patterns": [
        "*"
      ],
      "help": ""
    }
  ],
  "outputSpec": [
    {
      "name": "output_file",
      "class": "file",
      "patterns": [
        "*"
      ],
      "help": ""
    }
  ],
  "timeoutPolicy": {
    "*": {
      "hours": 48
    }
  },
  "runSpec": {
    "interpreter": "bash",
    "file": "src/medal_ceremony.sh",
    "systemRequirements": {
      "*": {
        "instanceType": "mem1_ssd1_x4"
      }
    },
    "execDepends": [
                {
			"name": "libbio-samtools-perl"
                }, {
			"name": "vcftools"
		}, {
			"name": "tabix"
		}, {
			"name": "Data::Compare",
			"package_manager": "cpan",
			"version": "1.25"
		}, {
			"name": "HTML::TableExtract",
			"package_manager": "cpan",
			"version": "2.11"
		}, {
			"name": "libdbi-perl"

		}],

  "bundledDepends": [
    {
      "name": "bundle_prod_clinical-classifier_GRCh37-lite.tar.gz",
      "id": {
             "$dnanexus_link": "file-F1j0FyQ0qxY9GG65G542pvj3"
            }
    }
  ],

    "distribution": "Ubuntu",
    "release": "12.04"
  },
  "access": {
    "network": [
      "*"
    ]
  }
}
