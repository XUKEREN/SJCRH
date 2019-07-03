#!/bin/env perl
# /home/medmonso/work/jz/clinical_medals/databases/ALSoD
# Alsod_Mutation_Data.cleaned.txt converted from
# http://alsod.iop.kcl.ac.uk
# - make small repairs to spreadsheet first: patch in missing NM_, etc.

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);

use MiscUtils qw(dump_die build_argv_list get_hash_option);
use FileUtils qw(find_binary);
use DelimitedFile;
use Reporter;
use DBI;
use SimpleVariantDB qw(
			F_VDB_CHR
			F_VDB_POS
			F_VDB_RA
			F_VDB_VA
			F_VDB_GENE
			F_VDB_TRANSCRIPT
			F_VDB_AA
			F_VDB_PMID
		     );

my %FLAGS;
my @clopts = (
	      "-prep=s",
	      # parse/extract sub-columns before loading

	      "-out=s",

	      "-load=s",
	      # load generated file into database
	      "-db-file=s",
	     );

my @ARGV_RAW = @ARGV;
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

prep_data() if $FLAGS{prep};
load_db() if $FLAGS{load};

sub load_db {
  my $f_db = $FLAGS{"db-file"} || die "-db-file";
  my $dbi = DBI->connect("dbi:SQLite:dbname=$f_db","","");
  my $infile_raw = $FLAGS{load} || die "-load";
  my $infile = abs_path($infile_raw);

  my %map;
  $map{F_VDB_CHR()} = "sj_chr";
  $map{F_VDB_POS()} = "sj_pos";
  $map{F_VDB_RA()} = "sj_ref_allele";
  $map{F_VDB_VA()} = "sj_alt_allele";
  $map{F_VDB_GENE()} = "Gene";
  $map{F_VDB_AA()} = "sj_aachange";

  my $db_name = "ALSoD";

  my %constants;
  my %properties;
#  $properties{"command_line"} = join " ", $0, map {$_ eq $infile_raw ? $infile : $_} @ARGV_RAW;

  my $svdb = new SimpleVariantDB(
				 "-dbi" => $dbi,
				);

  $svdb->init_table(
		    "-name" => $db_name,
		    "-file" => $infile,
		    "-map" => \%map,
		    "-map-only" => 1,
		    "-auto-configure" => 1,
		    "-load" => 1,
		    "-properties" => \%properties,
		    "-constants" => \%constants,
		   );


}

sub prep_data {
  # step 1:
  # run hgvs_split.pl on Excel-covered file to extract
  # annotations in discrete columns
  # TO DO: embed source data in __DATA__?
  my $f_in = $FLAGS{"prep"} || die "-prep";

  find_binary("hgvs_split.pl", "-die" => 1);
  my $cmd = sprintf 'hgvs_split.pl -file %s -split-c -field-hgvs-c HGVS_Nucleotide -field-hgvs-p HGVS_protein -genome GRCh37-lite -field-chr-pos "Location(Chr)" -field-ra "Seq. Original" -field-va "Seq. Mutated"', $f_in;

  system $cmd;
  die if $?;
}

