#!/bin/env perl
# Medal Ceremony: Project Windows 95 ("it sucks less"):
# http://www.tomshardware.com/news/windows-95-vista-sucks-less-win95,11153.html

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use constant TYPE_SOMATIC => "somatic";
use constant TYPE_GERMLINE => "germline";

use MiscUtils qw(dump_die build_argv_list get_hash_option);
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
	      "-db-file=s",
	      "-db-name=s",

	      "-type=s",
	      # germline/somatic
	      "-file=s",
	      # IARC flatfile, after running hgvs_split.pl on it

	      "-test-load",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $f_db = $FLAGS{"db-file"} || die "-db-file";
my $dbi = DBI->connect("dbi:SQLite:dbname=$f_db","","");
my $type = $FLAGS{type} || die "-type";
my $infile = $FLAGS{file} || die "-file";

my %map;
$map{F_VDB_POS()} = "sj_pos";
$map{F_VDB_RA()} = "sj_ref_allele";
$map{F_VDB_VA()} = "sj_alt_allele";
# from hgvs_split.pl
$map{F_VDB_AA()} = "ProtDescription";
# from raw file

my $db_name;
if ($type eq TYPE_GERMLINE) {
  $db_name = "iarc_tp53_germline";

  foreach my $f (qw(
		     Generation
		     Individual_ID
		  )) {
    $map{$f} = $f;
    # passthrough, only in germline
  }
} elsif ($type eq TYPE_SOMATIC) {
  $db_name = "iarc_tp53_somatic";
  $map{F_VDB_PMID()} = "PubMed";
  # from raw file, only in somatic
} else {
  die "-type must be germline or somatic";
}

#
#  load database:
#
my $svdb = new SimpleVariantDB(
			       "-dbi" => $dbi,
			      );

if ($FLAGS{"test-load"}) {
  my $rows = $svdb->get_rows_with_constants(
					    "-name" => $db_name,
					   );
  dump_die($rows->[0]);
}

#
#  set constant fields applying to both databases:
#
my %const;
$const{F_VDB_GENE()} = "TP53";
$const{F_VDB_CHR()} = 17;
$const{F_VDB_TRANSCRIPT()} = "NM_000546.5";
# http://p53.iarc.fr/p53Sequences.aspx
# - page doesn't contain .5
# - current version at this time
$svdb->set_constants(
		     "-name" => $db_name,
		     "-values" => \%const
		    );

$svdb->add_blankify("p.?");
# set invalid protein specification to blank

$svdb->init_table(
		  "-name" => $db_name,
		  "-file" => $infile,
		  "-map" => \%map,
		  "-map-only" => 1,
		  "-auto-configure" => 1,
		  "-load" => 1,
		 );
