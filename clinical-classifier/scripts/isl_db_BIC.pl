#!/bin/env perl
# Medal Ceremony: Project Windows 95 ("it sucks less"):
# http://www.tomshardware.com/news/windows-95-vista-sucks-less-win95,11153.html
# loader for NHGRI BIC (BRCA1/BRCA2)

use strict;
use warnings;
use Cwd qw(abs_path);

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use constant NHGRI_BRCA1_NM => "NM_007294.3";
use constant NHGRI_BRCA2_NM => "NM_000059.3";
# TO DO: REFERENCE?

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

	      "-file=s",
	      # BIC flatfile after running hgvs_split.pl

	      "-test-load",
	     );
my @ARGV_RAW = @ARGV;
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $f_db = $FLAGS{"db-file"} || die "-db-file";
my $dbi = DBI->connect("dbi:SQLite:dbname=$f_db","","");
my $infile_raw = $FLAGS{file} || die "-file";
my $infile = abs_path($infile_raw);

my %map;
$map{F_VDB_POS()} = "sj_pos";
$map{F_VDB_RA()} = "sj_ref_allele";
$map{F_VDB_VA()} = "sj_alt_allele";
# from hgvs_split.pl
$map{F_VDB_AA()} = "HGVS Protein";
# from raw file

foreach my $raw (
	       "Clinically Important",
	       "Depositor",
	       "Ethnicity",
	      ) {
  my $cooked = $raw;
  $cooked =~ s/\s+/_/g;
  $map{$cooked} = $raw;
  # passthrough
}

my $db_name;
my ($gene, $chrom, $nm);
if ($infile =~ /(brca[12])/) {
  $gene = uc($1);
  if ($gene eq "BRCA1") {
    $chrom = 17;
    $nm = NHGRI_BRCA1_NM;
  } elsif ($gene eq "BRCA2") {
    $chrom = 13;
    $nm = NHGRI_BRCA2_NM;
  }
  $db_name = "NHGRI_BIC_" . uc($1);
} else {
  die "$infile must contain brca1 or brca2";
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
  foreach my $r (@{$rows}) {
    #    dump_die($rows->[0]);
    printf "%s\n", join " ", map {">" . $_ . "<"} $r->{svdb_va};
  }
  exit(0);
}

my %properties;
$properties{"command_line"} = join " ", $0, map {$_ eq $infile_raw ? $infile : $_} @ARGV_RAW;
#dump_die(\%properties);
# hack but close enough
# TO DO: version

#
#  set constant fields applying to both databases:
#
my %constants;
$constants{F_VDB_GENE()} = $gene;
$constants{F_VDB_TRANSCRIPT()} = $nm;
$constants{F_VDB_CHR()} = $chrom;
# FIX ME: "p," entry

$svdb->add_blankify("p.?");
# set invalid protein specification to blank

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
