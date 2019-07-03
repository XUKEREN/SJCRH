#!/bin/env perl
# Medal Ceremony: Project Windows 95 ("it sucks less"):
# http://www.tomshardware.com/news/windows-95-vista-sucks-less-win95,11153.html
# loader for ASU TERT
# http://telomerase.asu.edu/diseases.html#tert
#
# TO DO:
# get genomic positions via mutalyzer?

use strict;
use warnings;
use Cwd qw(abs_path);

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list get_hash_option);
use FileUtils qw(read_simple_file);
use DelimitedFile;
use Reporter;
use DBI;
use HTMLUtils qw(parse_html_tables);
use SimpleVariantDB qw(
			F_VDB_GENE
			F_VDB_TRANSCRIPT
			F_VDB_AA
			F_VDB_CDS
			F_VDB_PMID
		     );

my %FLAGS;
my @clopts = (
	      "-db-file=s",

	      "-file=s",
	      # copy of HTML from
	      # http://telomerase.asu.edu/diseases.html

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
foreach my $f (
	       F_VDB_GENE,
	       F_VDB_TRANSCRIPT,
	       F_VDB_AA,
	       F_VDB_CDS,
	       F_VDB_PMID
	      ) {
  $map{$f} = $f;
  # passthrough
}

my $db_name = "ASU_TERT";

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
  my @h = sort keys %{$rows->[0]};
    printf "%s\n", join " / ", @h;
  foreach my $r (@{$rows}) {
#    printf "%s\n", join " ", map {">" . $_ . "<"} $r->{svdb_aachange};
    printf "%s\n", join " / ", map {$r->{$_}} @h;
  }
  exit(0);
}

my %properties;
$properties{"command_line"} = join " ", $0, map {$_ eq $infile_raw ? $infile : $_} @ARGV_RAW;
#dump_die(\%properties);
# hack but close enough
# TO DO: version

my ($rows, $nm) = parse_asu_tert($infile);

#
#  set constant fields:
#
my %constants;
$constants{F_VDB_GENE()} = "TERT";
$constants{F_VDB_TRANSCRIPT()} = $nm;

$svdb->add_blankify("p.(?)");
$svdb->add_blankify("c.(?)");
$svdb->add_blankify("n/a");

$svdb->init_table(
		  "-name" => $db_name,
		  "-rows" => $rows,
		  "-auto-configure" => 1,
		  "-load" => 1,
		  "-properties" => \%properties,
		  "-constants" => \%constants,
		 );

sub parse_asu_tert {
  #
  #  parse ASU data for TERT:
  #  http://telomerase.asu.edu/diseases.html
  #
  my ($asu_tert) = @_;

  local $/ = undef;
  open(IN, $asu_tert) || die;
  my $blob_all = <IN>;
  close IN;
  my @chunks = split /a name="/, $blob_all;
  my @hits = grep {/^tert/i} @chunks;
  die unless @hits == 1;
  my $blob = $hits[0];
  # extract just the TERT section

  @hits = $blob =~ /(NM_\d+)/g;
  die unless @hits == 1;
  my $nm = $hits[0];
  # isoform

#  my $tables = parse_html_tables("-file" => $asu_tert, "-keep-html" => 1);
  my $tables = parse_html_tables("-blob" => $blob, "-keep-html" => 1);
  die unless @{$tables} == 1;
  my $rows = $tables->[0];
  # variant rows

  my $GENE = "TERT";
  # hack

  my $last_domain;
  my @out_rows;

  foreach my $row (@{$rows}) {
    if (my $d = $row->{Domains}) {
      $last_domain = $d;
    } else {
      $row->{Domains} = $last_domain;
    }

    my $aa = $row->{"AA substitution"};
    if ($aa and $aa =~ /</) {
      # hack: strip HTML from AA (still needed for References to parse info)
      $aa =~ s/<[^>]+>//g;
      $row->{"AA substitution"} = $aa;
    }

    my $cds = $row->{Mutation} || "";
    if ($cds and $cds =~ /</) {
      $cds =~ s/<[^>]+>//g;
    }

    my @pmid = $row->{References} =~ /pubmed\/(\d+)/g;

    my %r;
    $r{F_VDB_PMID()} = join ",", @pmid;
    $r{F_VDB_AA()} = $aa || "";
    $r{F_VDB_CDS()} = $cds;
    $r{Presentation} = $row->{Presentation} || "";
    push @out_rows, \%r;
  }

  return (\@out_rows, $nm);
}
