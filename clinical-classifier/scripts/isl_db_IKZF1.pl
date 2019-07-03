#!/bin/env perl
# extract variant data from cut-and-pasta figure 1 table data from
# https://www.ncbi.nlm.nih.gov/pubmed/29681510

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);

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

use constant PMID => 29681510;
use constant GENE_NAME => "IKZF1";

my %FLAGS;
my @clopts = (
	      "-generate",
	      # extract data and generate flatfile
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

generate_flatfile() if $FLAGS{generate};
load_db() if $FLAGS{load};

sub load_db {
  my $f_db = $FLAGS{"db-file"} || die "-db-file";
  my $dbi = DBI->connect("dbi:SQLite:dbname=$f_db","","");
  my $infile_raw = $FLAGS{load} || die "-load";
  my $infile = abs_path($infile_raw);

  my %map;
  $map{F_VDB_POS()} = "WU_HG19_Pos";
  $map{F_VDB_RA()} = "ReferenceAllele";
  $map{F_VDB_VA()} = "MutantAllele";
  $map{F_VDB_AA()} = "AAChange";

  foreach my $f (qw(
		   number_of_cases_identified
		 functionally_damaging
		  )) {
    # passthrough
    $map{$f} = $f;
  }

  my %constants;
  $constants{F_VDB_CHR()} = 7;
  $constants{F_VDB_GENE()} = GENE_NAME();
  $constants{F_VDB_PMID()} = PMID();

  my %properties;

  my $svdb = new SimpleVariantDB(
				 "-dbi" => $dbi,
				);

  $svdb->init_table(
		    "-name" => GENE_NAME(),
		    "-file" => $infile,
		    "-map" => \%map,
		    "-map-only" => 1,
		    "-auto-configure" => 1,
		    "-load" => 1,
		    "-properties" => \%properties,
		    "-constants" => \%constants,
		   );


}

sub generate_flatfile {

  my $f_out = $FLAGS{out} || "IKZF1_pmid_29681510_table_1.tab";

  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   GeneName
					   Chr
					   pos_raw
					   ra_raw
					   ma_raw
					   AAChange
					   WU_HG19_Pos
					   ReferenceAllele
					   MutantAllele
					   indel_adjust_type
					   number_of_cases_identified
					   functionally_damaging
					)
				      ],
			 "-auto_qc" => 1,
			);

  while (my $line = <DATA>) {
    chomp $line;
    #  print STDERR "line=$line\n";
    my @f = split /\s+/, $line;
    my %r;
    $r{GeneName} = "IKZF1";
    $r{Chr} = 7;
    $r{AAChange} = $f[0];

    my ($pos, $ra, $ma) = @f[2,3,4];
    die "$pos $ra $ma" unless defined($pos) and defined($ra) and defined($ma);
    my $indel_adj_type = "";
    my $ra_raw = $ra;
    my $ma_raw = $ma;
    my $pos_raw = $pos;
    if (length($ra) != length($ma)) {
      # this is NOT VCF allele style!

      if (length($ra) == 2 and length($ma) == 1 and
	  substr($ra, -1) eq $ma) {
	# single-base deletion at target base, e.g. AG => G
	$ra = substr($ra, 0, 1);
	$ma = "-";
	$indel_adj_type = "mickey_mouse_this";
      } elsif (length($ra) == 2 and length($ma) == 1 and
	       substr($ra, 0, 1) eq $ma) {
	# single-base deletion of NEXT base, e.g. GA => G
	$ra = substr($ra, -1);
	$ma = "-";
	$pos++;
	#      die "$pos_raw $ra_raw $ma_raw $ra $ma $pos";
	$indel_adj_type = "mickey_mouse_next";
	# "What is this Mickey Mouse sh*t?" -Gunnery Sergeant Hartman
      } else {
	die "unhandled $pos $ra $ma";
      }
    }

    $r{ra_raw} = $ra_raw;
    $r{ma_raw} = $ma_raw;
    $r{pos_raw} = $pos_raw;

    $r{WU_HG19_Pos} = $pos;
    $r{ReferenceAllele} = $ra;
    $r{MutantAllele} = $ma;
    $r{number_of_cases_identified} = $f[6];
    $r{functionally_damaging} = $f[$#f];
    $r{indel_adjust_type} = $indel_adj_type;

    # only load records marked functionally damaging:
    my $v_damaging = $r{functionally_damaging};
    my $is_damaging;
    if ($v_damaging eq "yes") {
      $is_damaging = 1;
    } elsif ($v_damaging eq "no") {
      # ignore
    } else {
      die "unhandled damaging value $v_damaging";
    }

    $rpt->end_row(\%r) if $is_damaging;
  }
  $rpt->finish();
}


__DATA__
Pro18Thr (P18T) 50367245 C A missense 1 0 0 21.2 damaging tolerated no
Met31Val (M31V) 50367284 A G missense 2 3.1 3 105 3.2 3 105 18.7 damaging tolerated yes
Val52Leu (V52L) 50367347 G C missense 1 0 0 5.8 benign tolerated yes
Val53Met (V53M) 50367350 G A missense 4 3.0 3 104 8.0 3 104 9.1 probably damaging tolerated yes
Arg69His (R69H) 50444276 G A missense 1 2.0 3 104 6.5 3 105 15.7 benign tolerated no
Asp81Asn (D81N) 50444311 G A missense 1 9.5 3 106 3.2 3 105 28.4 damaging tolerated yes
Ser105Leu (S105L) 50444384 C T missense 1 9.5 3 106 0 16.5 probably damaging tolerated yes
Arg162Pro (R162P) 50450301 G C missense 1 0 0 27.6 damaging deleterious yes
His163Tyr (H183Y) 50450303 C T missense 1 0 0 27.5 damaging deleterious yes
Asp186fs (D186fs) 50450369 AG G frameshift 2 0 0 36.0 NA NA yes
Asp252Asn (D252N) 50459465 G A missense 3 1.0 3 104 3.2 3 105 14.8 damaging tolerated yes
Ser258Pro (S258P) 50459483 T C missense 1 0 0 15.1 benign tolerated yes
Met306* (M306*) 50467680 GA G nonsense 1 0 0 36.0 NA NA yes
Thr333Ala (T333A) 50467762 A G missense 1 1.0 3 104 1.0 3 104 17.4 damaging tolerated yes
Gly337Ser (G337S) 50467774 G A missense 10 1.2 3 103 3.1 3 103 3.3 benign tolerated yes
Met347Val (M347V) 50467804 A G missense 1 0 0 1.1 benign tolerated yes
Tyr348Cys (Y348C) 50467808 A G missense 2 0 0 18.3 damaging deleterious yes
Ala365Val (A365V) 50467859 C T missense 1 0 0 17.2 damaging tolerated yes
Cys394* (C394*) 50467947 C A nonsense 1 0 0 35.0 NA tolerated yes
Leu411Phe (L411P) 50467996 C T missense 1 0 0 19.7 damaging tolerated no
Pro420Gln (P420Q) 50468024 C A missense 1 0 0 13.2 benign tolerated no
Arg423Cys (R423C) 50468032 C T missense 1 0 0 10.1 damaging deleterious yes
His432Gln (H423Q) 50468061 C G missense 1 0 0 0.7 benign tolerated no
Ala434Gly (A434G) 50468066 C G missense 1 0 0 16.4 probably damaging deleterious yes
Leu449Phe (L449F) 50468110 C T missense 1 0 0 0.1 benign tolerated yes
Met459Val (M459V) 50468140 A G missense 1 0 0 0.2 benign tolerated yes
Met476Thr (M476T) 50468192 T C missense 1 0 0 22.7 damaging tolerated yes
Met518Lys (M518K) 50468318 T A missense 1 1.1 3 105 0 15.6 probably damaging deleterious no
