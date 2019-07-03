#!/bin/env perl
# convert ClinVar VCF distribution to batch-tabix compatible
# e.g.
# ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20161101.vcf.gz

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use FileUtils qw(find_binary);
use MiscUtils qw(dump_die build_argv_list get_hash_option);
use DelimitedFile;
use Reporter;

my $F_CHR = "Chr";
my $F_POS = "WU_HG19_Pos";
my $F_RA = "ReferenceAllele";
my $F_VA = "MutantAllele";
my $F_CLNSIG = "CLNSIG";
my $F_DBSNP = "dbSNP";

my %FLAGS;
my @clopts = (
	      "-prep-vcf=s",
	      "-genome=s",

	      "-merge-citations=s",
	      "-cooked-tab=s",
	      "-stdout",

	      "-cooked-summary=s",

	      "-citations=s",
	     );

GetOptions(
	   \%FLAGS,
	   @clopts
	  );

if ($FLAGS{"prep-vcf"}) {
  prep_vcf();
} elsif ($FLAGS{"merge-citations"}) {
  merge_citations();
} elsif (my $cf = $FLAGS{"cooked-summary"}) {
  # summarize contents of final file
  cooked_summary($cf);
} else {
  die "?\n";
}

sub prep_vcf {
  my $f_vcf = $FLAGS{"prep-vcf"} || die;
  my $f_citations = $FLAGS{"citations"} || die "-citations";
  my $genome = $FLAGS{genome} || die "-genome";

  foreach my $bin (qw(vcf2tab.pl mux.pl bambino2annovar.pl annovar2medals.pl)) {
    find_binary($bin, "-die" => 1);
  }

  my $allow_cache = 1;

  my $cmd;
  my $f_step1 = basename($f_vcf) . ".citations.tab";
  my $f_step2 = $f_step1 . ".annovar_merged.tab";
  my $f_step3 = $f_step2 . ".mp.tab";
  my $f_step4 = $f_step3 . ".tabix.tab.gz";

  if ($allow_cache ? not(-s $f_step1) : 1) {
    my $tmp = $f_step1 . ".tmp";
    $cmd = sprintf 'vcf2tab.pl -file %s -clinvar 2 -lite -sj-post -no-insertion-adjustment -stdout | clinvar2tabix.pl -merge-citations %s -cooked-tab - -stdout > %s', $f_vcf, $f_citations, $tmp;
    system($cmd);
    die "error $?" if $?;
    die unless -s $tmp;
    rename($tmp, $f_step1) || die;
  }
  die unless -s $f_step1;

  if ($allow_cache ? not(-s $f_step2) : 1) {
    $cmd = sprintf 'mux.pl -file %s -ram 10000 -template "bambino2annovar.pl -genome %s -file %%s -drop-invalid-input" -suffix annovar_merged.tab -count 5000 -clean glob -wait 30', $f_step1, $genome;
    system $cmd;
  }
  die unless -s $f_step2;

  if ($allow_cache ? not(-s $f_step3) : 1) {
    $cmd = sprintf 'annovar2medals.pl -genome %s -file %s -keep-all -lite', $genome, $f_step2;
    system $cmd;
  }
  die unless -s $f_step3;

  if ($allow_cache ? not(-s $f_step4) : 1) {
    $cmd = sprintf 'vcf2tab.pl -tabix %s -sj-post', $f_step3;
    system $cmd;
  }
  die unless -s $f_step4;

}

sub merge_citations {
  my $infile = $FLAGS{"cooked-tab"} || die;
  my $f_citations = $FLAGS{"merge-citations"};
  # "var_citations.txt" from clinvar tab-delimited distributions, e.g.
  # ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt

  #
  # parse citations:
  #
  my $df = new DelimitedFile("-file" => $f_citations,
			     "-headers" => 1,
			     );
  my %citations;
  my %sources;
  while (my $row = $df->get_hash()) {
    my $cvid = $row->{'VariationID'} || dump_die($row, "where is VariationID");
    # switched from using dbSNP rs ID to ClinVar ID as it's unambiguous.
    # e.g. both 3.37059088.C.G and 3.37059088.C.T are annotated as
    # rs63751707 (dbSNP shows multiple alternate alleles), however
    # ClinVar IDs are unique (90406 and 90407)
    my $source = $row->{citation_source} || die;
    my $id = $row->{citation_id} || die;
#    printf STDERR "save $cvid $source $id\n";
    $citations{$cvid}{$source}{$id} = 1;
    $sources{$source} = 1;
  }

  #
  # merge clinvar output w/citations:
  #
  $df = new DelimitedFile(
			  "-file" => $infile,
			  "-headers" => 1,
			 );
  my $outfile = basename($infile) . ".citations.tab";
  my @ro;
  if ($FLAGS{stdout}) {
    @ro = ("-fh" => *STDOUT);
  } else {
    @ro = ("-file" => $outfile);
  }
  my $rpt = $df->get_reporter(
			      @ro,
			      "-extra" => [
					   sort keys %sources
					  ],
  			      "-auto_qc" => 1,
			     );

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
#    my $rs = $row->{variant_id} || die;
    # from vcf2tab, VCF distribution contains only records w/dbSNP
    # UGH, FORMAT CHANGE: previous ClinVar distribution used the VCF
    # ID column (column 3) to hold the dbSNP ID.  Later versions
    # changed this to the ClinVar ID.

    my $cvid = $row->{ClinVar_Variation_ID} || die "no ClinVar_Variation_ID";

    foreach my $source (keys %sources) {
      $row->{$source} = "";
    }

    if (my $c = $citations{$cvid}) {
      foreach my $source (sort keys %{$c}) {
	my @ids = sort keys %{$c->{$source}};
	$row->{$source} = join ",", @ids;
      }
    } else {
#      printf STDERR "no data for cvid %s\n", $cvid;
    }
    $rpt->end_row($row);
  }
  $rpt->finish();
}

sub cooked_summary {
  my ($fn) = @_;
  my $df = new DelimitedFile("-file" => $fn,
			     "-headers" => 1,
			     );
  my $variant_count = 0;
  my %pmid;
  my %descs;
  while (my $row = $df->get_hash()) {
    $variant_count++;
    if (my $cd = $row->{CLNSIG_desc}) {
      $descs{$cd}++;
    }

    if (my $pmids = $row->{PubMed}) {
      foreach (split /,/, $pmids) {
	$pmid{$_} = 1;
      }
    }
  }

  printf STDERR "variants:%d  publications:%d\n", $variant_count, scalar keys %pmid;
  printf STDERR "single 5-tier calls:\n";

  my @single = qw(
		   Benign
		   Likely_benign
		   Uncertain_significance
		   Likely_pathogenic
		   Pathogenic
		);
  foreach (@single) {
    printf STDERR "  %s: %d\n", $_, $descs{$_};
  }




}
