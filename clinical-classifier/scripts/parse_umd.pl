#!/bin/env perl
# parse UMD database flatfiles into a single, simpler tab-delimited file
# MNE 2/2014

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use File::Basename;

use HTMLUtils qw(parse_html_tables);
use Reporter;

use PostprocessedVariantReport qw(
				   F_GENE
				   F_REFSEQ
				   F_AACHANGE
);

my %FLAGS;
GetOptions(\%FLAGS,
	   "-dir=s",
	   # directory of .htm files
	   # e.g. /nfs_exports/genomes/1/projects/ClinicalSeq/germline/umd.be/full_2014_02_21

	   "-out=s",
	  );

my $outfile = $FLAGS{out} || "umd_combined.tab";

my @files = glob(sprintf "%s/*.htm", $FLAGS{dir} || die "-dir");

my %gene2nm;
my @variant_files;
foreach my $fn (@files) {
  if ($fn =~ /gene_info/) {
    # map gene name to NM_
    parse_gene_info($fn, \%gene2nm);
  } else {
    push @variant_files, $fn;
  }
}

my $total = 0;

my @fields = (
	      "Old nomenclature",
	      "cDNA Nomenclature",
	      "exon",
	      "mRNA nomenclature",
	      "Protein nomenclature",
	      "Functional Domain",
	      "Biological significance",
	      "Validation by",
	      "Validation date",
	      "# records",
	     );


my $rpt = new Reporter(
		       "-file" => $outfile,
		       "-delimiter" => "\t",
		       "-labels" => [
				     F_GENE,
				     F_REFSEQ,
				     F_AACHANGE,
				     @fields,
				     "Biological_significance_cleaned",
				     "file"
				    ],
			 "-auto_qc" => 1,
		      );

my %SIGNIFICANCE_MAP = (
			"1 - Neutral" => "B",
			"Neutral" => "B",

			"2 - Likely Neutral" => "LB",
			"2 - Likely neutral" => "LB",
			"Likely Neutral" => "LB",

			"3 - UV" => "U",
			"UV" => "U",
			"Polymorphism" => "U",
			# ?

			"4 - Likely causal" => "LP",
			"Likely Causal" => "LP",

			"5 - Causal" => "P",
			"Causal" => "P",

		       );
# standardize observed codes to SJ 5-tier:
# B  = benign
# LB = likely benign
# U  = unknown
# LP = likely pathological
# P  = pathological

my $newlines_repaired = 0;

foreach my $vf (@variant_files) {
  printf STDERR "%s...\n", $vf;
  my $gene = get_gene_name($vf);
  my $tables = parse_html_tables(
				 "-file" => $vf,
				 "-require-header" => "Protein nomenclature",
				 # "-dump" => 1,
				);
  die "can't find (single) result table" unless @{$tables} == 1;

  foreach my $row (@{$tables->[0]}) {
    foreach my $f (@fields) {
      die unless exists $row->{$f};
      if (my $v = $row->{$f}) {
	if ($v =~ /\n/) {
	  $v =~ s/\n/ /g;
	  $newlines_repaired++;
	  $row->{$f} = $v;
	}
      }
    }

    $row->{"# records"} =~ s/\W//g;
    # remove lingering control characters

    my $cleaned = $row->{"Biological significance"};
    $cleaned = "" unless defined $cleaned;
    if ($cleaned) {
      $cleaned = $SIGNIFICANCE_MAP{$cleaned} || die "can't clean significance $cleaned";
    }
    $row->{Biological_significance_cleaned} = $cleaned;

    $row->{file} = basename($vf);
    $row->{F_GENE()} = $gene;
    $row->{F_REFSEQ()} = $gene2nm{$gene} || die;
    $row->{F_AACHANGE()} = $row->{"Protein nomenclature"} || "";

    $rpt->end_row($row);
  }
}
$rpt->finish();
printf STDERR "newlines repaired: %d\n", $newlines_repaired;

sub parse_gene_info {
  my ($fn, $gene2nm) = @_;
  open(IN, $fn) || die;
  my %hits;
  while (<IN>) {
    if (/(NM_\d+\.\d+)/) {
      $hits{$1} = 1;
    } elsif (/(NM) (_\d+\.\d+)/) {
      my $nm = $1 . $2;
      $hits{$nm} = 1;
    }
  }
  my @hits = keys %hits;
  die "can't find NM in $fn" unless @hits == 1;

  my $gene = get_gene_name($fn);
  $gene2nm->{$gene} = $hits[0];
}


sub get_gene_name {
  my ($fn) = @_;
  basename($fn) =~ /^(\w+?)_/ || die;
  return $1;
}
