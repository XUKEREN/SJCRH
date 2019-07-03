#!/bin/env perl
# some important SV genes are not present in refFlat database used
# by FusionBuilder, see:
#
# http://jira.stjude.org/browse/COMPBIO-2772
#
# create placeholder records in refFlat.txt format.
#
# MNE 11/2014

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use Reporter;
use DelimitedFile;

use TdtConfig;
use MiscUtils qw(dump_die);

my %FLAGS;

my @GENES;

GetOptions(\%FLAGS,
	   "-gene=s" => \@GENES,
	   # specify one or more gene symbols.  Must be as found in UCSC

	   "-genome=s",

	   "-hgnc=s",
	   # /nfs_exports/genomes/1/Homo_sapiens/HGNC/hgnc_complete_set.txt
	  );

my $genome = $FLAGS{genome} || die "-genome";
my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";

my $regions = $config_genome->{CLINCLS_CANCER_RELATED_GENES_FILE} || die;

my $hgnc = $FLAGS{hgnc};
unless ($hgnc) {
  my $config_species = TdtConfig::readConfig('species', "Homo_sapiens") || die "can't find hs species config";
  $hgnc = $config_species->{HGNC} || die "no HGNC config";
}
die "no HGNC" unless $hgnc;
printf STDERR "HGNC: %s\n", $hgnc;

die "specify -gene SYMBOL...\n" unless @GENES;

my %wanted_genes = map {$_, 1} @GENES;

#
#  find intervals for each gene:
#
my %gene2region;
open(IN, $regions) || die;
while (<IN>) {
  chomp;
  my ($gene, $interval) = split /\t/, $_;
  if ($wanted_genes{$gene}) {
    die "dup" if $gene2region{$gene};
    $gene2region{$gene} = $interval;
  }
}

#
#  find NG_ accessions for each gene:
#
my $df = new DelimitedFile("-file" => $hgnc,
			   "-headers" => 1,
			  );
# while (my $row = $df->next("-ref" => 1)) {  # headerless
my %gene2acc;
while (my $row = $df->get_hash()) {
  my $gene = $row->{"Approved Symbol"};
  $gene2acc{$gene} = $row->{"RefSeq IDs"} if $wanted_genes{$gene};
}

foreach my $gene (sort keys %wanted_genes) {
  my $interval = $gene2region{$gene} || die "no interval for $gene!";
  my $acc = $gene2acc{$gene} || die;
  $interval =~ s/\s+$//;
  $interval =~ /^(\d+):(\d+)\-(\d+)$/ || die "interval parse fail for $interval";
  my ($chrom, $start, $end) = ("chr" . $1, $2, $3);
  my $strand = "+";
  # always + since genomic sequence?
  # in any event transcript info is not available

  my @fields = (
		$gene,
		$acc,
		$chrom,
		$strand,
		$start,
		$end,
		$end,
		$end,
		"1",
		# exon count
		$start,
		$end,
		# exonStarts/exonEnds
	      );

  printf "%s\n", join "\t", @fields;
}
