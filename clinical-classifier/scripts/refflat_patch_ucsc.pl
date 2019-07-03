#!/bin/env perl
# some important SV genes are not present in refFlat database used
# by FusionBuilder, see:
#
# http://jira.stjude.org/browse/COMPBIO-2772
#
# fetch isoforms from UCSC and convert to refFlat.txt format.
#
# MNE 11/2014
#
# **** ABANDONED: now using entire NG_ region instead ****
#   => see refflat_patch_hgnc.pl


use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use Reporter;

use MiscUtils qw(dump_die);
# use DelimitedFile;
# use Reporter;

my %FLAGS;

my @GENES;

GetOptions(\%FLAGS,
	   "-gene=s" => \@GENES,
	   # specify one or more gene symbols.  Must be as found in UCSC

	   "-knownGene=s",
	   # /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/mRNA/UCSC/2014_11_19/knownGene.txt
	   # this should be a newer version for better chance of HGNC match

	   "-kgXref=s",
	   # /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/mRNA/UCSC/2014_11_19/kgXref.txt

	   "-debug",

	  );

die "specify -gene [symbol]...\n" unless @GENES;

my $kg = get_file_param("knownGene");
my $f_index = get_file_param("kgXref");

my %wanted_genes = map {$_, 1} @GENES;

open(IDX, $f_index) || die;
my @idx_fields = qw(
		     kgID
		     mRNA
		     spID
		     spDisplayID
		     geneSymbol
		     refseq
		     protAcc
		     description
		     rfamAcc
		     tRnaName
		  );
# see UCSC table browser

my %gene2index;
my %wanted_ids;

#
#  get index info for desired genes:
#
while (<IDX>) {
  chomp;
  my %r;
  my @f = split /\t/, $_, -1;
  die unless @f == @idx_fields;
  @r{@idx_fields} = @f;
  $r{raw_line} = $_;
  my $gene = $r{geneSymbol};
  if ($wanted_genes{$gene}) {
    push @{$gene2index{$gene}}, \%r;
    $wanted_ids{$r{kgID}} = 1;
  }
}
close IDX;

#
#  load knownGene rows for desired entries:
#
open(KG, $kg) || die;
while (<KG>) {
  chomp;
  my @f = split /\t/, $_, -1;
  my $id = $f[0];
  if ($wanted_ids{$id}) {
    $wanted_ids{$id} = \@f;
  }
}

#
#  process each gene:
#

my $debug_report = $FLAGS{debug};
my $rpt = new Reporter(
		       "-file" => "debug.tab",
		       "-delimiter" => "\t",
		       "-labels" => \@idx_fields,
		       "-auto_qc" => 1,
		      ) if $debug_report;

foreach my $gene (sort keys %wanted_genes) {
  if ($gene2index{$gene}) {
    my $idx_rows = $gene2index{$gene} || die "no index for $gene";

    my %saw;

    foreach my $idx_row (@{$idx_rows}) {
      $rpt->end_row($idx_row) if $debug_report;
      my $id = $idx_row->{kgID} || die;

      my $row = $wanted_ids{$id};
      die unless ref $row;

      my $thing = join "_", @{$row}[1..9];
      # check for pure-duplicate annotations

      my @f = ($gene, @{$row}[0 .. 9]);
      # change to refFlat format:
      #  - 10 columns
      #  - gene symbol first

      if (1) {
	# report NCBI mRNA accession rather than UCSC ID
	$f[1] = $idx_row->{mRNA} if $idx_row->{mRNA};
#	dump_die($idx_row);
      }




      printf "%s\n", join "\t", @f;
    }
  } else {
    printf STDERR "ERROR: no rows found for %s!\n", $gene;
  }
}
$rpt->finish() if $debug_report;


sub get_file_param {
  my ($p) = @_;
  my $file = $FLAGS{$p} || die "-" . $p;
  die "where is $file" unless -s $file;
  return $file;
}
