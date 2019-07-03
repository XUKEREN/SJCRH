#!/bin/env perl
# generate promoter regions for germline-reviewable genes
# cast a wide net including 5' UTR plus a window upstream
# (e.g. canonical PTEN promoter significantly overlaps 5' UTR)
#
# MNE 4/2016
#
# TO DO:
# - manual interval override from dale

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;
use List::Util qw(min max);

use MiscUtils qw(dump_die build_argv_list get_hash_option);
use DelimitedFile;
use Reporter;
use SJPreferredIsoform;
use GeneSymbolMapper;
use RefFlatFile;
use Reporter;
use FileUtils qw(today_file);

use TdtConfig;

use constant GENOME_CONFIG_PARAM => "CLINCLS_GL_REVIEWABLE_PROMOTER_REGIONS";

#my $PROMOTER_FLANK = 5000;
# how much sequence upstream of the 5' UTR to use
my $PROMOTER_FLANK = 2500;
# 4/27/2016:
# Dale: [JZ] said she was good with 2.5k upstream and then up
# until the start of cds.

my %GENES_IGNORE = map {$_, 1} qw(
				   SRY
				);
# SRY is in PAR which is masked in GRCh37-lite and so filtered
# from our final refFlat

my %MANUAL_GENE_ACCESSION = (
			     "CCDC26" => "NR_130920",
);
# currently-installed HGNC is missing details for some genes.
# new version has this info, but *format has changed* so update
# will require rewrite.

my %CUSTOM_INTERVALS;
$CUSTOM_INTERVALS{APC} = {
			  "gene" => "APC",
			  "accession" => "NM_000038",
			  "chrom" => "chr5",
			  "start" => 112026616,
			  "end" => 112081649
			 };

my %FLAGS;
my @clopts = (
	      "-genome=s",
	      "-reviewable=s",
	      "-refflat=s",
	      "-gene=s",
	      "-tartan",
	      "-out=s",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $genome = $FLAGS{genome} || die "-genome";

if ($FLAGS{tartan}) {
  tartan_import();
  exit(0);
}


my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
#dump_die($config_genome);

my $f_refflat = $FLAGS{refflat} || $config_genome->{REFSEQ_REFFLAT} || die;

my $rff = new RefFlatFile();
$rff->strip_sharp_annotations(1);
$rff->canonical_references_only(1);
$rff->parse_file(
		 "-refflat" => $f_refflat,
		 "-type" => "refgene",
		);

my $f_pref = $config_genome->{GENE_TRANSCRIPT_MATRIX} || die;
printf STDERR "preferred: %s\n", $f_pref;
my $sjpi = new SJPreferredIsoform("-file" => $f_pref);
my $gsm = new_gsm();
my $pref_genes = $sjpi->get_genes;
foreach my $sym (@{$pref_genes}) {
  $gsm->add_gene("-gene" => $sym);
}

my $f_reviewable = $FLAGS{reviewable} || die "-reviewable";
my $df = new DelimitedFile(
			   "-file" => $f_reviewable,
			   "-headers" => 1,
			  );

my $gene = $FLAGS{gene};
my $outfile = get_outfile();

my $rpt = new Reporter(
		       "-file" => $outfile,
		       "-delimiter" => "\t",
		       "-labels" => [
				     qw(
					 gene
					 accession
					 chrom
					 start
					 end
				      )
#					 strand
#					 downstream_distance
				    ],
		       "-auto_qc" => 1,
		      );

while (my $row = $df->get_hash()) {
  my $symbols = get_ordered_symbols($row);

  my $ignore;
  foreach (@{$symbols}) {
    $ignore = 1 if $GENES_IGNORE{$_};
  }
  next if $ignore;

  if ($gene) {
    next unless grep {$_ eq $gene} @{$symbols};
  }

  #
  # find preferred isoform for each gene:
  #
  my $nm;
  my $method;

  my $idx = 0;
  foreach my $sym (@{$symbols}) {
    # direct lookup of primary and (sometimes) alternate symbol
    if ($nm = $sjpi->get_preferred_isoform($sym)) {
      $method = $idx == 0 ? "primary" : "alternate";
      last;
    }
    $idx++;
  }

  my $has_custom;
  foreach (@{$symbols}) {
    if (my $r = $CUSTOM_INTERVALS{$_}) {
      $rpt->end_row($r);
      $has_custom = 1;
      last;
    }
  }
  next if $has_custom;

  unless ($nm) {
    my $idx = 0;
    foreach my $sym (@{$symbols}) {
      if (my $mapped = $gsm->resolve("-symbol" => $sym)) {
	$method = $idx == 0 ? "gsm_primary" : "gsm_alternate";
	$nm = $sjpi->get_preferred_isoform($mapped) || die;
	last;
      }
      $idx++;
    }
  }

#  printf STDERR "ERROR: no NM for %s\n", $row->{gene} unless $nm;
  unless ($nm) {
    if ($nm = get_hgnc_nm($symbols->[0])) {
      $method = "hgnc";
    }
  }

  unless ($nm) {
    $method = "manual" if $nm = $MANUAL_GENE_ACCESSION{$symbols->[0]};
  }

  printf STDERR "%s: %s (%s)\n",
    $row->{gene}, ($nm || "MISSING"), ($method || "ERROR");

  # extract 5' UTR + upstream window:

  if ($nm) {
    my $hits = $rff->find_by_accession($nm);
    unless ($hits) {
      printf STDERR "ERROR: no refFlat entry for %s, ", $nm;
      my $nm_hgnc = get_hgnc_nm($symbols->[0]);
      $hits = $rff->find_by_accession($nm_hgnc) || die "hgnc lookup of $nm_hgnc failed";
      printf STDERR "using HGNC %s\n", $nm_hgnc;
    }
    foreach my $row (@{$hits}) {
      #
      # may be mappings to more than one location (e.g. C4A)
      #
      my $txStart = ($row->{txStart} || die) + 1;
      my $txEnd = $row->{txEnd} || die;
      my $cdsStart = ($row->{cdsStart} || die) + 1;
      my $cdsEnd = $row->{cdsEnd} || die;
      # +1: convert from interbase to in-base
      my $strand = $row->{strand} || die;
      my ($start, $end, $downstream_distance);

#      dump_die($row, join(" ", $row->{name}, $txStart, $txEnd, $cdsStart, $cdsEnd), 1);

      if ($strand eq "+") {
	$start = $txStart - $PROMOTER_FLANK;
	$end = $txStart + $PROMOTER_FLANK;
	my $max = $cdsStart - 1;
	$end = $max if $end > $max;
	$downstream_distance = $end - $txStart;
      } elsif ($strand eq "-") {
	$start = $txEnd - $PROMOTER_FLANK;
	$end = $txEnd + $PROMOTER_FLANK;
	my $min = $cdsEnd + 1;
	$start = $min if $start < $min;
	$downstream_distance = $start - $txEnd;
      } else {
	die;
      }

      my %r;
      $r{gene} = $symbols->[0];
      $r{accession} = $nm;
      $r{chrom} = $row->{chrom} || die;
      $r{start} = $start;
      $r{end} = $end;
      $r{strand} = $strand;
      $r{downstream_distance} = $downstream_distance;
      $rpt->end_row(\%r);
    }
  }
}
$rpt->finish();


sub get_ordered_symbols {
  my ($row) = @_;
  my @raw = map {$row->{$_} || die "no $_"} qw(gene alternate);
  my @symbols;
  my %saw;
  foreach my $sym (@raw) {
    push @symbols, $sym unless $saw{$sym};
    $saw{$sym} = 1;
  }
  return \@symbols;
}


sub new_gsm {
  my $config_species = TdtConfig::readConfig('species', "Homo_sapiens") || die "can't find hs species config";
  my $hgnc = $FLAGS{hgnc};
  my $gene_info = $FLAGS{"gene-info"};

  unless ($hgnc) {
    $hgnc = $config_species->{HGNC} || die "no HGNC config";
  }
  printf STDERR "HGNC: %s\n", $hgnc;
  unless ($gene_info) {
    $gene_info = $config_species->{ENTREZ_GENE_GENEINFO} || die "no config ENTREZ_GENE_GENEINFO";
  }

  my $gsm = new GeneSymbolMapper(
				 "-hgnc_file" => $hgnc,
				 "-eg_file" => $gene_info,
#				 "-hgnc_synonym_enable" => 0,
				 # disable for e.g. FAH/FANCA
				);

  return $gsm;
}

sub get_hgnc_nm {
  my ($sym) = @_;

  my $hgnc = $gsm->hgnc() || die;
  my $hgnc_hits = $hgnc->find("-symbol" => $sym, "-approved" => 1);
  die "need exactly 1 hgnc match" unless @{$hgnc_hits} == 1;

  die "no refseq field" unless exists $hgnc_hits->[0]->{"RefSeq (supplied by NCBI)"};
  my $nm_hgnc = $hgnc_hits->[0]->{"RefSeq (supplied by NCBI)"};
  $nm_hgnc = undef if $nm_hgnc =~ /NG_/;
  # ignore for accessions that won't be in refFlat

  return $nm_hgnc;

}

sub tartan_import {
  my $outfile = get_outfile();
  die "where is $outfile" unless -s $outfile;
  my $cmd = sprintf 'tartan_import_helper.pl -genome %s -single-import -new-file %s -param %s', $genome, $outfile, GENOME_CONFIG_PARAM;
  system $cmd;
}

sub get_outfile {
  my $outfile = $FLAGS{out} || today_file("reviewable_promoter_regions_%s.tab");
  return $outfile;
}
