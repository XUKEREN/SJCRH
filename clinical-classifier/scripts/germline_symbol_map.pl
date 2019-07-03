#!/bin/env perl
# map germline reviewable/reportable symbols to different annotation
# namespaces (GENE_EXON_REGION or RefFlat)
# MNE 3/2016

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list get_hash_option);
use DelimitedFile;
use Reporter;
use GeneSymbolMapper;
use TdtConfig;
use FileUtils qw(write_simple_file);

my $F_GENE = "HGNC Name";
my $F_ALT = "Alternate Name Provided";

my %FLAGS;
my @clopts = (
	      "-genome=s",
	      "-reportable=s",
	      "-reviewable=s",
	      "-target=s",

	      "-refflat=s",
	      "-chr2gene=s",

	      "-tartan-import-all",
	      "-refflat-config",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my %expected_missing = map {$_, 1} qw(
DGCR
IGAD1
IGH
IGHA1
IGHA2
IGHE
IGHG1
IGHG2
IGHG3
IGHG4
IGHM
IGK
IGKC
IGL
SRY
TRA
TRAC
TRB
TRD
TRG
NUTM2B
CCDC26
				    );
# genomic/non-coding, e.g. only NG_/NR_ records

my $genome = $FLAGS{genome} || die "-genome";
my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";

if ($FLAGS{"tartan-import-all"}) {
  tartan_import();
  exit(0);
}

my $target = $FLAGS{target} || die "-target";
die "-target must be GENE_EXON_REGION|refflat" unless $target eq "GENE_EXON_REGION" or "refflat";
my $reportable = $FLAGS{reportable} || die "-reportable";
my $reviewable = $FLAGS{reviewable} || die "-reviewable";

my $rows_reportable = load_list($reportable);
my $rows_reviewable = load_list($reviewable);
gene_qc($rows_reportable, $rows_reviewable);

my $gsm = get_gsm();
if ($target eq "GENE_EXON_REGION") {
  my $ger_dir = $config_genome->{GENE_EXON_REGION_DIR} || die;
  printf STDERR "GENE_EXON_REGION: %s\n", $ger_dir;
  $gsm->populate_gene_exon_region("-dir" => $ger_dir);
} elsif ($target eq "refflat") {
  my $refflat;
  if ($FLAGS{"refflat-config"}) {
    $refflat = $config_genome->{REFSEQ_REFFLAT};
  } else {
    $refflat = $FLAGS{"refflat"} || die "-refflat FILE (must be the final version to be used by FusionBuilder)\n";
  }
  $gsm->populate_refflat("-refflat" => $refflat);
} elsif ($target eq "chr2gene") {
  my $chr2gene_dir = $FLAGS{chr2gene} || $config_genome->{CHR2GENE} || die;
  $gsm->populate_chr2gene("-dir" => $chr2gene_dir);
} else {
  die "target not implemented"
}

translate_list($gsm, $rows_reportable, "reportable");
translate_list($gsm, $rows_reviewable, "reviewable");

sub translate_list {
  my ($gsm, $rows, $tag) = @_;

  my @mapped;

  my $count_primary = 0;
  my $count_alt = 0;
  my $count_mapped = 0;
  my $count_missing = 0;
  my $count_missing_unexpected = 0;

  foreach my $row (@{$rows}) {
    my $gene = $row->{$F_GENE} || die "no gene";
    my $alt = $row->{$F_ALT};
    $alt = "" if $alt and $alt eq $gene;
    foreach ($gene, $alt) {
      s/^\s+//;
      s/\s+$//;
    }

    my $target_symbol;
    if ($gsm->contains($gene)) {
      $target_symbol = $gene;
      $count_primary++;
    } elsif ($alt and $gsm->contains($alt)) {
      printf STDERR "using alternative symbol %s for %s\n", $alt, $gene;
      $target_symbol = $alt;
      $count_alt++;
    } elsif (my $mapped = $gsm->resolve("-symbol" => $gene)) {
      printf STDERR "%s gene %s maps to %s\n", $tag, $gene, $mapped;
      $target_symbol = $mapped;
      $count_mapped++;
    }
    if ($target_symbol) {
      push @mapped, $target_symbol;
    } else {
      unless ($expected_missing{$gene}) {
	printf STDERR "ERROR: %s not found!\n", $gene;
	$count_missing_unexpected++;
      }
      $count_missing++;
    }
  }

  printf STDERR "genes:%d primary:%d alt:%d mapped:%d missing:%d unexpected_missing:%d\n\n",
    scalar (@{$rows}), $count_primary, $count_alt, $count_mapped, $count_missing, $count_missing_unexpected;

  my $outfile = get_outfile($tag, $target);

  write_simple_file(\@mapped, $outfile);

}

sub get_outfile {
  my ($tag, $target) = @_;
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  my $outfile = sprintf 'germline_%s_to_%s_%d_%02d_%02d.txt', $tag, $target, 1900 + $year, $mon + 1, $mday;
  return $outfile;
}


sub gene_qc {
  my ($reportable, $reviewable) = @_;

  my %reviewable = map {($_->{$F_GENE} || dump_die($_)), 1} @{$reviewable};
  my %reportable = map {($_->{$F_GENE} || dump_die($_)), 1} @{$reportable};

  my $errors;
  foreach my $g (sort keys %reportable) {
    unless ($reviewable{$g}) {
      printf STDERR "ERROR: reportable gene %s not in reviewable\n", $g;
      $errors = 1;
    }
  }
  die if $errors;

}

sub load_list {
  my ($infile) = @_;
  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  my @rows;
  while (my $row = $df->get_hash()) {
    dump_die($row, "no $F_GENE in $infile") unless $row->{$F_GENE};
    dump_die($row, "no $F_ALT in $infile") unless exists $row->{$F_ALT};

    push @rows, $row;
  }
  return \@rows;
}

sub get_gsm {
  my $config_species = TdtConfig::readConfig('species', "Homo_sapiens") || die "can't find hs species config";
  my $hgnc = $FLAGS{hgnc};
  my $gene_info = $FLAGS{"gene-info"};

  unless ($hgnc) {
    $hgnc = $config_species->{HGNC} || die "no HGNC config";
  }
  unless ($gene_info) {
    $gene_info = $config_species->{ENTREZ_GENE_GENEINFO} || die "no config ENTREZ_GENE_GENEINFO";
  }

  my $gsm = new GeneSymbolMapper(
				 "-hgnc_file" => $hgnc,
				 "-eg_file" => $gene_info,
				 "-hgnc_synonym_enable" => 0,
				 # disable for e.g. FAH/FANCA
				);

  return $gsm;
}

sub tartan_import {

  my %pair2param;
  $pair2param{GENE_EXON_REGION}{reviewable} = "CLINCLS_GERMLINE_REVIEWABLE_GENES_GER";
  $pair2param{GENE_EXON_REGION}{reviewable} = "CLINCLS_GERMLINE_REVIEWABLE_GENES_GER";



  foreach my $target (qw(GENE_EXON_REGION refflat chr2gene)) {
    foreach my $tag (qw(reviewable reportable)) {
      my $outfile = get_outfile($tag, $target);

      my $thing = $target;
      $thing = "GER" if $thing eq "GENE_EXON_REGION";

      my $param = sprintf 'CLINCLS_GERMLINE_%s_GENES_%s',
	uc($tag), uc($thing);

      my $cmd = sprintf 'tartan_import_helper.pl -genome %s -param %s -single-import -new-file %s', $genome, $param, $outfile;

      printf "%s\n", $cmd;

    }
  }

}
