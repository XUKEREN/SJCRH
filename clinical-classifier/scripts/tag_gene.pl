#!/bin/env perl
# add gene annotations to flatfiles

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
use TdtConfig;
use GeneAnnotation;
use RefFlatFile;
use GeneSymbolMapper;

my $F_CHR;
my $F_START;
my $F_END;

my $OUT_F_GENES = "Genes";
my $STRIP_SHARP = 1;

my @GENE2INTERVAL;

my %FLAGS;
my @clopts = (
	      "-file=s",
	      "-files=s",
	      "-genome=s",
	      "-refflat=s",
	      "-conserting",
	      "-f-chr=s" => \$F_CHR,
	      "-f-start=s" => \$F_START,
	      "-f-end=s" => \$F_END,
	      "-strip-sharp=i" => \$STRIP_SHARP,
	      "-fq",

	      "-gene2interval=s" => \@GENE2INTERVAL,
	      "-verbose",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

if (@GENE2INTERVAL) {
  gene2interval();
  exit(0);
}

my $infiles = build_argv_list("-flags" => \%FLAGS,
			    "-single" => "file",
			    "-set" => "files");

my $f_refflat = get_refflat();

if ($FLAGS{conserting}) {
  $F_CHR = "chrom";
  $F_START = "loc.start";
  $F_END = "loc.end";
}

foreach my $f ($F_CHR, $F_START) {
  die "chr/pos not defined, missing field $f" unless defined $f;
}
$F_END = $F_START unless $F_END;

my $ga = new GeneAnnotation(
			    "-style" => "refgene_flatfile",
			    "-refgene_flatfile" => $f_refflat,
			    "-ignore_non_coding" => 0
			   );

my %uniq;
foreach my $infile (@{$infiles}) {
  my $outfile = ($FLAGS{fq} ? $infile : basename($infile)) . ".genes.tab";
  printf STDERR "in:%s out:%s\n", $infile, $outfile;
  die "duplicate outfile $outfile for $infile" if $uniq{$outfile};
  $uniq{$outfile} = $infile;

  if (-s $outfile) {
    printf STDERR "outfile %s exists\n", $outfile;
    next;
  }

  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			    );
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   $OUT_F_GENES
					  ],
			      "-auto_qc" => 1,
			     );

  while (my $row = $df->get_hash()) {


    foreach my $f ($F_CHR, $F_START, $F_END) {
      dump_die($row, "missing field $f") unless exists $row->{$f};
    }

    my $chrom = $row->{$F_CHR} || die;
    my $start = $row->{$F_START} || die;
    my $end = $row->{$F_END} || die;

    $ga->find(
	      "-reference" => $chrom,
	      "-start" => $start,
	      "-end" => $end
	     );
    my $genes_genomic = $ga->results_genes_genomic_order();

    my @genes;
    my %saw;
    foreach my $gene (@{$genes_genomic}) {
      next if $saw{$gene};

      if ($STRIP_SHARP and not($gene =~ /^_loc/)) {
	$gene =~ s/_locPar$//;
	$gene =~ s/_loc[A-Z]+$//;
      }

      push @genes, $gene;
      $saw{$gene} = 1;
    }
    $row->{$OUT_F_GENES} = join ",", @genes;

    $rpt->end_row($row);

  }
  $rpt->finish();
}

sub gene2interval {
  my $f_refflat = get_refflat();

  my $rf = new RefFlatFile();
  $rf->strip_sharp_annotations(1);
  $rf->parse_file(
		  "-refflat" => $f_refflat,
		  "-type" => "refgene",
		 );

  my @headers = qw(
		    gene_query
		    chr
		    start
		    end
		    gene_found
		    accession
		 );

  my $rpt = new Reporter(
			 "-fh" => \*main::STDOUT,
			 "-delimiter" => "\t",
			 "-labels" => \@headers,
			);

  # TO DO:
  # - preferred isoforms
  # - option to collapse multiple mappings if within XXX on the same chrom?

  my $gsm;

  my @expanded;
  foreach (@GENE2INTERVAL) {
    push @expanded, split /,/, $_;
  }

  foreach my $gene (@expanded) {
    my $hits = $rf->find_by_gene($gene);
    unless ($hits) {
      unless ($gsm) {
	$gsm = get_gsm();
	foreach my $row (@{$rf->rows}) {
	  $gsm->add_gene("-gene" => $row->{gene});
	}
      }
      if (my $alt = $gsm->resolve("-symbol" => $gene)) {
	$hits = $rf->find_by_gene($alt);
      }
    }
    if ($hits) {
      foreach my $hit (@{$hits}) {
	my %r;
	$r{gene_query} = $gene;
	$r{gene_found} = $hit->{gene} || die;
	$r{accession} = $hit->{name} || die;
	$r{chr} = $hit->{chrom} || die;
	$r{start} = $hit->{txStart} + 1;
	$r{end} = $hit->{txEnd};
	$rpt->end_row(\%r);
      }
    } else {
      my %r;
      foreach (@headers) {
	$r{$_} = "";
      }
      $r{gene_query} = $gene;
      $rpt->end_row(\%r);
    }
  }
  $rpt->finish();

}

sub get_refflat {
  my $f_refflat = $FLAGS{refflat};
  unless ($f_refflat) {
    my $genome = $FLAGS{genome} || die "-genome";
    my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
    $f_refflat = $config_genome->{REFSEQ_REFFLAT} || die "no REFSEQ_REFFLAT";
  }
  die "no refflat" unless $f_refflat;
  printf STDERR "refflat: %s\n", $f_refflat if $FLAGS{verbose};
  return $f_refflat;
}

sub get_gsm {
  my $species = "Homo_sapiens";
  my $config_species = TdtConfig::readConfig('species', $species) || die "can't find config or $species";

  my $gsm = new GeneSymbolMapper(
				 "hgnc_file" => ($config_species->{HGNC} || die),
				 "eg_file" => ($config_species->{ENTREZ_GENE_GENEINFO} || die),
				);
  return $gsm;
}
