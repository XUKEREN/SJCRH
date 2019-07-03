#!/bin/env perl
# build intervals list for oncogenes
# MNE 3/2019

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
use TSOncoDB;
use TdtConfig;
use FileUtils qw(read_simple_file write_simple_file);
use GeneSymbolMapper qw(new_gsm_lite);
use RefFlatFile;
use Reporter;

my $BUFFER_UPSTREAM = 1000;

my %FLAGS;
my @clopts = (
	      "-genome=s",
	      "-build-raw-gene-list",
	      "-build-clean-gene-list=s",

	      "-build-intervals=s",
	      "-out=s",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

build_raw_gene_list() if $FLAGS{"build-raw-gene-list"};
build_clean_gene_list($FLAGS{"build-clean-gene-list"}) if $FLAGS{"build-clean-gene-list"};
build_intervals() if $FLAGS{"build-intervals"};

sub build_raw_gene_list {
  # will likely only work for GRCh37-lite as medal ceremony isn't
  # formally supported on 38 yet
  #
  # from email to Steve Rice:
  #
  # 1) SJ “recurrent” genes from JZ’s COSMIC analysis which have
  # historically been semi-synonymous with oncogenes (but not always,
  # recurrence does not necessarily mean gain of function, sometimes
  # it’s loss). 2) amplification genes from medal ceremony somatic CNV
  # config, which are presumably gain-of-function 3) genes from medal
  # ceremony SV configuration (presumably oncogenic fusions) 4) named
  # oncogenes from OncoKB database

  my %all;

  # any annotation database source, which covers items 1 and 4 above:
  my $tsdb = new TSOncoDB();
  foreach my $g (@{$tsdb->get_oncogenes()}) {
    $all{$g} = 1;
  }

  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";

  # item 2 above:
  my $f_somatic_cnv = $config_genome->{"CLINCLS_CNV_CHECK_FILE"} || die;
  open(TMP, $f_somatic_cnv) || die;
  while (<TMP>) {
    my @f = split /\t/, $_;
    my ($g, $type) = @f;
    $all{$g} = 1 if $type eq "Amp";
  }

  foreach my $var (qw(
		       CLINCLS_SV_CHECK_FILE
		       CLINCLS_SV_CHECK_MANUAL_FILE
		    )) {
    my $f = $config_genome->{$var} || die;
    open(TMP, $f) || die;
    while(<TMP>) {
      my ($g1, $g2, @stuff) = split /\t/, $_;
      foreach my $g ($g1, $g2) {
	next if $g eq "pair1" or $g eq "pair2";
	# headers in 2nd file
	$all{$g} = 1 if $g;
	# might be single or a pair
      }
    }
  }

  write_simple_file([ sort keys %all ], "genes_raw.txt");
}

sub build_clean_gene_list {
  my ($infile) = @_;
  my $genes_raw = read_simple_file($infile);

  my %reject;
  foreach my $g (@{$genes_raw}) {
    my @others = grep {$_ ne $g} @{$genes_raw};
    my $gsm = new_gsm_lite();
    foreach my $o (@others) {
      $gsm->add_gene("-gene" => $o);
    }
    if (my $hit = $gsm->find($g)) {
      # ambiguity
      my @two = ($g, $hit);
      my $hgnc = $gsm->hgnc();
      my @approved;
      foreach my $check (@two) {
	push @approved, $check if $hgnc->find("-symbol" => $check,
					      "-approved" => 1);
      }
      die unless @approved == 1;
      # only 1 should be an approved symbol

      my ($reject) = grep {$_ ne $approved[0]} @two;
      $reject{$reject} = 1;
    }
  }

  # TO DO:
  # also modernize symbols, e.g TRA@ => TRA

  my $gsm = new_gsm_lite();
  my $hgnc = $gsm->hgnc();
  my $approved = $hgnc->get_approved_symbols();
  foreach my $g (@{$approved}) {
    $gsm->add_gene("-gene" => $g);
  }

  my @cleaned;
  foreach my $g (@{$genes_raw}) {
    unless ($reject{$g}) {
      my $hit = $gsm->find($g);
      if ($g eq $hit) {
	push @cleaned, $g;
      } else {
	# modernize
	push @cleaned, $hit;
      }
    }
  }

  write_simple_file(\@cleaned, "genes_cleaned.txt");
}

sub build_intervals {

  my %manual;
  while (<DATA>) {
    chomp;
    next if /^#/;
    my @f = split /,/;
    die unless @f == 5;
    my ($genome,$gene,$chr,$start,$end) = @f;
    $manual{$genome}{$gene} = [ $chr, $start, $end ];
  }

  my $genes = read_simple_file($FLAGS{"build-intervals"});
  my $gsm = new_gsm_lite();
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $f_refflat = $config_genome->{REFSEQ_REFFLAT};

  my $rf = new RefFlatFile(
			   "-strip_sharp_annotations" => 1,
			   "-canonical_references_only" => 1,
#			   "-missing_genes_ok" => 1
			  );
  $rf->parse_file(
		  "-refflat" => $f_refflat,
		  "-type" => "refgene",
		 );
  foreach my $r (@{$rf->rows}) {
    $gsm->add_gene("-gene" => $r->{gene});
  }

  my $outfile = $FLAGS{out} || sprintf 'oncogene_intervals_%s.tab', $genome;

  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   gene
					   chr
					   start
					   end
					   length
					)
				      ],
			 "-auto_qc" => 1,
			);


  foreach my $g_query (@{$genes}) {
    if (my $g_refflat = $gsm->find($g_query)) {
      my ($chr, $start, $end) = $rf->find_gene_interval(
							"-gene" => $g_refflat,
							"-single" => 1,
							"-buffer" => 0,
							"-buffer-upstream" => $BUFFER_UPSTREAM,
							"-multi-ok" => 1,

						       );
      if (ref $chr) {
	foreach my $result (@{$chr}) {
	  my %r;
	  $r{gene} = $g_refflat;
	  @r{qw(chr start end)} = @{$result};
	  $r{length} = ($r{end} - $r{start}) + 1;
	  $rpt->end_row(\%r);
	}
      } else {
	# single simple result
	my %r;
	$r{gene} = $g_refflat;
	$r{chr} = $chr;
	$r{start} = $start;
	$r{end} = $end;
	$r{length} = ($end - $start) + 1;
	$rpt->end_row(\%r);
      }
    } elsif (my $m = $manual{$genome}{$g_query}) {
      my %r;
      my ($chr, $start, $end) = @{$m};
      $start -= $BUFFER_UPSTREAM;
      $end += $BUFFER_UPSTREAM;
      # no strand info, so add to both for safety
      $r{gene} = $g_query;
      $r{chr} = $chr;
      $r{start} = $start;
      $r{end} = $end;
      $r{length} = ($end - $start) + 1;
      $rpt->end_row(\%r);
    } else {
      printf STDERR "ERROR: no refFlat data for %s\n", $g_query;
    }
  }
  $rpt->finish();

}

__DATA__
# manually-specified genomic info from e.g.
# https://www.ncbi.nlm.nih.gov/gene/6955 ("Genomic context" section)
GRCh38,TRA,chr14,21621904,22552132
GRCh38,TRB,chr7,142299011,142813287
GRCh38,TRD,chr14,22422546,22466577
GRCh38,IGH,chr14,105586437,106879844
GRCh38,IGK,chr2,88857361,90235368
GRCh38,IGL,chr22,22026076,22922913
