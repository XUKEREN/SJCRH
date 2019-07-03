#!/bin/env perl
# generate annotations for use w/CLINCLS_GL_REPORTABLE_GENE_ANNOTATION_FILE

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

my %FLAGS;
my @clopts = (
	      "-annotations=s",

	      "-genome=s",
	      "-refflat=s",
	      "-ger=s",
	      # manual overrides: taken from genome config by default
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $f_annot = $FLAGS{annotations} || die "-annotations";

my $df = new DelimitedFile("-file" => $f_annot,
			   "-headers" => 1,
			  );

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $outfile = sprintf "GermlineReportableGeneAnnotation.txt.%d_%02d_%02d", 1900+$year, $mon + 1, $mday;

my $rpt = new Reporter(
		       "-file" => $outfile,
		       "-delimiter" => "\t",
		       "-labels" => [
				     qw(
					 Gene
					 TruncationGold
				      )
				    ],
		       "-auto_qc" => 1,
		      );


# while (my $row = $df->next("-ref" => 1)) {  # headerless
my %saw;
my %gene2row;

while (my $row = $df->get_hash()) {
  my $gene = $row->{"HGNC Name"} || die;
  my @things;

  $gene =~ s/^\s+//;
  $gene =~ s/\s+$//;

  foreach my $f (qw(class1 class2)) {
    my $v = $row->{$f} || next;
    if ($gene eq "EGFR" and $v eq "CMPLX") {
      $v = "LoF";
      push @things, "GoF";
      # complex case, track both
    } elsif ($gene eq "SHOC2") {
      $v = "GoF" if $v eq "DISREG";
      next if $v eq "aberrant regulation";
      # for compatibility w/earlier config file, which was set to N
    }

    $v = "GoF" if $v eq '"assume LoF is bad, but possible GoF"';
    push @things, $v;
  }

  my %things = map {$_, 1} @things;

  my $truncation_gold;
  if (@things == 1) {
    if ($things[0] eq "GoF") {
      $truncation_gold = "N";
    } elsif ($things[0] eq "LoF") {
      $truncation_gold = "Y";
    } else {
      die "unknown value " . $things[0];
    }
  } elsif (@things == 2 and $things{LoF} and $things{GoF}) {
    $truncation_gold = 2;
    # new status indicating both TS and oncogene
    # (should help generate more sophisticated germline CNV config)
  } else {
    die @things;
  }

  die "duplicate $gene" if $saw{$gene};
  $saw{$gene} = 1;

  my %r;
  $r{Gene} = $gene;
  $r{TruncationGold} = $truncation_gold;
  $gene2row{$gene} = \%r;

  $rpt->end_row(\%r);
}

#
#  check symbols vs. refFlat and GENE_EXON_REGION and add
#  duplicate entries for alternate symbols
#
my $gsm_refflat = new_gsm();
my $refflat = $FLAGS{"refflat"} || die "-refflat";
# TO DO: get from genome config when rdev fixed
$gsm_refflat->populate_refflat("-refflat" => $refflat);

my %cloned;

foreach my $gene (sort keys %gene2row) {
  if ($gsm_refflat->contains($gene)) {
    # identical, OK
  } elsif (my $gene_rf = $gsm_refflat->resolve("-symbol" => $gene)) {
    printf STDERR "reportable symbol %s matches refFlat symbol %s\n",
      $gene, $gene_rf;
    my %clone = %{$gene2row{$gene}};
    $clone{Gene} = $gene_rf;
    $cloned{$gene_rf} = 1;
    $rpt->end_row(\%clone);
  }
}

my $gsm_ger = new_gsm();
my $ger_dir = $FLAGS{ger} || die "-ger";
$gsm_ger->populate_gene_exon_region("-dir" => $ger_dir);
foreach my $gene (sort keys %gene2row) {
  if ($gsm_ger->contains($gene)) {
    # identical, OK
  } elsif (my $gene_rf = $gsm_ger->resolve("-symbol" => $gene)) {
    printf STDERR "reportable symbol %s matches GENE_EXON_REGION symbol %s\n",
      $gene, $gene_rf;
    my %clone = %{$gene2row{$gene}};
    $clone{Gene} = $gene_rf;
    $rpt->end_row(\%clone) unless $cloned{$gene_rf};
  }
}

$rpt->finish();


sub new_gsm {
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


