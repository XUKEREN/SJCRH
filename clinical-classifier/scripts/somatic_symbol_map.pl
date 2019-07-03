#!/bin/env perl
# describe me

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list get_hash_option);
use FileUtils qw(read_simple_file);
use DelimitedFile;
use Reporter;
use TdtConfig;
use GeneSymbolMapper;

my %FLAGS;
my @clopts = (
	      "-genome=s",
	      "-map-list-to-refflat=s",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $genome = $FLAGS{genome} || die "-genome";
my $CONFIG_GENOME = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";

map_list() if $FLAGS{"map-list-to-refflat"};

sub map_list {
  my $gsm = get_gsm();
  my $refflat = $CONFIG_GENOME->{REFSEQ_REFFLAT} || die "no refflat";
  $gsm->populate_refflat("-refflat" => $refflat);

  my $infile = $FLAGS{"map-list-to-refflat"} || die;
  my $outfile = sprintf "%s.%s", basename($infile), $genome;

  my $gene_list = read_simple_file($infile);
  my $wf = new WorkingFile($outfile);
  my $fh = $wf->output_filehandle();

  my $count_same = 0;
  my $count_changed = 0;
  foreach my $gene (@{$gene_list}) {
    my $gene_out;
    if ($gsm->contains($gene)) {
      $gene_out = $gene;
      $count_same++;
    } elsif (my $mapped = $gsm->resolve("-symbol" => $gene)) {
      printf STDERR "map %s => %s\n", $gene, $mapped;
      $gene_out = $mapped;
      $count_changed++;
    } else {
      die "ERROR: can't resolve $gene";
    }
    printf $fh "%s\n", $gene_out;
  }
  $wf->finish();

  printf "identical:%d changed:%d\n", $count_same, $count_changed;


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
