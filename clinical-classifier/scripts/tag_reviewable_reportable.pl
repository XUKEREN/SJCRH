#!/bin/env perl
# tag germline-reportable and germline-reviewable genes
# MNE 12/2018

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list get_hash_option);
use FileUtils qw(read_simple_file);
use DelimitedFileHP;
use Reporter;
use TdtConfig;
use GeneSymbolMapper qw(new_gsm_lite);

my %FLAGS;
my @clopts = (
	      "-file=s",
	      "-genome=s",
	      "-f-gene=s",
	      "-out=s",

	      "-strip=s",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $infile = $FLAGS{file} || die "-file";

my $genome = $FLAGS{genome} || die "-genome";
my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
my $f_reportable = $config_genome->{CLINCLS_GERMLINE_REPORTABLE_GENES} || die;
my $f_gene = $FLAGS{"f-gene"} || die "-f-gene";

my $gsm_reviewable = new_gsm_lite();
my $f_reviewable = $config_genome->{CLINCLS_CANCER_RELATED_GENES_FILE} || die;
open(IN, $f_reviewable) || die;
while (<IN>) {
  chomp;
  my @f = split /\t/, $_;
  die unless @f == 2;
  $gsm_reviewable->add_gene("-gene" => $f[0]);
}

my $gsm_reportable = new_gsm_lite();
my $gl = read_simple_file($f_reportable) || die;
foreach my $gene (@{$gl}) {
  $gsm_reportable->add_gene("-gene" => $gene);
}


my $df = new DelimitedFileHP(
			     "-file" => $infile,
			     "-headers_extra" => [
						  "is_reviewable_gene",
						  "is_reportable_gene"
						 ]
			    );

my $strip_to_reportable;
my $strip_to_reviewable;
if (my $st = $FLAGS{strip}) {
  if ($st eq "reportable") {
    $strip_to_reportable = 1;
  } elsif ($st eq "reviewable") {
    $strip_to_reviewable = 1;
  } else {
    die "-strip must be reviewable or reportable";
  }
}

my $outfile = $FLAGS{out};
unless ($outfile) {
  $outfile = basename($infile);
  if ($strip_to_reviewable) {
    $outfile .= ".reviewable.tab";
  } elsif ($strip_to_reportable) {
    $outfile .= ".reportable.tab";
  } else {
    $outfile .= ".reviewable_reportable.tab";
  }
}

$df->write_init("-file" => $outfile);
$df->prepare_query("-fields" => [ $f_gene ]);
my $dropped = 0;
my $match_same = 0;
my $match_alt = 0;
my $count = 0;
my %alt;

while ($df->next_row()) {
  printf STDERR "%d (%d/%d)...\n", $count, $match_same, $match_alt if ++$count % 25000 == 0;
  my $gene_raw = $df->get_query()->[0];
  my $is_reviewable = 0;
  my $is_reportable = 0;

  if ($gene_raw) {
    my @genes = split /,/, $gene_raw;
    # Annovar+ sometimes reports multiple genes in case of ambiguity
    foreach my $gene (@genes) {
      foreach my $thing ([$gsm_reportable, \$is_reportable],
			 [$gsm_reviewable, \$is_reviewable]
			) {
	my ($gsm, $flag_ref) = @{$thing};

	if (my $found = $gsm->find($gene)) {
	  $$flag_ref = 1;
	  if ($gene eq $found) {
	    $match_same++;
	  } else {
	    $match_alt++;
	    unless ($alt{$gene}{$found}) {
	      printf STDERR "alt match %s => %s\n", $found, $gene;
	      $alt{$gene}{$found} = 1;
	    }
	  }
	}
      }
    }

    my $usable;
    die "reportable gene not reviewable" if $is_reportable and not($is_reviewable);
    # grievous eror

    if ($strip_to_reportable) {
      $usable = $is_reportable;
    } elsif ($strip_to_reviewable) {
      $usable = $is_reviewable;
    } else {
      $usable = 1;
    }

    if ($usable) {
      my %extra;
      $extra{is_reportable_gene} = $is_reportable;
      $extra{is_reviewable_gene} = $is_reviewable;
      $df->write_row("-extra" => \%extra);
    } else {
      $dropped++;
    }
  }
}
$df->write_finish();

printf STDERR "match_same:%d match_alt:%d dropped:%d\n", $match_same, $match_alt, $dropped;


