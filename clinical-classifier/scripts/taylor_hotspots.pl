#!/bin/env perl

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

my $F_TAYLOR_PDF_S3 = "s3.txt";
# text extracted from PDF of supplementary table 3
# http://www.nature.com/nbt/journal/v34/n2/full/nbt.3391.html#contrib-auth

#my $F_TAYLOR_MAF = "minimalist_test_maf.txt";
# missing many variants
my $F_TAYLOR_MAF = "pancan_unfiltered.maf";
# https://github.com/taylor-lab/hotspots/

my %FLAGS;
my @clopts = (
	      "-build-false-positives",
	     );

GetOptions(
	   \%FLAGS,
	   @clopts
	  );

if ($FLAGS{"build-false-positives"}) {
  build_false_positives();
}

sub build_false_positives {
  open(IN, $F_TAYLOR_PDF_S3) || die;

  my %fp;
  # limited data from S3 PDF info
  while (<IN>) {
    my $raw = $_;
    next if /,/;
    # header
    my @f = split /\s+/, $_;
    my ($gene, $codon) = @f;
    my $reason = $f[$#f];
    $codon =~ s/_.*//;
    # splice
    next if $codon =~ /^c\./;
    # need protein
    if ($codon =~ /^[A-Z\*]\d+$/) {
      die "duplicate" if $fp{$gene}{$codon};
      $fp{$gene}{$codon} = $reason || die;
      printf STDERR "track %s %s\n", $gene, $codon;
#      die "$gene $codon";
    } else {
      printf STDERR "skipping codon:%s raw:%s\n", $codon, $raw;
    }
  }

  # find equivalent rows in input MAF
  my $df = new DelimitedFile("-file" => $F_TAYLOR_MAF,
			     "-headers" => 1,
			     );
  my $outfile = "taylor_s3_false_positives.maf";
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   qw(
					       false_positive_reason
					    )
					  ],
  			      "-auto_qc" => 1,
			     );

  my %saw;
  while (my $row = $df->get_hash()) {
    my $gene = $row->{Hugo_Symbol} || next;
    # not always present
    if (my $aa = $row->{HGVSp_Short}) {
      my $aa_raw = $aa;
      $aa =~ s/^p\.//;
      if ($aa =~ /^([A-Z\*]\d+)/) {
	my $key = $1;
	printf STDERR "search %s %s\n", $gene, $key;
	if (my $reason = $fp{$gene}{$key}) {
	  $saw{$gene}{$key} = 1;
	  $row->{false_positive_reason} = $reason;
	  $rpt->end_row($row);
	}
      }
    }
  }
  $rpt->finish();

  foreach my $gene (sort keys %fp) {
    foreach my $codon (sort keys %{$fp{$gene}}) {
      if ($saw{$gene}{$codon}) {
	printf STDERR "success for %s %s\n", $gene, $codon;
      } else {
	printf STDERR "ERROR: can't find MAF data for %s %s\n", $gene, $codon;
      }
    }
  }



}

