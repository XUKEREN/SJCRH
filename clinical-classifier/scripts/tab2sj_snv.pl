#!/bin/env perl
# convert tab-delimited text to SJ format for SNVs
# MNE 1/2015

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die);
use DelimitedFile;
# use Reporter;

my %FLAGS;

GetOptions(\%FLAGS,
	   "-file=s",
	   # single file

	   "-type=i",
	   # various formats encountered
    );

my @FILES = $FLAGS{file} || die "-file";

my $file_type = $FLAGS{type} || die "-type";

convert_type_1() if $file_type == 1;

sub convert_type_1 {

  foreach my $infile (@FILES) {
    my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
    my $outfile = basename($infile) . ".sj.tab";
    my $rpt = get_rpt($df, $outfile);
    while (my $row = $df->get_hash()) {
#      dump_die($row, "Debug", 1);
      $row->{GeneName} = $row->{"Gene.refGene"} || die;
      $row->{Chr} = $row->{CHROM} || die;
      $row->{WU_HG19_Pos} = $row->{POS} || die;
      $row->{ReferenceAllele} = $row->{REF} || die;

      my $function = $row->{"ExonicFunc.refGene"} || die;
      my $class;
      if ($function eq "nonsynonymous_SNV") {
	$class = "missense";
      } elsif ($function eq "stopgain") {
	$class = "nonsense";
      } else {
	die "unhandled function $function";
      }
      $row->{Class} = $class;

      my $aac = $row->{"AAChange.refGene"} || die;
      my @things = split /:/, $aac;
      die unless @things == 5;
      my ($gene, $nm, $thing, $codon, $protein) = @things;
      die unless $nm =~ /^NM_/;

      $row->{mRNA_acc} = $nm;
      $row->{AAChange} = $protein;

      my @alts = split ";", $row->{ALT} || die;
      dump_die($row, "multiple alternate alleles", 1) if @alts > 1;
      foreach my $alt (@alts) {
	$row->{MutantAllele} = $alt;
	$rpt->end_row($row);
      }


    }
    $rpt->finish();
  }

}

sub get_rpt {
  my ($df, $outfile) = @_;

  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   qw(
					       GeneName
					       Chr
					       WU_HG19_Pos
					       ReferenceAllele
					       MutantAllele
					       Class
					       AAChange
					       mRNA_acc
					    )
					  ]
			     );
  return $rpt;
}
