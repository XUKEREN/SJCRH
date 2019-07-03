#!/bin/env perl
# add HGNC annotation to tab-delimited file

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
use GeneSymbolMapper qw(new_gsm_lite);

my %FLAGS;
my @clopts = (
	      "-file=s",
	      "-gene-column=s",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $infile = $FLAGS{file} || die "-file";
my $f_gene = $FLAGS{"gene-column"} || die "-gene-column";
my $gsm = new_gsm_lite();

my $hgnc_genes = $gsm->get_hgnc_approved();
foreach my $g (@{$hgnc_genes}) {
  $gsm->add_gene("-gene" => $g);
}

my $df = new DelimitedFile(
			   "-file" => $infile,
			   "-headers" => 1,
			  );

my $outfile = basename($infile) . ".hgnc.tab";
my $rpt = $df->get_reporter(
			    "-file" => $outfile,
			    "-extra" => [
					 qw(
					     HUGO
					     is_approved
					     is_blacklisted
					  )
					],
			    "-auto_qc" => 1,
			   );

while (my $row = $df->get_hash()) {
  die "where is $f_gene" unless exists $row->{$f_gene};
  my $sym = $row->{$f_gene};
  my $hugo = $gsm->find($sym) || "";
  $row->{HUGO} = $hugo;
  my $is_approved;
  my $blacklist;
  if ($hugo) {
    $is_approved = 1;
    $blacklist = "";
  } else {
    $is_approved = 0;
    $blacklist = $gsm->hgnc->is_blacklisted($sym) || "";
#    dump_die($row, $blacklist) if $blacklist;
  }

  $row->{is_blacklisted} = $blacklist;
  $row->{is_approved} = $is_approved;

  $rpt->end_row($row);
}

$rpt->finish();
