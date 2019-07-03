#!/bin/env perl
# query annotations in tumor suppressor/oncogene database

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
use GeneSymbolMapper qw(new_gsm_lite);
use FileUtils qw(read_simple_file);

my %FLAGS;
my @clopts = (
	      "-file=s",
	      "-f-gene=s",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

query();

sub query {
  my $ts_onco_db = new TSOncoDB(
				"-gsm" => new_gsm_lite()
			       );
  my $infile = $FLAGS{file} || die "-file";
  my $f_gene = $FLAGS{"f-gene"} || die "-f-gene";

  my $outfile = basename($infile) . ".ts_onco.tab";

  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   qw(
					   in_ts_onco_db
					   is_LoF
					   is_GoF
					    )
					  ],
  			      "-auto_qc" => 1,
			     );
  while (my $row = $df->get_hash()) {
    my $gene = $row->{$f_gene} || dump_die($row, "no $f_gene");
    my $in_db;
    if ($ts_onco_db->find_row("-gene" => $gene)) {
      $in_db = 1;
      $row->{is_LoF} = $ts_onco_db->is_lof("-gene" => $gene) || 0;
      $row->{is_GoF} = $ts_onco_db->is_gof("-gene" => $gene) || 0;
    } else {
      $in_db = 0;
      $row->{is_LoF} = $row->{is_GoF} = "";
    }
    $row->{in_ts_onco_db} = $in_db;
    $rpt->end_row($row);
  }
  $rpt->finish();

}

