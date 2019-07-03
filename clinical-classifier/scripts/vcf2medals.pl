#!/bin/env perl
# meta wrapper script for all steps to run a VCF through medal ceremony
# MNE 8/2015
#
# vcftools
# tabix
#
# annovar
# snv-annovar
# clinical-classifier   *** NOTE: on research, DEV ENVIRONMENT REQUIRED!
#
# to do:
# - add tags to outfile stack files to e.g. flag for cleanup

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list);
use OutfileStack;
use FileUtils qw(find_binary);
# use DelimitedFile;
# use Reporter;

my %FLAGS;
GetOptions(
	   \%FLAGS,
	   "-genome=s",

	   "-file=s",
	   "-files=s",

	   "-germline-reviewable",
	  );
# TO DO: -cluster

my $files = build_argv_list(
			    "-flags" => \%FLAGS,
			    "-single" => "file",
			    "-set" => "files"
			   );

my $genome = $FLAGS{genome} || die "-genome";
my $germline_reviewable = $FLAGS{"germline-reviewable"};
printf STDERR "filtering to germline-reviewable only?: %s\n", $germline_reviewable ? "y" : "n";

foreach my $bin (qw(
		     tabix
		     vcf2tab.pl
		     bambino2annovar.pl
		     annotate_variation.pl
		     double_medal.pl
		  )) {
  # verify module loads
  find_binary($bin, "-die" => 1);
}


foreach my $vcf (@{$files}) {
  my $os = new OutfileStack(
			    "-start_file" => $vcf,
			    "-use_basename" => 1
			   );

  $os->add_level("-suffix" => ".cooked.tab.gz");
  unless ($os->current_valid()) {
    # convert to tab-delimited w/special Annovar-prep filtering
    my $cmd = sprintf 'vcf2tab.pl -annovar -file %s -gz', $vcf;
    run_cmd($cmd);
  }

  if ($germline_reviewable) {
    # optionally filter to germline-reviewable regions only
    $os->add_level("-suffix" => ".reviewable.tab");

    unless ($os->current_valid()) {
      my $cmd = sprintf 'filter_region.pl -germline-reviewable -genome %s -file %s',
	$genome,
	  $os->get_previous_file();
      run_cmd($cmd);
    }
  }

  $os->add_level("-suffix" => ".annovar_merged.tab");
  unless ($os->current_valid()) {
    # Annovar+
    my $cmd = sprintf 'bambino2annovar.pl -genome %s -file %s -scratch', $genome, $os->get_previous_file();
    run_cmd($cmd);
  }

  $os->add_level("-suffix" => ".mp.tab");
  unless ($os->current_valid()) {
    # medal prep
    my $cmd = sprintf 'annovar2medals.pl -genome %s -file %s', $genome, $os->get_previous_file();
    run_cmd($cmd);
  }

  $os->add_level("-suffix" => ".medals.tab.medals.tab");
  unless ($os->current_valid()) {
    # medals
    my $cmd = sprintf 'double_medal.pl -file %s', $os->get_previous_file();
    run_cmd($cmd);
  }


}

sub run_cmd {
  my ($cmd) = @_;
  printf STDERR "%s: running %s\n", scalar(localtime), $cmd;
  system $cmd;
  die "$cmd exited with $?" if $?;
}
