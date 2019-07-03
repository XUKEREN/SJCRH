#!/bin/env perl
# ExAC coverage utilities

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

#use lib "/nfs_exports/apps/internal/script_bin/packages/edmonson_perl_lib/";
#use DevelopmentPath;
# optionally replaces production path with development version

use Getopt::Long;

use MiscUtils qw(dump_die);
use TabixPrep qw(tabix_concatenate);
use TARTANUtils qw(tartan_genome_put_helper);
# use DelimitedFile;
# use Reporter;

my $EXAC_COVERAGE_FTP_ROOT = 'ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/coverage/';
my @AVAILABLE_CHROMS = (1 .. 22, qw(X Y));
use constant TARTAN_PARAM => "EXAC_COVERAGE_DIR";

my %FLAGS;
GetOptions(\%FLAGS,
	   "-mirror",
	   "-tabix",
	   "-tartan-put=s",
	   "-genome=s",
	  );

if ($FLAGS{mirror}) {
  mirror_coverage_files();
} elsif ($FLAGS{tabix}) {
  build_tabix();
} elsif ($FLAGS{"tartan-put"}) {
  tartan_put();
} else {
  die;
}

sub mirror_coverage_files {
  foreach my $chr (@AVAILABLE_CHROMS) {
    my $bn = sprintf 'Panel.chr%s.coverage.txt.gz', $chr;
    unless (-s $bn) {
      my $remote = $EXAC_COVERAGE_FTP_ROOT . "/" . $bn;
      my $cmd = "wget $remote";
      system $cmd;
      if ($?) {
	unlink $bn;
	die "$cmd failed with $?";
      }
    }
  }
}

sub build_tabix {
  my @infiles;
  foreach my $chr (@AVAILABLE_CHROMS) {
    my $bn = sprintf 'Panel.chr%s.coverage.txt.gz', $chr;
    if (-s $bn) {
      push @infiles, $bn;
    } else {
      die "where is $bn";
    }
  }

  my $outfile = "Panel.coverage.txt.gz";

  tabix_concatenate(
		    "-in" => \@infiles,
		    "-out" => $outfile,
		    "-sorted" => 1,
		    "-header_chr" => '#chrom',
		    "-header_start" => 'pos',
		    "-ping" => 100000,
#		    "-max" => 10,
		   );

}

sub tartan_put {
  my $out_dir = $FLAGS{"tartan-put"};

  my $cmd = sprintf 'tartan_import_helper.pl -genome %s -put %s -param %s',
    $FLAGS{genome}, $out_dir, TARTAN_PARAM;
  system $cmd;
}

