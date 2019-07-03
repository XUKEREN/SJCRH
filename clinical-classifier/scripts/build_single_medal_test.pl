#!/bin/env perl
# build unique list of variants for medal regression tests
# MNE 10/2015

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

#use lib "/nfs_exports/apps/internal/script_bin/packages/edmonson_perl_lib/";
#use DevelopmentPath;
# optionally replaces production path with development version

use Getopt::Long;

use MiscUtils qw(dump_die);
use DelimitedFile;
use Reporter;
use Counter;

my %FLAGS;
GetOptions(
	   \%FLAGS,
	   "-dir=s",
	  );

my $dir = $FLAGS{dir} || die;
my @files = glob($dir . "/*");
die unless @files;

my $df = new DelimitedFile("-file" => $files[0],
			   "-headers" => 1,
			  );

my @out_headers;

foreach my $h (@{$df->headers_raw}) {
  $h = "WU_HG19_Pos" if $h eq "WU_HG18_Pos";
}

my $outfile = "merged.tab";
my $rpt = new Reporter(
		       "-file" => $outfile,
		       "-delimiter" => "\t",
		       "-labels" => $df->headers_raw(),
		       "-auto_qc" => 1,
		      );

my %saw;
my $c = new Counter(\@files);
foreach my $infile (@files) {
  $c->next($infile);
  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			    );
  while (my $row = $df->get_hash) {
    unless ($row->{WU_HG19_Pos}) {
      $row->{WU_HG19_Pos} = $row->{WU_HG18_Pos} || die;
    }
    my $key = join ".", @{$row}{qw(Chr WU_HG19_Pos ReferenceAllele MutantAllele)};
    $rpt->end_row($row) unless $saw{$key};
    $saw{$key} = 1;
  }
}
$rpt->finish;



