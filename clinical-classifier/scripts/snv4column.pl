#!/bin/env perl
# split a SNV4 annotation into individual columns

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

my $F_CHR = "Chr";
my $F_POS = "WU_HG19_Pos";
my $F_RA = "ReferenceAllele";
my $F_VA = "MutantAllele";

my %FLAGS;
my @clopts = (
	      "-field=s",
	      # field label

	      "-snv4",
	      "-snv5",
	      # including sample name

	      "-file=s",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $infile = get_hash_option(\%FLAGS, "file");
my $outfile = basename($infile) . ".alleles.tab";

my $df = new DelimitedFile(
			   "-file" => $infile,
			   "-headers" => 1,
			  );

my $rpt = $df->get_reporter(
			    "-file" => $outfile,
			    "-extra" => [
					 $F_CHR,
					 $F_POS,
					 $F_RA,
					 $F_VA
					],
			    "-auto_qc" => 1,
			   );

# while (my $row = $df->next("-ref" => 1)) {  # headerless
my $fname = $FLAGS{field} || die "-field";
while (my $row = $df->get_hash()) {
  my $v = $row->{$fname} || die;
  my @f = split /\./, $v;
  if ($FLAGS{snv5}) {
    shift @f;
  } elsif ($FLAGS{snv4}) {
    # no
  } else {
    die;
  }
  die unless @f == 4;

  @{$row}{$F_CHR, $F_POS, $F_RA, $F_VA} = @f;

  foreach my $f ($F_RA, $F_VA) {
    $row->{$f} = "-" if $row->{$f} eq "_" or $row->{$f} eq "";
    # Zhaoming
  }

  $rpt->end_row($row);

}

$rpt->finish();

