#!/bin/env perl
# HGMD database utilities

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
use GenomeUtils qw(reverse_complement);

my %FLAGS;
my @clopts = (
	      "-cook=s",
	      # digest HGMD mutation position info into SNV4 style
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

cook_variants() if $FLAGS{cook};

sub cook_variants {
  my $infile = $FLAGS{cook} || die;
  my $outfile = $infile . ".cooked.tab";
  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   qw(
					       Chr
					       WU_HG19_Pos
					       ReferenceAllele
					       MutantAllele
					    )
					  ],
  			      "-auto_qc" => 1,
			     );

  while (my $row = $df->get_hash()) {
    my $type = $row->{base} || die;
    my $strand = $row->{strand} || die;
    my $ra = $row->{wildbase};
    my $va = $row->{mutbase};
    my $chr = $row->{chromosome} || die;
    my $pos = $row->{coordstart} || die;
    if ($strand eq "-") {
      foreach ($ra, $va) {
	$_ = reverse_complement($_) if $_;
      }
    }

    if ($type eq "D") {
      # deletion
      die unless $ra;
      die if $va;
      $va = "-";
    } elsif ($type eq "I") {
      # insertion
      die unless $va;
      die if $ra;
      $ra = "-";
    } elsif ($type eq "M" or $type eq "S") {
      # substitution, splice
      die unless $ra and $va and length($ra) == length($va);
    } elsif ($type eq "X") {
      # MNV or complex
      die unless $ra and $va;
    } elsif ($type eq "R") {
      # regulatory (?)
      die unless $ra and $va;
    } else {
      die "unhandled type $type";
    }

    $row->{Chr} = $chr;
    $row->{WU_HG19_Pos} = $pos;
    $row->{ReferenceAllele} = $ra;
    $row->{MutantAllele} = $va;
    $rpt->end_row($row);
  }
  $rpt->finish();



}
