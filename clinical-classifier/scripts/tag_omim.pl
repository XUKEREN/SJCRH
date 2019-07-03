#!/bin/env perl
# tag OMIM database data
# MNE 8/2018

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list get_hash_option unique_ordered_list);
use DelimitedFile;
use Reporter;
use OMIM;

my %FLAGS;
my @clopts = (
	      "-file=s",
	      "-out=s",
	      "-gene-column=s",
	      "-chr-column=s",
	      "-no-chr",

	      "-glob=s",

	      "-omim-prep",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

if ($FLAGS{"omim-prep"}) {
  omim_prep();
} else {
  tag_omim();
}

sub tag_omim {
  my $infile = $FLAGS{file} || die "-file";
  my $outfile = $FLAGS{out} || basename($infile) . ".omim.tab";
  my $f_gene = $FLAGS{"gene-column"} || die "-gene-column";
  my $f_chr = $FLAGS{"chr-column"};
  unless ($f_chr or $FLAGS{"no-chr"}) {
    die "specify -chr-column COLUMN for chromosome, or -no-chr\n";
  }


  my $omim = new OMIM();

  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   qw(
					       OMIM_ID
					       OMIM_gene
					       OMIM_phenotypes
					    )

					  ],
  			      "-auto_qc" => 1,
			     );

  while (my $row = $df->get_hash()) {
    die "where is $f_gene" unless exists $row->{$f_gene};
    my $omim_id = "";
    my $omim_gene = "";
    my $omim_pheno = "";

    if (my $gene = $row->{$f_gene}) {
      if (my $hits = $omim->find("-gene" => $gene)) {

	if (@{$hits} > 1 and $f_chr and my $user_chr = $row->{$f_chr}) {
	  $user_chr =~ s/^chr//i;
	  my @filtered;
	  foreach my $r (@{$hits}) {
	    my $db_chr = $r->{Chromosome} || die;
	    $db_chr =~ s/^chr//i;
	    push @filtered, $r if $db_chr eq $user_chr;
	  }
	  if (@filtered) {
	    $hits = \@filtered;
	  } else {
	    dump_die($hits->[0], "removed all matches w/user chr $user_chr", 1);
	  }
	}

	$omim_id = listify($hits, "Mim Number");
	$omim_gene = listify($hits, "Approved Symbol");
	$omim_pheno = listify($hits, "Phenotypes");
      }
    }
    $row->{OMIM_ID} = $omim_id;
    $row->{OMIM_gene} = $omim_gene;
    $row->{OMIM_phenotypes} = $omim_pheno;

    $rpt->end_row($row);
  }
  $rpt->finish();
}

sub listify {
  my ($rows, $column) = @_;
  return join ";", @{unique_ordered_list([map {$_->{$column}} @{$rows}])};
  # use semicolon rather than a command, because semicolon is used to
  # separate entries when there are multiple phenotypes in a record
}


sub omim_prep {
  # tweak OMIM file distributions for simple tab-delimited use:
  # - strip extraneous header lines
  # - remove leading spaces from first column name
  # - strip trailing comment lines
  my @f_omim = glob($FLAGS{"glob"} || die "-glob");

  foreach my $f_omim (@f_omim) {
    my $of = $f_omim . ".tab";
    my $wf = new WorkingFile($of);
    my $fh = $wf->output_filehandle();
    open(IN, $f_omim) || die;
    my $in_body;
    my $body_lines = 0;
    while (<IN>) {
      if ($in_body) {
	if (/^#/) {
	  # trailing comments, quit
	  last;
	} else {
	  print $fh $_;
	  $body_lines++;
	}
      } else {
	# header
	if (/^#/) {
	  if (/\t/) {
	    # first comment line with tabs is the header line
	    s/^\#\s+//;
	    # remove leading comment sign and whitespace
	    print $fh $_;
	    $in_body = 1;
	  }
	} else {
	  die "too far";
	}
      }
    }
    die "ERROR: no body lines written" unless $body_lines;
    $wf->finish();
  }

}
