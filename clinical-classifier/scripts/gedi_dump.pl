#!/bin/env perl
# build collections of GeDI variants for use in clinical classification.
# This is a general list of PCGP variants rather than the curated
# set of recurrent variants.

use strict;
use warnings;
use Getopt::Long;

use DBTools qw(get_dbi_gedi export_query_to_flatfile);
use TdtConfig;
use MiscUtils qw(dump_die);
use DelimitedFile;
use TabixPrep;
use Variant qw(INDEL_CHAR);

use constant TARTAN_PARAM => 'CLINCLS_GEDI_TABIX_DIR';

my %FLAGS;
my @clopts = (
	      "-genome=s",
	      "-dump-old",
	      "-legacy-somatic-to-tabix",
	      # convert old legacy files to single tabix file

	      "-tartan-index=s",

	      # TO DO:
	      # - query db for new data and write to tabix
	      #   (i.e. the new standard going forward)
             );
GetOptions(
           \%FLAGS,
           @clopts
          );

if (my $out_dir = $FLAGS{"tartan-index"}) {
  my $genome = $FLAGS{genome} || die "-genome";

  my $cmd = sprintf 'tartan_import_helper.pl -genome %s -put %s -param %s',
    $genome, $out_dir, TARTAN_PARAM;

  system $cmd;

  exit(0);
}


dump_old() if $FLAGS{"dump-old"};
legacy_somatic_to_tabix() if $FLAGS{"legacy-somatic-to-tabix"};

sub dump_old {
  # obsolete: old method
  my $dbi = get_dbi_gedi("-type" => "research");

  my @views = (
	       "indel_find_t1",
	       # VALID SOMATIC indels
	       "snp_find_t1",
	       # VALID SOMATIC SNVs, tier 1
	       "snp_find_t2",
	       # VALID SOMATIC SNVs, tier 2
	       "snp_find_t3",
	       # VALID SOMATIC SNVs, tier 3
	      );

  foreach my $view (@views) {
    printf STDERR "exporting %s...\n", $view;
    export_query_to_flatfile("-dbi" => $dbi,
			     "-table" => $view);
  }
  $dbi->disconnect();
}

sub legacy_somatic_to_tabix {
  #
  #  convert legacy flatfiles to single tabix file for use w/batch annotation
  #
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find\
 config for $genome";

  my @rows;
  foreach my $var (
		   "CLINCLS_GEDI_SNP_FIND_T1_FILE",
		   "CLINCLS_GEDI_SNP_FIND_T2_FILE",
		   "CLINCLS_GEDI_SNP_FIND_T3_FILE",
                   "CLINCLS_GEDI_INDEL_FIND_T1_FILE"
		  ) {
    my $file = $config_genome->{$var} || die "no file for $var";

    my $df = new DelimitedFile(
			       "-file" => $file,
			       "-headers" => 1,
			      );
    while (my $row = $df->get_hash()) {
      push @rows, $row;
    }
  }

  rows2tabix(
	     "-rows" => \@rows,
	     "-somatic" => 1
	    );
}

sub rows2tabix {
  my (%options) = @_;
  my $rows = $options{"-rows"} || "-rows";
  my $somatic_mode = $options{"-somatic"};

  my $outfile;
  if ($somatic_mode) {
    $outfile = sprintf 'gedi_somatic_indelpos_cooked.tab.gz';
  } else {
    die;
  }

  my @required_fields = qw(
			    chromosome
			    pos
			    reference_allele
			    non_reference_allele
			    origin
			    official_project
			 );

  my $tp = new TabixPrep(
                         "-outfile" => $outfile,
                         "-headers" => \@required_fields,
                         "-header_chr" => "chromosome",
                         "-header_start" => "pos",
                        );

  foreach my $row (@{$rows}) {
    foreach my $f (@required_fields) {
      dump_die($row, "missing field $f") unless exists $row->{$f};
    }

    my $origin = $row->{origin} || die;
    die unless $origin eq "GERMLINE" or $origin eq "SOMATIC";
    next if $somatic_mode and $origin ne "SOMATIC";

    my ($chr, $pos, $ra, $va) = @{$row}{qw(chromosome pos reference_allele non_reference_allele)};
    foreach ($ra, $va) {
      s/\-//g;
      # NOTE: some gedi complex coordinates appear broken, e.g.
      # 5.35874572.------TTAC.TCGCCCTGCA
      # ...should be 35874568?
    }
    foreach ($ra, $va) {
      $_ = INDEL_CHAR unless $_;
    }

    my $v = new Variant();
    $v->import_generic(
		       "-reference-name" => $chr,
		       "-base-number" => $pos,
		       "-reference-allele" => $ra,
		       "-variant-allele" => $va
		      );
    if ($v->is_insertion()) {
      # for insertions, cook from bambino to standard format
      # as this is the coordinate system expected by batch tabix.
      $pos--;
    }

    my %r = %{$row};
    $r{pos} = $pos;
    $r{reference_allele} = $ra;
    $r{non_reference_allele} = $va;
    # update to reflect insertion correction, if required

    $tp->add_row("-row" => \%r);

  }
  $tp->finish();

}
