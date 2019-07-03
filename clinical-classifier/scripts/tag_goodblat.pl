#!/bin/env perl
# "goodblat" list tabix building and tagging

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list get_hash_option log_message);
use DelimitedFile;
use Reporter;
use TabixPrep;
use TdtConfig;
use Variant;
use TabixBatchAnnotation;

my $F_TABIX_CHR = "Chr";
my $F_TABIX_POS = "WU_HG19_Pos";
my $F_TABIX_RA = "ReferenceAllele";
my $F_TABIX_VA = "MutantAllele";
my $F_TABIX_GOODBLAT = "GoodBLATFlag";
# standardize

my $F_IN_GOODBLAT = "BLAT_check";

my $F_OUT_GOODBLAT = "GoodBLATFlag";
# tag field output

use constant GOODBLAT => "GoodBLAT";
use constant BADBLAT => "BadBLAT";

my $F_CHR;
my $F_POS;
my $F_RA;
my $F_VA;

my $QUEUE_SIZE = 1000;
my $TABIX_INDEL_EQUIVALENCE_DISTANCE = 25;
my $FIELD_GOODBLAT_TABIX = "__goodblat_tabix";
use constant TABIX_BATCH_SIZE => 500;

my %FLAGS;
my @clopts = (
	      "-genome=s",
	      "-build-tabix",
	      "-goodblat=s",

	      "-file=s",
	      "-tabix-goodblat=s",
	      "-f-chr=s" => \$F_CHR,
	      "-f-pos=s" => \$F_POS,
	      "-f-ra=s" => \$F_RA,
	      "-f-va=s" => \$F_VA,

	      "-good-only",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $genome = $FLAGS{genome} || die "-genome";
my $CONFIG_GENOME = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";

if ($FLAGS{"build-tabix"}) {
  # convert goodblat flatfile to tabix
  build_tabix();
} else {
  # tag
  tag_file();
}

sub build_tabix {
  my $infile = $FLAGS{goodblat};
  # /nfs_exports/genomes/1/projects/WHOLEGENOME/Eval-CAP-NGS-A-2015/BucketIntermediate/SiteEval/xma/tagger_v3/blat/blat.good.bad.new
  unless ($infile) {
    printf STDERR "using config system file, may be out of date\n";
    $infile = $CONFIG_GENOME->{BLATLIST};
  }

  my $outfile = sprintf "%s.gz", basename($infile);

  my @headers = (
		 $F_TABIX_CHR,
		 $F_TABIX_POS,
		 $F_TABIX_RA,
		 $F_TABIX_VA,
		 $F_TABIX_GOODBLAT
		);
  # should we digest to just positions only??

  my $tp = new TabixPrep(
			 "-outfile" => $outfile,
			 "-headers" => \@headers,
			 "-header_chr" => $F_TABIX_CHR,
			 "-header_start" => $F_TABIX_POS,
			 #		       "-input_sorted" => 1,
			);

  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			    );
  while (my $row = $df->get_hash()) {
    my @snv4 = split /\./, $row->{SNV4} || die;
    die unless @snv4 == 4;
    my %r;
    my ($chr, $pos, $ra, $va) = @snv4;
    my $goodblat = $row->{$F_IN_GOODBLAT} || die;

    my $usable = 1;

    if ($usable) {
      $chr =~ s/^chr//i;
      $r{$F_TABIX_CHR} = $chr;
      $r{$F_TABIX_POS} = $pos;
      $r{$F_TABIX_RA} = $ra;
      $r{$F_TABIX_VA} = $va;

      my $flag;
      if ($goodblat eq GOODBLAT) {
	$flag = 1;
      } elsif ($goodblat eq BADBLAT) {
	$flag = 0;
      } else {
	die "unhandled value $goodblat";
      }
      $r{$F_TABIX_GOODBLAT} = $flag;

      $tp->add_row("-row" => \%r);
    }
  }

  $tp->finish();
}

sub tag_file {
  die "fields not defined" unless $F_CHR and $F_POS and $F_RA and $F_VA;
  my $infile = $FLAGS{file} || die "-file";

  my $outfile = basename($infile) . ".goodblat.tab";

  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			    );

  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   $F_OUT_GOODBLAT,
					  ],
			      "-clobber" => 1,
			      "-auto_qc" => 1,
			     );
  my @queue;
  my $flush = sub {
    if (@queue) {
      add_batch_goodblat(
			"-rows" => \@queue,
		       );
      foreach my $row (@queue) {
	my $usable = 1;
	$usable = 0 if $FLAGS{"good-only"} and not($row->{$F_OUT_GOODBLAT});
	$rpt->end_row($row) if $usable;
      }
    }
    @queue = ();
  };

  while (my $row = $df->get_hash()) {
    foreach my $f ($F_RA, $F_VA) {
      my $v = $row->{$f};
      $row->{$f} = "-" if $v eq "" or $v eq " ";
    }

    push @queue, $row;
    &$flush() if @queue >= $QUEUE_SIZE;
  }
  &$flush();
  $rpt->finish();
}

sub get_tabix_goodblat {
  my $f_tabix = $FLAGS{"tabix-goodblat"} || die "-tabix-goodblat";
  printf STDERR "GoodBad: %s\n", $f_tabix;
  return get_batch_tabix($f_tabix);
}

sub get_batch_tabix {
  my ($f_tabix) = @_;
  return new TabixFile(
		       "-file" => $f_tabix,
		       "-indel_wiggle_bases" => $TABIX_INDEL_EQUIVALENCE_DISTANCE
		      );
}


sub add_batch_goodblat {
  my (%options) = @_;
  my $rows = $options{"-rows"} || die "-rows";
  my $start_time = time;
  my $label = "GoodBLAT";

  log_message("batch $label start");
  my $tabix = get_tabix_goodblat();
  my $f_row_key = "_user_row";

  my @query_variants;
  foreach my $r (@{$rows}) {
    my $v = get_variant_from_row($r);
    $v->{$f_row_key} = $r;
    push @query_variants, $v;
  }

  my $tba = new TabixBatchAnnotation(
				     "-tabix" => $tabix,
				     "-twobit" => ($CONFIG_GENOME->{TWOBIT} || die),
				     "-split_count" => TABIX_BATCH_SIZE,
				     "-f_tabix_chr" => $F_TABIX_CHR,
				     "-f_tabix_pos" => $F_TABIX_POS,
				     "-f_tabix_ref_allele" => $F_TABIX_RA,
				     "-f_tabix_var_allele" => $F_TABIX_VA,

				     "-user_row_key" => $f_row_key,

				     "-store_hits" => $FIELD_GOODBLAT_TABIX,
#				     "-store_site" => $FIELD_CLINVAR_TABIX_SITE,
				    );

  $tba->query(
	      "-query" => \@query_variants,
	     );

  foreach my $row (@{$rows}) {
    my $v_goodbad = "";

    if (my $hits = $row->{$FIELD_GOODBLAT_TABIX}) {
      dump_die($row, "multiple hits") if @{$hits} > 1;

      my $usable = 1;
      if ($usable) {
	$v_goodbad = $hits->[0]->{$F_TABIX_GOODBLAT};
	die unless defined $v_goodbad;
      }
    }

    $row->{$F_OUT_GOODBLAT} = $v_goodbad;
  }

  log_message(sprintf "batch %s annotation for %d rows took %d",
	  $label, scalar(@{$rows}), time - $start_time);


}


sub get_variant_from_row {
  my ($row) = @_;
  my $chr = $row->{$F_CHR} || dump_die($row, "no $F_CHR");
  my $pos = $row->{$F_POS} || dump_die($row, "no $F_POS");
  my $ra = $row->{$F_RA};
  my $va = $row->{$F_VA};
  foreach ($ra, $va) {
    dump_die($row, "undef reference or variant allele") unless defined $_;
    $_ = "-" unless $_;
  }

  my $v = new Variant();
  $v->import_generic(
		     "-reference-name" => $chr,
		     "-base-number" => $pos,
		     "-reference-allele" => $ra,
		     "-variant-allele" => $va
		    );
  # don't need to bother with insertion position tweaking since
  # goodbad list is SNV only
  return $v;
}

