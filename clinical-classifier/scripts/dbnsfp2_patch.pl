#!/bin/env perl
# patch dbNSFP v2 with a few annotations from dbNSFP v3.5a so we can
# hold off on full parsing upgrade for now.
# WARNING: matching is sometimes ambiguous!
# MNE 2/2018

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
use TdtConfig;
use DelimitedFileHP;
use TabixPrep;

my $CACHE_SIZE = 100000;
my $F_ENSEMBL_TRANSCRIPT = "Ensembl_transcriptid";
use constant DBNSFP_NULL => ".";

my %FLAGS;
my @clopts = (
	      "-genome=s",
	      "-max=i",

	      "-dbnsfp2=s",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my @F_CADD = qw(
		 CADD_raw
		 CADD_phred
	      );
my @F_REVEL = qw(
		 REVEL_score
	      );
# exclude "rankscore" fields for now to save space

my @F_COPY = (@F_CADD, @F_REVEL);
# in some cases we have to drill down to subsets of fields to
# help resolve ambiguity

my $genome = $FLAGS{genome} || die;
die "untested genome" unless $genome eq "GRCh37-lite";

my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
#dump_die($config_genome);

my $f_nsfp2 = $FLAGS{dbnsfp2};
unless ($f_nsfp2) {
  ($f_nsfp2) = glob($config_genome->{DBNSFP_TABIX_DIR} . "/*.gz");
}

my ($f_nsfp3) = glob($config_genome->{DBNSFP3_TABIX_DIR} . "/*.gz");

printf STDERR "nsfp2: %s\n", $f_nsfp2;
printf STDERR "nsfp3: %s\n", $f_nsfp3;

my $F_CHR2 = "#chr";
my ($F_POS2) = $config_genome->{DBNSFP_TABIX_POS_FIELD} || die;
my $F_CHR3 = "hg19_chr";
my $F_POS3 = "hg19_pos(1-based)";

my $F_REF = "ref";
my $F_ALT = "alt";

my $f_tabix_out = basename($f_nsfp2) . ".patched.gz";

my $df = new DelimitedFileHP(
			     "-file" => $f_nsfp2,
			     "-headers_extra" => \@F_COPY,
			    );
$df->write_init(
		"-file" => $f_tabix_out,
		"-compress" => "bgzip",
	       );

my ($last_chr, $query_start, $query_end);
my %cache;

my $max = $FLAGS{max};
my $rows_read=0;

my %counts;

while ($df->next_row()) {
  my $this_chr = $df->get_value($F_CHR2);
  my $this_pos = $df->get_value($F_POS2);

  my $need_query;
  $need_query = 1 unless $this_chr eq ($last_chr || "");
  $need_query = 1 unless $this_pos <= ($query_end || 0);
  $last_chr = $this_chr;

  if ($need_query) {
    %cache = ();
    my $loaded = 0;
    $query_end = $this_pos + $CACHE_SIZE;
    my $cmd = sprintf 'tabix -h %s %s:%d-%d', $f_nsfp3, $this_chr, $this_pos, $query_end;
    my $fh = new FileHandle();
    $fh->open($cmd . "|");

    my $df3 = new DelimitedFileHP(
				  "-fh_in" => $fh
				 );
    while ($df3->next_row) {
      my $key = get_key($df3, $F_CHR3, $F_POS3);
      my %excerpt = map {$_, $df3->get_value($_)} @F_COPY, $F_ENSEMBL_TRANSCRIPT;
      push @{$cache{$key}}, \%excerpt;
      $loaded++;
    }

    printf STDERR "query=%s loaded=%d\n", $cmd, $loaded;
  }

  my $key = get_key($df, $F_CHR2, $F_POS2);
  my %extra;

  if (my $hits = $cache{$key}) {

    if (@{$hits} == 1) {
      $counts{exactly_one_raw}++;
    }

    if (@{$hits} > 1) {
      # 1.865715.A.T: multiple rows, but all columns identical
      my $identical = identical_check($hits, \@F_COPY);
      if ($identical) {
	# since all rows identical, just use the first one
	$hits = [ $hits->[0] ];
	$counts{rescue_identical}++;
      }
    }

    if (@{$hits} > 1) {
      my %tid_2 = map {$_, 1} split /;/, $df->get_value($F_ENSEMBL_TRANSCRIPT);
      my $wanted;
      foreach my $hit (@{$hits}) {
	foreach my $tid_3 (split /;/, $hit->{$F_ENSEMBL_TRANSCRIPT}) {
	  if ($tid_2{$tid_3}) {
	    # able to identify the specific matching row via
	    # ENSEMBL transcript ID
	    #	      dump_die($df->get_hash, "debug 2 $tid_3", 1);
	    #	      dump_die($hit, "debug 3 $tid_3");
	    #	      printf STDERR "%s: disambiguate via %s\n", $key, $tid_3;
	    $wanted = $hit;
	    last;
	  }
	}
      }
      if ($wanted) {
	# disard non-matching rows
	$hits = [ $wanted ];
	$counts{rescue_transcript_id}++;
      }
    }

    if (@{$hits} == 0) {
      die "unpossible";
    } elsif (@{$hits} == 1) {
      # 1:1 row mapping, simply copy all fields
      foreach my $f (@F_COPY) {
	$extra{$f} = $hits->[0]->{$f};
      }
#      $counts{exactly_one_cooked}++;
    } else {
      #
      #  we still have more than one possible matching row in dbNSFP v3.
      #  might be able to resolve ambiguity within SUBSETS of the fields.
      #

      foreach my $f_subset (\@F_CADD, \@F_REVEL) {
	if (identical_check($hits, $f_subset)) {
	  # all values for these fields identical, so can use any row
	  foreach my $f (@{$f_subset}) {
	    $extra{$f} = $hits->[0]->{$f};
	  }
	  $counts{"subset_rescue_identical_" . $f_subset->[0]}++;
	} else {
	  my $set_non_null = get_set_non_null($hits, $f_subset);

	  if (@{$set_non_null} == 1) {
	    # only one row has non-null data for the set, so just use that
	    foreach my $f (@{$f_subset}) {
	      $extra{$f} = $set_non_null->[0]->{$f};
	    }
	    $counts{"subset_rescue_non_null_" . $f_subset->[0]}++;
	  } else {
	    # can't resolve ambiguity; just use data from the first record.
	    my $set = @{$set_non_null} ? $set_non_null : $hits;
	    # prefer non-null rows, if any
	    foreach my $f (@{$f_subset}) {
	      $extra{$f} = $set->[0]->{$f};
	    }
	    $counts{"subset_ambig_" . $f_subset->[0]}++;
	  }
	}
      }

    }

#     foreach my $f (@F_COPY) {
#       my $list = unique_ordered_list([grep {$_ ne "."} map {$_->{$f}} @{$hits}]);
#       # "." = dbNSFP null value
#       if (@{$list} > 1) {
# #	dump_die($df->get_hash(), "multi values for $f for $key: " . join ",", @{$list});
# 	$counts{ambiguous}++;
#       }
#       my $v = $list->[0];
#       # if ambiguous, just pick the first one.
#       $v = "" unless defined $v;
#       $extra{$f} = $v;
#     }
  } else {
    # no matches
    foreach my $f (@F_COPY) {
      $extra{$f} = "";
    }
    $counts{no_data}++;
  }

  $df->write_row("-extra" => \%extra);

#  dump_die($df->get_hash, $need_query);
#  my $uid = $df->get_value("MutationAssessor_UniprotID");

  $rows_read++;
  last if $max and $rows_read >= $max;
}

$df->write_finish();

printf STDERR "total:%d\n", $rows_read;
foreach my $f (sort keys %counts) {
  printf STDERR "  %s: %d\n", $f, $counts{$f};
}

my $tf = new TabixFile(
		       "-file" => $f_tabix_out,
		       "-index" => 1,
		       "-f_chr" => $F_CHR2,
		       "-f_start" => $F_POS2,
		       "-f_end" => $F_POS2,
		      );

sub get_key {
  my ($df, $f_chr, $f_pos) = @_;
  return join ".", map {$df->get_value($_)} $f_chr, $f_pos, $F_REF, $F_ALT;
}

sub identical_check {
  # are all values for a set of keys unique for a set of hashes?
  my ($hits, $fields) = @_;

  my $identical = 1;

  foreach my $f (@{$fields}) {
    my %unique;
    foreach my $h (@{$hits}) {
      $unique{$h->{$f}} = 1;
    }
    if (scalar keys %unique > 1) {
      $identical = 0;
      last;
    }
  }

  return $identical;
}


sub get_set_non_null {
  # are all values for a set of keys unique for a set of hashes?
  my ($hits, $fields) = @_;

  my @non_null;

  foreach my $h (@{$hits}) {
    my $has_anything;
    foreach my $f (@{$fields}) {
      my $v = $h->{$f};
      if ($v eq DBNSFP_NULL or $v eq "") {
      } else {
	$has_anything = 1;
	last;
      }
    }
    push @non_null, $h if $has_anything;
  }
  return \@non_null;
}
