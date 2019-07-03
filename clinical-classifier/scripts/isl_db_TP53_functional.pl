#!/bin/env perl
# TP53 functional data
# TO DO: pubmed references for individual subsets?

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);

use MiscUtils qw(dump_die build_argv_list get_hash_option median);
use FileUtils qw(find_binary);
use DelimitedFile;
use Reporter;
use DBI;
use SimpleVariantDB qw(
			F_VDB_CHR
			F_VDB_POS
			F_VDB_RA
			F_VDB_VA
			F_VDB_GENE
			F_VDB_TRANSCRIPT
			F_VDB_AA
			F_VDB_PMID
		     );

use constant DB_NAME => "TP53_functional";

use constant CALL_FUNCTIONAL => "F";
use constant CALL_NON_FUNCTIONAL => "N";
use constant CALL_UNKNOWN => "U";

use constant THRESHOLD_FUNCTIONAL => 0.25;
use constant THRESHOLD_NON_FUNCTIONAL => 0.75;

my %FLAGS;
my @clopts = (
	      "-prep=s",
	      # parse/extract sub-columns before loading

	      "-out=s",

	      "-load=s",
	      # load generated file into database
	      "-db-file=s",
	      "-patch=s",

	      "-cache",
	     );

my @ARGV_RAW = @ARGV;
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

prep_data() if $FLAGS{prep};
load_db() if $FLAGS{load};

sub load_db {
  my $f_db = $FLAGS{"db-file"} || die "-db-file";
  my $dbi = DBI->connect("dbi:SQLite:dbname=$f_db","","");
  my $infile = $FLAGS{load} || die "-load";

  my %map;
  $map{F_VDB_CHR()} = "sj_chr";
  # actually a constant
  $map{F_VDB_POS()} = "sj_pos";
  $map{F_VDB_RA()} = "sj_ref_allele";
  $map{F_VDB_VA()} = "sj_alt_allele";

  foreach my $f_flat (qw(
		   WAF1_Act.patched
		   MDM2_Act.patched
		   BAX_Act.patched
		   _14_3_3_s_Act.patched
		   AIP_Act.patched
		   GADD45_Act.patched
		   NOXA_Act.patched
		   p53R2_Act.patched
		   Giac_A549_WT_Nut_norm
		   Giac_A549_Null_Nut_norm
		   Giac_A549_Null_Eto_norm.inverted
		   Kolt_RFS_H1299_norm

		   functional_call
		  )) {
    my $f_db = $f_flat;
    $f_db =~ s/\..*$//;
    # periods not allowed in SQL table names

    $map{$f_db} = $f_flat;
  }

  my %constants;
  $constants{F_VDB_GENE()} = "TP53";

  my %properties;

  my $svdb = new SimpleVariantDB(
				 "-dbi" => $dbi,
				);

  $svdb->init_table(
		    "-name" => DB_NAME,
		    "-file" => $infile,
		    "-map" => \%map,
		    "-map-only" => 1,
		    "-auto-configure" => 1,
		    "-load" => 1,
		    "-properties" => \%properties,
		    "-constants" => \%constants,
		   );

}

sub prep_data {
  # step 1:
  # run hgvs_split.pl on Excel-covered file to extract
  # annotations in discrete columns
  my $f_in = $FLAGS{"prep"} || die "-prep";
  my $f_patch = $FLAGS{patch} || die "-patch";

  find_binary("hgvs_split.pl", "-die" => 1);

  my $f_hgvs = $f_in . ".hgvs_split.tab";
  if (($FLAGS{cache} and -s $f_hgvs) ? 0 : 1) {
    # step 1: extract HGVS fields and correct alleles
    my $cmd = sprintf 'hgvs_split.pl -field-hgvs HG19_Variant -file %s -split -extract-gpos -extract-hgvs-chr -rc -genome GRCh37-lite', $f_in;
    # FIX ME: other genomes
    print STDERR "$cmd\n";
    system $cmd;
  } else {
    printf STDERR "skipping HGVS conversion, exists\n";
  }

  # step 2:
  # patch in normalized TP53 values from second spreadsheet
  die unless -s $f_hgvs;

  my $df = new DelimitedFile("-file" => $f_patch,
			     "-headers" => 1,
			     );
  my %patch;
  while (my $row = $df->get_hash()) {
    my $aa = $row->{AA} || die;
    # manual edit of 1st column name
    die if $patch{$aa};
    $patch{$aa} = $row;
  }

  $df = new DelimitedFile("-file" => $f_hgvs,
			  "-headers" => 1,
			 );
  my $f_patched = basename($f_hgvs) . ".patched.tab";

  my @f_orig = qw(
		   WAF1_Act
		   MDM2_Act
		   BAX_Act
		   _14_3_3_s_Act
		   AIP_Act
		   GADD45_Act
		   NOXA_Act
		   p53R2_Act
		);

  my $rpt = $df->get_reporter(
			      "-file" => $f_patched,
			      "-extra" => [
					   (map {$_ . ".patched"} @f_orig),
					   "normalization_unpatchable",
					   "Giac_A549_Null_Eto_norm.inverted",
					   "functional_call",
					  ],
  			      "-auto_qc" => 1,
			     );

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    my $aa = $row->{Protein_p1_TP53} || die;
    $aa =~ s/^p\.//;

    my $normalization_unpatchable = "";

    foreach my $f (@f_orig) {
      my $v_orig = $row->{$f};
      die unless defined $v_orig;
      my $v_patched = $v_orig;

      my $f_patched = $f;
      $f_patched =~ s/_Act$//;

      $f_patched =~ s/^_//;
      $f_patched = "AIP1" if $f_patched eq "AIP";
      # hacktacular

#      print STDERR "orig $v_orig\n";

      #
      # patch promoter activity with normalized data:
      #
      if ($v_orig eq "No_Protein" or
	  $v_orig eq "No_Data" or
	  $v_orig eq "Not_Relevant") {
	# ignore
      } elsif ($aa =~ /^[A-Z]\d+[A-Z]+$/) {
	# AA substitution
	my $patch_row = $patch{$aa} || die "no match for $aa";
	dump_die($patch_row, "no column for aa=$aa f=$f_patched") unless exists $patch_row->{$f_patched};
	$v_patched = $patch_row->{$f_patched};
	$v_patched = $v_orig if $v_patched eq "";
	# use original value if blank, e.g. P89H
	if ($v_patched =~ /E\-/) {
	  my $new = sprintf '%f', $v_patched;
#	  die "$v_patched $new";
	  $v_patched = $new;
	  # convert exponential notation to float
	}
	dump_die($patch_row, "no value for $f_patched") unless defined $v_patched;
      } else {
	# other AA value
	$normalization_unpatchable = $aa;
      }

#      printf STDERR "$f_patched $v_patched\n";
      $row->{$f . ".patched"} = $v_patched;
    }
    $row->{normalization_unpatchable} = $normalization_unpatchable;

    #
    # invert values for GiaC experiment 3:
    # "High value define as functional and Low value as non functional for experiment 3 (Giac_A549_Null_Eto_norm)"
    #
    my $raw = clean_value($row->{Giac_A549_Null_Eto_norm});
    my $inverted = $raw;
    if ($raw =~ /^[\d\.]+$/) {
      $inverted = 1 - $raw;
    } elsif ($raw eq "No_Data") {
      # ignore
    } else {
      die "unhandled value $raw";
    }
    $row->{"Giac_A549_Null_Eto_norm.inverted"} = $inverted;

    #
    #  generate digested functional call:
    #
    my %calls;
    foreach my $field_set (
			   [ map {$_ . ".patched"} @f_orig ],
			   # promoter activitiy

			   [ "Giac_A549_WT_Nut_norm",
			     "Giac_A549_Null_Nut_norm",
			     "Giac_A549_Null_Eto_norm.inverted" ],
			   # Giacomelli data, using inverted values for
			   # experiment 3

			   [ "Kolt_RFS_H1299_norm" ],
			   # Kolter data
			  ) {

      my @raw = map {clean_value($row->{$_})} @{$field_set};
      my @usable;
      foreach my $v (@raw) {
	if ($v eq "No_Protein" or $v eq "No_Data" or $v eq "Not_Relevant" or $v eq "-") {
	} elsif ($v =~ /^[\d\.]+$/) {
	  push @usable, $v;
	} else {
	  dump_die($row, "unhandled field $v");
	  die "unhandled value \"$v\"";
	}
      }

      my $call;
      if (@usable) {
	my $median = median(\@usable);
	if ($median <= THRESHOLD_FUNCTIONAL) {
	  $call = CALL_FUNCTIONAL;
	} elsif ($median <= THRESHOLD_NON_FUNCTIONAL) {
	  $call = CALL_UNKNOWN;
	} else {
	  $call = CALL_NON_FUNCTIONAL;
	}
#      } else {
#	$call = CALL_UNKNOWN;
	# rather than calling unknown, just ignore because sometimes
	# datasets aren't populated
	# FIX ME: maybe don't use the data in this case?
      }
      $calls{$call}++ if $call;
    }

    my $final_call;
    if (scalar keys %calls == 0) {
      # e.g. no sets have data
      $final_call = CALL_UNKNOWN;
    } elsif (scalar keys %calls == 1) {
      ($final_call) = keys %calls;
    } else {
      # two or more calls observed.
      ($final_call) = sort {$calls{$b} <=> $calls{$a}} keys %calls;
      my $final_call_count = $calls{$final_call};
      my @equivalent = grep {$_ == $final_call_count} values %calls;
      if (@equivalent > 1) {
	# call does not have dominant value
	dump_die(\%calls, "tie of " . scalar @equivalent, 1);
	if (0 and @equivalent == 3) {
	  # 3 way tie
	  $final_call = CALL_UNKNOWN;
	  # OK??
	} else {
	  my %tie;
	  foreach my $c (keys %calls) {
	    if ($calls{$c} == $final_call_count) {
	      $tie{$c} = 1;
	    }
	  }

	  my %call2rank;
	  $call2rank{CALL_FUNCTIONAL()} = 1;
	  $call2rank{CALL_UNKNOWN()} = 2;
	  $call2rank{CALL_NON_FUNCTIONAL()} = 3;

	  my ($winner) = sort {$call2rank{$a} <=> $call2rank{$b}} keys %tie;

	  print STDERR "winner=$winner\n";
	  $final_call = $winner;
	}
      }
    }
    die unless $final_call;

    $row->{functional_call} = $final_call;

    $rpt->end_row($row);
  }
  $rpt->finish();
}

sub clean_value {
  # unquote and convert European-style "," decimal to US-style "."
  my ($v) = @_;
  $v =~ s/^\"(.*)\"$/$1/;
  $v =~ s/,/./;
  $v =~ s/\s+$//;
  return $v;
}
