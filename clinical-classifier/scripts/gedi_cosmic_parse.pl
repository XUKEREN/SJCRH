#!/bin/env perl
# split CDS/AA annotations into component elements
# MNE 6/2014

# BEWARE:
# - bases may need strand correction! (COSMIC)
# - CDS annotations using +/- intronic base annotation
# - standardization of insertion site base #

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die);
# use DelimitedFile;
# use Reporter;
use NucleotideSubstitutionParser;
use AAParser;
use Reporter;

my %FLAGS;
GetOptions(
	   \%FLAGS,
	   "-file=s",
	   # input file

	   "-out=s",
	   # output file (optional)

	   "-debug-aa=s",
	   "-debug-cds=s",
	  );

my $fn = $FLAGS{file} || die "-file";
open(IN, $fn) || die;

my $outfile = $FLAGS{out} || sprintf '%s.split.tab', basename($fn);

my $rpt = new Reporter(
		       "-file" => $outfile,
		       "-delimiter" => "\t",
		       "-labels" => [
				     qw(
					 field1
					 field2
					 cds_raw
					 aa_raw
					 cds_category
					 cds_start
					 cds_end
					 cds_coding_intronic
					 cds_ref_base
					 cds_var_base
					 aa_category
					 aa_codon_start
					 aa_codon_end
					 aa_codon_reference
					 aa_codon_variant
				      )
				    ],
			 "-auto_qc" => 1,
			);

my $wf = new WorkingFile($outfile);
my $fh = $wf->output_filehandle();

printf STDERR "WARNING: strand data needed for COSMIC\n";

my $debug_aa = $FLAGS{"debug-aa"};
my $debug_cds = $FLAGS{"debug-cds"};

while (<IN>) {
  chomp;
  my @f = split /\t/, $_;
  die "4 fields expected" unless @f;
  my ($field1, $field2, $cds_raw, $aa_raw) = @f;

#  print STDERR "DEBUG FILTER\n"; next unless $cds_raw eq "c.603+1G>A";
#  print STDERR "DEBUG FILTER\n"; next unless $cds_raw eq "c.3123+1G>GA";
#  print STDERR "DEBUG FILTER\n"; next unless $cds_raw eq "c.2203-2C>A";
#  print STDERR "DEBUG FILTER\n"; next unless $cds_raw eq "c.655-12_655-9delaacc";
#  print STDERR "DEBUG FILTER\n"; next unless $cds_raw eq "c.375+1_375+18del18";
#  print STDERR "DEBUG FILTER\n"; next unless $cds_raw eq "c.1325_1335+12del23";
#  print STDERR "DEBUG FILTER\n"; next unless $aa_raw eq "p.L637_D638delLD";

  if ($debug_aa) {
    print STDERR "DEBUG FILTER AA $debug_aa\n";
    next unless $aa_raw eq $debug_aa;
  }

  if ($debug_cds) {
    print STDERR "DEBUG FILTER CDS $debug_cds\n";
    next unless $cds_raw eq $debug_cds;
  }

  #
  #  parse nucleotide annotation:
  #
  my $np = new NucleotideSubstitutionParser();
  # $np->auto_strand_fix($strand);
  # FIX ME: THIS MAY BE NEEDED!

  my $cds_category = "cds_parse_exception";
  my $cds_start = "";
  my $cds_end = "";
  my $cds_ref_base = "";
  my $cds_var_base = "";
  my $cds_coding_intronic = 0;

  if ($np->parse($cds_raw)) {
    $cds_coding_intronic = $np->is_coding_intronic;
    if ($np->is_substitution()) {
      # simple substitution
      $cds_category = "simple_substitution";
      $cds_start = $np->start();
      $cds_end = $np->end();
      $cds_ref_base = $np->reference_sequence();
      $cds_var_base = $np->variant_sequence();
    } elsif ($np->is_deletion()) {
      # simple deletion
      $cds_category = "simple_deletion";
      $cds_start = $np->start();
      $cds_end = $np->end();
      my $len = ($cds_end - $cds_start) + 1;
      unless ($np->event_length == $len) {
	printf STDERR "WARNING: inconsistency between location len %d and event size %d in %s\n", $len, $np->event_length, $cds_raw;
	# can't just die, data not always consistent, e.g. 7_378del1016
      }
      $cds_ref_base = $np->reference_sequence;
      if ($cds_ref_base and
#	  length($cds_ref_base) == $len and
# might be inconsistency if intronic event
	  $cds_ref_base =~ /^[ACGTN]+$/i) {
	# OK, keep
      } else {
	# unknown/broken: simulate
	$cds_ref_base = "N" x $len;
      }
      $cds_var_base = "-";
    } elsif ($np->is_insertion()) {
      # simple insertion
      $cds_category = "simple_insertion";
      $cds_start = $np->start();
      $cds_end = $np->end();
      $cds_ref_base = "-";
      my $elen = $np->event_length() || die "no elen for $cds_raw";

      $cds_var_base = $np->variant_sequence();
      unless ($cds_var_base) {
	# sometimes not provided, e.g. "c.3907_3908ins10"
	$cds_var_base = "N" x $elen;
      }
      die $cds_raw unless length($cds_var_base) == $elen;
    } elsif ($np->is_complex_indel()) {
      $cds_category = "complex_indel";
      $cds_start = $np->start();
      $cds_end = $np->end();
      $cds_ref_base = $np->reference_sequence();
      $cds_var_base = $np->variant_sequence();
    } else {
      print STDERR "not handling parsed event $cds_raw\n";
      die;
    }
  } else {
    printf STDERR "WARNING: skipping unparsable COSMIC event %s\n", $cds_raw;
  }

  #
  #  parse AA annotation:
  #
  my $aap = new AAParser();

  my $aa_category = "aa_parse_exception";
  my $is_substitution = "";
  my $aa_codon_start = "";
  my $aa_codon_end = "";
  my $aa_codon_reference = "";
  my $aa_codon_variant = "";

  if ($aap->parse_substitution($aa_raw)) {
    $aa_category = "aa_substitution";
    $aa_codon_start = $aa_codon_end = $aap->codon_number;
    $aa_codon_reference = $aap->codon_reference;
    $aa_codon_variant = $aap->codon_variant;
  } elsif ($aap->parse($aa_raw)) {
    $aa_category = "aa_other";
    $aa_codon_start = $aap->codon_start;
    $aa_codon_end = $aap->codon_end;
    $aa_codon_reference = $aap->codon_reference;
    $aa_codon_variant = $aap->codon_variant;
  } else {
#    printf STDERR "unparsable AA %s\n", $aa_raw;
  }

  my %r;
  $r{field1} = $field1;
  $r{field2} = $field2;
  $r{cds_raw} = $cds_raw;
  $r{aa_raw} = $aa_raw;

  $r{cds_category} = $cds_category;
  $r{cds_start} = $cds_start;
  $r{cds_end} = $cds_end;
  $r{cds_coding_intronic} = $cds_coding_intronic;
  $r{cds_ref_base} = $cds_ref_base;
  $r{cds_var_base} = $cds_var_base;
  $r{aa_category} = $aa_category;
  $r{aa_codon_start} = $aa_codon_start;
  $r{aa_codon_end} = $aa_codon_end;
  $r{aa_codon_reference} = $aa_codon_reference;
  $r{aa_codon_variant} = $aa_codon_variant;

  $rpt->end_row(\%r);
}

$rpt->finish();
