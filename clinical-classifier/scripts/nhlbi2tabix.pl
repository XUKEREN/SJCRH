#!/bin/env perl
# convert NHLBI data to tabix format, for either genome present in the file
# OBSOLETE: see nhlbi_vcf2tab.pl which pre-cooks output
# MNE 8/2015

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use File::Basename;
use TabixPrep;

use MiscUtils qw(dump_die build_argv_list);
use GenomeUtils qw(cook_chromosome_name);
use TARTANUtils qw(tartan_genome_put_helper);

use constant LIFTOVER_FAILURE_CODE => -1;

my $INCLUDE_INDELS = 0;

my %FLAGS;
my @clopts = (
	      "-genome=s",
	      "-nhlbi-dir=s",

	      "-max-rows=i",
	      "-glob-filter=s",
	      "-include-indels=i" => \$INCLUDE_INDELS,

	      "-tartan-index=s",

	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

die "OBSOLETE: use nhlbi_vcf2tab.pl";

if (my $out_dir = $FLAGS{"tartan-index"}) {
  tartan_genome_put_helper(
			   "-genome" => $FLAGS{genome},
			   "-subdir" => "NHLBI_tabix",
			   "-out" => $out_dir,
			  );
  exit(0);
}

printf STDERR "including indels: %s\n", $INCLUDE_INDELS ? "yes" : "no";

my $dir = $FLAGS{"nhlbi-dir"} || die "-nhlbi-dir";
my $genome = $FLAGS{"genome"} || die "-genome";

my @nhlbi_files = glob($dir . "/*.txt");
my $glob_filter = $FLAGS{"glob-filter"};
@nhlbi_files = grep {/$glob_filter/} @nhlbi_files if $glob_filter;
die "no files" unless @nhlbi_files;

my $rpt;

my $F_POSITION;
if ($genome eq "GRCh37-lite") {
  $F_POSITION = "#base(NCBI.37)";
} elsif ($genome eq "GRCh38") {
  $F_POSITION = "GRCh38Position";
} else {
  die "position column not defined for genome $genome";
}
die "no position column" unless $F_POSITION;
my $outfile = sprintf "NHLBI_tabix_%s.tab.gz", $genome;

my $OUT_F_CHROM = "Chromosome";
my $OUT_F_POSITION = "Position";

my $max_rows = $FLAGS{"max-rows"};

my $tp;

my %unusable;
foreach my $file (@nhlbi_files) {
  printf STDERR "scanning %s...\n", $file;
  my $headers_raw = get_headers($file);
  my @headers_out = (@{$headers_raw}, $OUT_F_CHROM, $OUT_F_POSITION);

  unless ($tp) {
    $tp = new TabixPrep(
			"-headers" => \@headers_out,
			"-outfile" => $outfile,
			"-header_chr" => $OUT_F_CHROM,
			"-header_start" => $OUT_F_POSITION
		       );
  }

  open(IN, $file) || die;
  my $rows = 0;
  while (my $line = <IN>) {
    next if $line =~ /^#/;
    chomp $line;
    if ($max_rows and ++$rows > $max_rows) {
      print STDERR "debug, stopping after $max_rows\n";
      last;
    }
    my @f = split /\s+/, $line;
    die sprintf "fields:%d headers:%d in file $file line $line", scalar(@f), scalar(@{$headers_raw}) unless @f == scalar @{$headers_raw};
    my %r;
    @r{@{$headers_raw}} = @f;
#    dump_die(\%r, "header mismatch") if $check;
#    dump_die(\%r);

    my $pos_raw = $r{$F_POSITION} || die "no field $F_POSITION";
    if ($pos_raw eq LIFTOVER_FAILURE_CODE) {
      $unusable{bad_liftover}++;
    } else {
      @f = split /:/, $pos_raw;
      dump_die(\%r, "invalid position $pos_raw") unless @f == 2;
      my ($chrom, $pos) = @f;
      if (cook_chromosome_name($chrom, "-undef-unknown" => 1)) {

	my $usable = 1;
	my $alleles = $r{Alleles};
	my $has_indel;
	if ($alleles =~ />/) {
	  # newer release with format e.g. A>T;A>C
	  my @alleles = split /;/, $alleles;
	  foreach my $a (@alleles) {
	    my @f = split />/, $a;
	    die unless @f == 2;
	    $has_indel = 1 if length($f[0]) != length($f[1]);
	  }
	}
	if (not($INCLUDE_INDELS) and $has_indel) {
	  $unusable{indel}++;
	} else {
	  $r{$OUT_F_CHROM} = $chrom;
	  $r{$OUT_F_POSITION} = $pos;
	  $tp->add_row("-row" => \%r);
	}
      } else {
	$unusable{non_canonical_chrom}++;
      }
    }
  }
}
$tp->finish();

if (%unusable) {
  printf STDERR "unusable records summary:\n";
  foreach (sort keys %unusable) {
    printf STDERR "  %s: %d\n", $_, $unusable{$_};
  }
}

sub get_headers {
  my ($file) = @_;
  open(IN, $file) || die;
  my $header_line;
  while (<IN>) {
    chomp;
    if (/^\#/) {
      # header section
      if (/^##/) {
	# extra lines, ignore
      } else {
	die if $header_line;
	$header_line = $_;
      }
    } else {
      last;
    }
  }
  close IN;
  die unless $header_line;
  my @headers = split /\s+/, $header_line;
  return \@headers;
}

