#!/bin/env perl
# prepare UCSC annotation flatfile for use with tabix
# MNE 6/2015
#
# TO DO:
# - embed genome build into filename?  wget path?
# - full auto tartan mode

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list);
use FileUtils qw(universal_open find_binary);
use TemporaryFileWrangler;
use TARTANUtils qw(tartan_genome_put_helper);
use TabixPrep;
use GenomeUtils qw(reverse_complement cook_chromosome_name);

my $F_TABIX_CHROM = "Chr";
my $F_TABIX_START = "WU_HG19_Pos";
# use the same field as vcf2tab.pl in case e.g. hg38 version
# is generated from NCBI VCF distribution

my $ENABLE_EXTENDED_FIELDS = 1;

my %FLAGS;
GetOptions(
	   \%FLAGS,

	   "-snp=s",
	   # UCSC genome annotation database snpXXX flatfile,
	   # e.g. snp142.txt.gz

	   "-headers=s",
	   # manually constructed header line file,
	   # 1 label per column (see table browser):
	   # http://genome.ucsc.edu/cgi-bin/hgTables?command=start

	   "-tartan-index=s",
	   # output dir
	   "-genome=s",

	   "-cook",
	   # digest/QC into simpler batch-tabix compatible format
	   "-max=i",
	   "-ping=i",
	  );

die "obsolete, use NCBI VCF distribution instead";

if (my $out_dir = $FLAGS{"tartan-index"}) {
  tartan_genome_put_helper(
			   "-genome" => $FLAGS{genome},
			   "-subdir" => "dbSNP_tabix",
			   "-out" => $out_dir,
			  );
  exit(0);
}


find_binary("tabix", "-die" => 1);
my $f_snp = $FLAGS{snp} || die "-snp";
my $f_headers = $FLAGS{headers};
my $tfw = new TemporaryFileWrangler(
#				   "-verbose" => 1
				   );
unless ($f_headers) {
  my $hl = <DATA>;
  chomp $hl;
  my @f = split /\s+/, $hl;
  $f_headers = $tfw->get_tempfile("-append" => ".headers.txt");
  open(HTMP, ">" . $f_headers) || die;
  printf HTMP "%s\n", join "\t", @f;
  close HTMP;
}

# QC: check header file, compare column counts
open(TMP, $f_headers) || die;
my $hl = 0;
my @headers;
while (<TMP>) {
  die "ERROR: header line must begin with # for compatibility w/tabix" unless /^\#/;
  die "ERROR: header file must contain exactly one line" if ++$hl > 1;
  chomp;
  @headers = split /\t/, $_;
}

my $fh = universal_open($f_snp) || die;
my $line = <$fh>;
my @f = split /\t/, $line;

die sprintf "ERROR: data file contains %d columns, header contains %d columns",
scalar @f, scalar @headers unless @f == @headers;

if ($FLAGS{cook}) {
  build_cooked_dbsnp();
  exit(0);
} else {
  die "raw tabix build obsolete, specify -cook option";
}

my %h2i;
my $i = 0;
%h2i = map {$_ => ++$i} @headers;
#dump_die(\%h2i);

my $i_chrom = $h2i{chrom} || die "no header for chrom";
my $i_start = $h2i{chromStart} || die "no header for chromStart";
my $i_end = $h2i{chromEnd} || die "no header for chromEnd";

my $outfile = basename($f_snp) . ".headered.gz";
my $outfile_tmp = $outfile . ".tmp";

my $cmd = sprintf 'zcat %s %s --stdout --force | bgzip > %s', $f_headers, $f_snp, $outfile_tmp;
run_command($cmd);
rename($outfile_tmp, $outfile) || die;

$cmd = sprintf 'tabix -s %d -b %d -e %d %s',
  $i_chrom, $i_start, $i_end, $outfile;
# positions are interbase, so tabix "0-based" option isn't quite right.
# just use defaults, which for a 1-based query will still match
# because the end coordinate is essentially 1-based.

run_command($cmd);

sub run_command {
  my ($cmd) = @_;
  printf STDERR "%s: running %s...\n", scalar(localtime), $cmd;
  system $cmd;
  die "$cmd exited with $?" if $?;
}

sub build_cooked_dbsnp {
  my $fh = universal_open($f_snp) || die;

  my $outfile = basename($f_snp);
  $outfile =~ s/\.gz$//;
  $outfile .= ".cooked.gz";

  my @h_out = (
		 $F_TABIX_CHROM,
		 $F_TABIX_START,
		 qw(
		     ReferenceAllele
		     MutantAllele
		     rs
		  )
		);

  if ($ENABLE_EXTENDED_FIELDS) {
    push @h_out, qw(
		     strand
		     class
		     exceptions
		     );
  }

  my $tp = new TabixPrep(
			 "-outfile" => $outfile,
			 "-headers" => \@h_out,
			 "-header_chr" => $F_TABIX_CHROM,
			 "-header_start" => $F_TABIX_START,
			 "-header_end" => $F_TABIX_START,
			);

  my %reject;
  my $rows = 0;
  my $max = $FLAGS{max};
  my $rejected = 0;
  my $ping = $FLAGS{ping};
  while (my $line = <$fh>) {
    chomp $line;

    $rows++;
    last if $max and $rows > $max;

    my %ucsc;
    @ucsc{@headers} = split /\t/, $line;
#    dump_die(\%ucsc, "start", 1);

    my $chr = $ucsc{chrom};

    my $base_ref = $ucsc{refNCBI};
    foreach my $base_variant (get_variant_bases(\%ucsc)) {
      dump_die(\%ucsc, "variant base is ref") if $base_variant eq $base_ref;

      my $start = $ucsc{chromStart};
      # UCSC interbase

      printf STDERR "%s.%d: %d\n", $chr, $start, $rows if $ping and $rows % $ping == 0;

      $start++ unless $base_ref eq "-";
      # leave insertions the base # before the inserted sequence,
      # for compatibility with batch tabix.  All other positions
      # converted to in-base / 1-based.

      my %r;
      $r{$F_TABIX_CHROM} = $chr;
      $r{$F_TABIX_START} = $start;
      $r{ReferenceAllele} = $base_ref;
      $r{MutantAllele} = $base_variant;
      $r{rs} = $ucsc{name} || die;

      # debug options, maybe take these out eventually:
      $r{class} = $ucsc{class};
      $r{exceptions} = $ucsc{exceptions};
      $r{strand} = $ucsc{strand};

      my $usable = 1;

      unless (cook_chromosome_name($chr, "-undef-unknown" => 1)) {
	$reject{"noncanonical_chromosome__" . $chr}++;
	$usable = 0;
      }

      if ($base_variant eq "LENGTHTOOLONG") {
	$reject{"LENGTHTOOLONG"}++;
	$usable = 0;
      }

      if ($usable) {
	$tp->add_row("-row" => \%r);
      } else {
	$rejected++;
      }
    }
  }
  $tp->finish();

  if ($rejected) {
    printf STDERR "rejected %d variants:\n", $rejected;
    foreach my $k (sort keys %reject) {
      printf STDERR "  %s: %d\n", $k, $reject{$k};
    }
  }

}

sub get_variant_bases {
  # return genomic variant bases for a SNV (correcting for strand)
  my ($hit) = @_;
  my $observed = $hit->{observed} || die;
  my $strand = $hit->{strand} || die;
  my @entries = split /\//, $observed;
  my $reference_base = uc($hit->{refNCBI} || die);
  my $saw_reference;
  my @variant_bases;
  foreach my $entry (@entries) {
    if ($strand eq "-") {
      $entry = reverse_complement($entry);
    } elsif ($strand ne "+") {
      die;
    }

    $entry = uc($entry);
    if ($entry eq $reference_base) {
      $saw_reference = 1;
    } else {
      push @variant_bases, $entry;
    }
  }

#  rs144571919

#  dump_die($hit, sprintf("WARNING: didn't find reference %s in %s (%s)!",
#			 $reference_base, $observed, $hit->{name}), 1)
#      unless $saw_reference;
  # turn off: happens frequently for ?complex? indels

  # rs144571919:
  # reference is a C, observed bases are A/C, strand is +...WTF?
  return @variant_bases;
}


# headers retrieved from http://genome.ucsc.edu/cgi-bin/hgTables (snp* table)
__DATA__
#bin    chrom   chromStart      chromEnd        name    score   strand  refNCBI refUCSC observed        molType class   valid   avHet   avHetSE func    locType weight  exceptions      submitterCount  submitters      alleleFreqCount alleles alleleNs        alleleFreqs     bitfields
