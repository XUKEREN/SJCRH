#!/bin/env perl
# convert NHLBI VCF distribution to tabix format, precomputing
# genotype frequencies
# TO DO:
# - hg38
# - tartan
# MNE 6/2016

use Getopt::Long;
use File::Basename;
use TabixPrep;

use MiscUtils qw(dump_die build_argv_list);
use TARTANUtils qw(tartan_genome_put_helper);

my %FLAGS;
my @clopts = (
	      "-genome=s",
	      "-vcf-glob=s",

	      "-tartan-index=s",
	      # output dir
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );


if (my $out_dir = $FLAGS{"tartan-index"}) {
  tartan_genome_put_helper(
			   "-genome" => $FLAGS{genome},
			   "-subdir" => "NHLBI_tabix",
			   "-out" => $out_dir,
			  );
  exit(0);
}


my @vcfs = glob($FLAGS{"vcf-glob"} || die "-vcf-glob");
die "no VCFs" unless @vcfs;

my $F_TABIX_CHROM = "Chr";
my $F_TABIX_START = "Pos";
my $F_TABIX_RA = "ReferenceAllele";
my $F_TABIX_VA = "MutantAllele";

my @headers = (
	       $F_TABIX_CHROM,
	       $F_TABIX_START,
	       qw(
		   ReferenceAllele
		   MutantAllele
		   genotype_frequency
		)
	      );

my $tp = new TabixPrep(
		       "-outfile" => "nhlbi_cooked.tab.gz",
		       "-headers" => \@headers,
		       "-header_chr" => $F_TABIX_CHROM,
		       "-header_start" => $F_TABIX_START,
		      );

foreach my $vcf (@vcfs) {
  my $outfile = basename($vcf) . ".cooked.tab";
  unless (-s $outfile) {
    my $cmd = sprintf 'vcf2tab.pl -nhlbi -file %s -no-insertion-adjustment', $vcf;
    printf STDERR "%s\n", $cmd;
    system $cmd;
    die "$cmd exited with $?" if $?;
  }

  my $df = new DelimitedFile("-file" => $outfile,
			     "-headers" => 1,
			     );
  while (my $row = $df->get_hash()) {
    my %r = %{$row};
    $r{$F_TABIX_START} = $row->{Pos};
    $r{$F_TABIX_RA} = $row->{Chr_Allele};
    $r{$F_TABIX_VA} = $row->{Alternative_Allele};
    $tp->add_row("-row" => \%r);
  }
}
$tp->finish();
