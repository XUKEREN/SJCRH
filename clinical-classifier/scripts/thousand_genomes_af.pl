#!/bin/env perl
# extract alleles and frequencies from 1,000 Genomes VCF
# - more managable output size
# - split and adjust compound alleles
# MNE 4/2016

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list get_hash_option log_message);
use FileUtils qw(universal_open);
use DelimitedFile;
use Reporter;
use TabixPrep;
use Cluster;
use FileHandle;
use TARTANUtils qw(tartan_genome_put_helper);

my $RAM = 2000;

my %FLAGS;

my @F_OUT = qw(
	     Chr
	     Pos
	     ReferenceAllele
	     MutantAllele
	     AF
	    );


my @clopts = (
	      "-dir=s",
	      # /datasets/public/genomes/hsapiens/hg19/SNPS/thousandgenomes/release_20130502/vcf_with_sample_level_annotation/
	      "-vcf=s",
	      # single file

	      "-step1",
	      # process individual chrom files

	      "-step2",
	      # merge to final file

	      "-max=i",
	      # debug
	      "-ping=i",
	      "-now",

	      "-genome=s",
	      "-tartan-put=s",

	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

if ($FLAGS{step1}) {
  run_step_1();
} elsif ($FLAGS{step2}) {
  run_step_2();
} elsif ($FLAGS{vcf}) {
  convert_vcf();
} elsif ($FLAGS{"tartan-put"}) {
  tartan_put();
}

sub run_step_1 {
  my $max = $FLAGS{"max"};
  my $vcf_dir = $FLAGS{dir} || die "-dir";
  my @gz = glob($vcf_dir . "/*.gz");

  foreach my $vcf (@gz) {
    my $cmd = sprintf '%s -vcf %s', $0, $vcf;
    $cmd .= sprintf ' -max %d', $FLAGS{max} if $FLAGS{max};
    $cmd .= sprintf ' -ping %d', $FLAGS{ping} if $FLAGS{ping};

    my $outfile = get_outfile($vcf);

    my $c = new Cluster(
			"-outfile" => $outfile,
			"-project" => "PCGP",
			#			  "-tracking_dir" => $CURRENT_DIR,
		       );

    $c->node_class("");
    #    $c->node_class("blade");
    # disable default idataplex requirement (good idea?)
    # large_mem OK, not sure about verna though

    #	$c->app("jdk6u30");
    # we do need java, but pre-load via module
#    $c->app("picard-lsf");
    # see if this helps w/killiness

    $c->memory_reserve_mb($RAM);
    $c->memory_limit_mb($RAM);
    $c->command($cmd);

    if ($FLAGS{now}) {
      if (-s $outfile) {
	printf STDERR "%s exists\n", $outfile;
      } else {
	printf STDERR "running %s\n", $cmd;
	system($cmd);
      }
    } else {
      $c->run();
    }
  }

}

sub new_tp {
  my ($outfile) = @_;
  die unless $outfile;
  my $tp = new TabixPrep(
			 "-outfile" => $outfile,
			 "-headers" => \@F_OUT,
			 "-header_chr" => "Chr",
			 "-header_start" => "Pos",
			);
#  $tp->input_sorted(1);
  # no longer true since we are adjusting alleles

  return $tp;
}

sub get_outfile {
  my ($vcf) = @_;
  return sprintf '%s.af.tab.gz', basename($vcf);
}

sub convert_vcf {
  my $vcf = $FLAGS{vcf} || die;
  my $outfile = get_outfile($vcf);

  my $cmd = sprintf 'vcf2tab.pl -file %s -1000 -stdout', $vcf;
  $cmd .= sprintf ' -max %d', $FLAGS{max} if $FLAGS{max};
  $cmd .= '|';

  my $tp = new_tp($outfile);
  my $fh = new FileHandle();
  $fh->open($cmd);
  my $df = new DelimitedFile("-fh" => $fh, "-headers" => 1);
  my $ping = $FLAGS{ping};
  my $count = 0;
  while (my $row = $df->get_hash()) {
    my %r;
    $r{Chr} = $row->{Chr} || die;
    $r{Pos} = $row->{WU_HG19_Pos} || die;
    $r{ReferenceAllele} = $row->{ReferenceAllele} || die;
    $r{MutantAllele} = $row->{MutantAllele} || die;
    $r{AF} = $row->{AF};
    dump_die($row, "no AF") unless defined $r{AF};
    # can be 0, e.g. 8       156419  rs370310067     G       A
    $tp->add_row("-row" => \%r);

    if ($ping and ++$count % $ping == 0) {
      log_message(sprintf "processed %d, at %s.%s", $count, @r{qw(Chr Pos)});
    }

  }
  $tp->finish();
}

sub run_step_2 {
  #
  # merge results from individual chromosomes into final file
  #
  my $vcf_dir = $FLAGS{dir} || die "-dir";
  my @gz = glob($vcf_dir . "/*.gz");

  my $outfile = "thousand_genomes_AF.tab.gz";

  my @processed;
  foreach my $vcf (@gz) {
    my $outfile = get_outfile($vcf);
    die "where is $outfile" unless -s $outfile;
    push @processed, $outfile;
  }

  my $tp = new_tp($outfile);
  $tp->input_sorted(1);
  my $ping = $FLAGS{ping};
  my $count = 0;
  foreach my $f (@processed) {
    printf STDERR "%s...\n", $f;
    my $df = new DelimitedFile(
			       "-file" => $f,
			       "-headers" => 1,
			     );
    $df->headers_raw->[0] =~ s/^#//;
    # undo tabix comment
    while (my $row = $df->get_hash()) {
      $tp->add_row("-row" => $row);
      if ($ping and ++$count % $ping == 0) {
	log_message(sprintf "processed %d, at %s.%s", $count, @{$row}{qw(Chr Pos)});
      }
    }
  }
  $tp->finish();
}

sub tartan_put {
  my $out_dir = $FLAGS{"tartan-put"};
  tartan_genome_put_helper(
			   "-genome" => $FLAGS{genome},
			   "-subdir" => "thousand_genomes_AF",
			   "-out" => $out_dir,
			  );
}
