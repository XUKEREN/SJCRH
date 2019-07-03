#!/bin/env perl
# convert flatfiles to tabix
# MNE

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use File::Basename;
use Cluster;
use CommandLineRebuilder;

use FileUtils qw(find_binary newer_than universal_open);
use MiscUtils qw(dump_die build_argv_list);
use TdtConfig;
use TabixPrep;
use TabixFile;
use DelimitedFile;
use Counter;

#my $F_POS = "pos(1-coor)";
# for 2.x, the only coordinates (i.e. hg19)
# for 3.0a, hg38 coordinates and the new default sort order.
# hg19 coordinates are "hg19_pos(1-based)"

my %FLAGS;
my @clopts = (
	      "-file=s",
	      "-stdin",
	      "-glob=s",

	      "-max=i",
	      "-f-chr=s",
	      "-f-pos=s",
	      "-force",
	      "-sort-needed=i",
	      "-out=s",

	      "-vcf",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

$FLAGS{"file"} = "-" if $FLAGS{stdin};

#my $infile = $FLAGS{file} || die "-file";
my @infiles;
if (my $f = $FLAGS{file}) {
  @infiles = $f;
} elsif (my $glob = $FLAGS{glob}) {
  @infiles = glob($glob);
  die "no files match $glob" unless @infiles;
}

if ($FLAGS{vcf}) {
  $FLAGS{"f-chr"} = "#CHROM";
  $FLAGS{"f-pos"} = "POS";
}

my $F_CHR = $FLAGS{"f-chr"} || die "specify -f-chr";
my $F_POS = $FLAGS{"f-pos"} || die "specify -f-pos";

my $sort_needed = $FLAGS{"sort-needed"};
die "-sort-needed [0|1]" unless defined $sort_needed;

my $outfile;

my $STDIN_MODE = $FLAGS{stdin} || $infiles[0] eq "-";

if ($STDIN_MODE) {
  $outfile = $FLAGS{out} || die "STDIN mode requires -out";
  die "outfile must end in .gz" unless $outfile =~ /\.gz$/;
} elsif (@infiles == 1) {
  my $bn = basename($infiles[0]);
  my $extra = $F_POS;
  $extra =~ s/\W/_/g;
  $extra =~ s/_$//;
  $outfile = sprintf '%s.%s.tabix.gz', $bn, $extra;
} else {
  $outfile = $FLAGS{out} || die "multiple-file input mode requires -out";
}
die "outfile must end in .gz" unless $outfile =~ /\.gz$/;

#
# get column names:
#
my @dfo = (
	   "-file" => $infiles[0],
	   "-headers" => 1
	  );
if ($FLAGS{vcf}) {
  push @dfo, "-skip_double_comments" => 1;
}

my $df = new DelimitedFile(@dfo);
my $headers = $df->headers_raw();
#die "1st header entry not commented" unless $headers->[0] =~ /^#/;

my $max_lines = $FLAGS{max};

my $needed;
if ($STDIN_MODE) {
  $needed = 1;
} else {
  $needed = -s $outfile ? 0 : 1;
  $needed = 1 if $FLAGS{force};
  foreach my $infile (@infiles) {
    $needed = 1 if newer_than($infile, $outfile);
  }
}

if ($needed) {
  if ($sort_needed) {
    #
    #  coordinates sorting required (i.e. for secondary position
    #  in alternate genome)
    #
    my $tp = new TabixPrep(
			   "-outfile" => $outfile,
			   "-header_chr" => $F_CHR,
			   "-header_start" => $F_POS,
			  );
#    $tp->init_headers_from_file($infiles[0]);
    $tp->headers($headers);
    if ($FLAGS{vcf}) {
      $tp->vcf_headers($df->double_comments());
    }

    my $c = new Counter(\@infiles);
    foreach my $infile (@infiles) {
      $c->next($infile);
      my $df;
      if ($STDIN_MODE) {
	$df = new DelimitedFile(
				"-file" => $infile,
			       );
	$df->headers_raw($headers);
      } else {
	$df = new DelimitedFile(
				"-file" => $infile,
				"-headers" => 1,
			       );
      }
      my $lines = 0;
      while (my $row = $df->get_hash()) {
        $tp->add_row("-row" => $row);
	last if $max_lines and ++$lines >= $max_lines;
      }
    }
    $tp->finish();
  } else {
    #
    #  primary coordinates are already sorted; much faster to build.
    #
    my $first = 1;
    # my $c = new Counter(\@infiles);
    my $tmp = $outfile . ".tmp";
    open(OUT, sprintf "|bgzip > %s", $tmp) || die;
    my $c = new Counter(\@infiles);
    foreach my $infile (@infiles) {
      $c->next($infile);
      my $fh = universal_open($infile);
      if ($first) {
	$first = 0;
      } else {
	# skip 2nd instance of header line in subsequent files
	my $header = <$fh>;
      }
      my $lines = 0;
      while (<$fh>) {
	chomp;
	s/\r$//;
	# remove DOS-style carriage returns
	print OUT $_ . "\n";
	last if $max_lines and ++$lines >= $max_lines;
      }
      close $fh;
    }
    close OUT;
    die if $?;
    rename($tmp, $outfile) || die;
  }
}

unless ($FLAGS{vcf}) {
  # VCF already done
  my $tf = new TabixFile(
			 "-file" => $outfile,
			 "-index" => 1,
			 "-f_chr" => $F_CHR,
			 "-f_start" => $F_POS,
			);
}

