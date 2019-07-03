#!/bin/env perl
# retrieve medal run files from TARTAN index for regression testing

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use File::Find;
use File::Basename;
use File::Copy;
use File::Path;
use File::Copy;
use Getopt::Long;

use MiscUtils qw(dump_die);
use FileUtils qw(read_simple_file write_simple_file);
use DelimitedFile;
use Reporter;
use Counter;

my %FLAGS;
GetOptions(
	   \%FLAGS,
	   "-type=s",
	   # cnv/sv/cnv

	   "-germline",
	   "-somatic",

	   "-wc",
	   "-rc",

	   "-build-single-variant-list=s",


	   "-diff-annotations",
	   "-dev=s",
	   "-prod=s",

	  );

if (my $fn = $FLAGS{"build-single-variant-list"}) {
  build_single_variant_list($fn);
  exit(0);
} elsif ($FLAGS{"diff-annotations"}) {
  diff_annotations();
  exit(0);
}


# TO DO:
# replace with faster tartan_finder.pl calls, e.g.
# /home/medmonso/work/medals/databases/gedi/legacy/find_raw_files.sh

my @roots = (
	     '/cgs01/clingen/prod/tartan/index/data/ClinicalPilot/ClinicalPilot/',
	     '/cgs01/clingen/prod/tartan/index/data/Clinical/2015/',
	     '/cgs01/clingen/prod/tartan/index/data/Clinical/2016/',
	    );
my $type = $FLAGS{type} || die "-type";

my $subdir;
my $output_column_remove;
if ($type eq "cnv") {
  my $sub = $FLAGS{germline} ? "germline" : $FLAGS{somatic} ? "somatic" : die "-germline/-somatic";
  $subdir = sprintf 'cnv-%s-classified', $sub;
  $output_column_remove = "Number ofGenes";
  # sic
} elsif ($type eq "sv") {
  $subdir = sprintf 'sv-somatic-classified';
  # somatic only
  $output_column_remove = "GSBClass";
} elsif ($type eq "snv") {
  my $sub = $FLAGS{germline} ? "germline" : $FLAGS{somatic} ? "somatic" : die "-germline/-somatic";
  $subdir = sprintf 'snv-%s-classified', $sub;
  if ($FLAGS{somatic}) {
    $output_column_remove = "dbSNP";
  } elsif ($FLAGS{germline}) {
    $output_column_remove = "paneldecision";
  } else {
    die;
  }
} elsif ($type eq "indel") {
  # no subdir
  my $sub = $FLAGS{germline} ? "germline" : $FLAGS{somatic} ? "somatic" : die "-germline/-somatic";
  $subdir = sprintf 'indel-%s-classified', $sub;
  if ($FLAGS{germline}) {
    $output_column_remove = "paneldecision";
  } else {
    die "only germline supported";
  }
} elsif ($FLAGS{germline}) {
  die "type must be snv/cnv/sv";
}

die "no subdir" unless $subdir;
my $cache_file = $subdir . ".cache";

my @fq;

if ($FLAGS{rc}) {
  # use cached results
  my $set = read_simple_file($cache_file);
  die unless @{$set};
  @fq = @{$set};
} else {
  my $matching_dirs = 0;

  my $callback = sub {
    if (-d $_ and $_ eq $subdir) {
      # faster to find directory of interest rather than following
      # symlinks and looking for medal files.
      $matching_dirs++;
      my $dir = join "/", $File::Find::dir, $_;
      my @files;
      if ($type eq "snv") {
	if ($FLAGS{somatic}) {
	  @files = glob($dir . "/*_smcls.txt");
	} elsif ($FLAGS{germline}) {
	  @files = glob($dir . "/*_glcls.txt");
	} else {
	  die;
	}
      } elsif ($type eq "indel") {
	@files = glob($dir . "/*_smcls_glcls.txt");
      } else {
	@files = glob($dir . "/*medals.tab");
      }
      if (@files) {
	push @fq, @files;
	printf STDERR "found in %s, total %d\n", $dir, scalar @fq;
      } else {
	printf STDERR "no medal files in %s\n", $dir;
      }

      printf "dirs:%d files:%d\n", $matching_dirs, scalar @fq;

    }
  };

  foreach my $root (@roots) {
    printf STDERR "scanning %s\n", $root;
    find($callback, $root);
  }


  if ($FLAGS{wc}) {
    write_simple_file(\@fq, $cache_file);
 }
}

my %saw;
foreach my $fq (@fq) {
#  printf STDERR "input: %s\n", $fq;
  my $bn = basename($fq);
  my $tag;
  if ($fq =~ /\/(WHOLE_GENOME|EXOME)\//)  {
    # sometimes same basename for multiple sequencing types
    $tag = $1;
  }
  $bn .= "." . $tag if $tag;
  $bn .= ".medals.tab" unless $bn =~ /medals\.tab$/;
  # final SNV files (Matt's scripts?) don't use the medals.tab suffix

  die "duplicate $bn" if $saw{$bn};
  $saw{$bn} = 1;

  # strip output columns to create input file:
  my $copy_infile = sprintf '%s/input/%s', $subdir, $bn;
  if ($copy_infile =~ s/\.medals.tab$//) {
    # OK
  } else {
    die "can't strip suffix for $copy_infile";
  }

  # copy medal'd results file:
  my $copy_outfile = sprintf '%s/output/%s', $subdir, $bn;
  create_dir_for_file($copy_outfile);
  copy($fq, $copy_outfile) || die "copy $fq to $copy_outfile failed: $!";

  create_dir_for_file($copy_infile);
  my $cmd = sprintf 'report_excerpt.pl -in %s -out %s -exclude-after "%s"',
    $copy_outfile, $copy_infile, $output_column_remove;
  system $cmd;
  die "ERROR $? running $cmd" if $?;

  if (($type eq "snv" or $type eq "indel") and
      $FLAGS{germline}) {
    # special handling required to undo column renaming
    my $tmpfile = $copy_infile . ".tmp";
    open(IN, $copy_infile) || die;
    my $first = <IN>;
    my $mod = $first;
    $mod =~ s/Somatic_//g;
    open(OUT, ">" . $tmpfile) || die;
    print OUT $mod;
    while (<IN>) {
      print OUT;
    }
    close IN;
    close OUT;

    copy($tmpfile, $copy_infile) || die;
    unlink $tmpfile;
  }

}

sub create_dir_for_file {
  my ($fn) = @_;
  my $dir = dirname($fn);
  unless (-d $dir) {
    mkpath($dir) || die;
  }
}

sub build_single_variant_list {
  my ($flist) = @_;
  my $files = read_simple_file($flist);

  my @wanted_fields = qw(
			  GeneName
			  Chr
			  WU_HG19_Pos
			  ReferenceAllele
			  MutantAllele
			  Class
			  AAChange
			  mRNA_acc
		       );

  my %saw;
  my $wrote = 0;
  my $skipped = 0;

  my $c = new Counter($files);

  my %track;
  # bucket by chrom and position for:
  # (a) more efficient medaling/file caching
  # (b) avoid concentrations of indels which are slower to process/mux

  foreach my $file (@{$files}) {
    my $df = new DelimitedFile("-file" => $file,
			       "-headers" => 1,
			      );
    while (my $row = $df->get_hash()) {
      foreach my $f (@wanted_fields) {
	if (not($row->{WU_HG19_Pos}) and $row->{WU_HG18_Pos}) {
	  $row->{WU_HG19_Pos} = $row->{WU_HG18_Pos};
	}

	dump_die($row, "missing $f") unless exists $row->{$f};

	my $chr = $row->{Chr} || die;
	my $pos = $row->{WU_HG19_Pos} || die;
	my $ra = $row->{ReferenceAllele} || die;
	my $va = $row->{MutantAllele} || die;

	my $snv4 = join ".", $chr, $pos, $ra, $va;
	if ($saw{$snv4}) {
	  $skipped++;
	} else {
	  $saw{$snv4} = 1;
	  push @{$track{$chr}{$pos}}, $row;
	}
      }
    }
    $c->next($file);
  }

  my $outfile = "unique.tab";
  my $rpt = new Reporter(
                         "-file" => $outfile,
                         "-delimiter" => "\t",
                         "-labels" => \@wanted_fields,
                         "-auto_qc" => 1,
                        );
  foreach my $chr (sort keys %track) {
    foreach my $pos (sort {$a <=> $b} keys %{$track{$chr}}) {
      foreach my $row (@{$track{$chr}{$pos}}) {
	$wrote++;
	$rpt->end_row($row);
      }
    }
  }

  printf STDERR "wrote:%d skipped:%d\n", $wrote, $skipped;
  $rpt->finish();
}

sub diff_annotations {
  my $f_dev = $FLAGS{dev} || die;
  my $f_prod = $FLAGS{prod} || die;
  my %base_fields = map {$_, 1} qw(
				    GeneName
				    Chr
				    WU_HG19_Pos
				    ReferenceAllele
				    MutantAllele
				    Class
				    AAChange
				    mRNA_acc
				 );

  my $rows_dev = load_rows($f_dev);
  my $rows_prod = load_rows($f_prod);

  my %f_dev = map {$_, 1} keys %{$rows_dev->[0]};
  my %f_prod = map {$_, 1} keys %{$rows_prod->[0]};
  my %f_all = (%f_dev, %f_prod);

  foreach my $f (sort keys %f_all) {
    next if $base_fields{$f};

    unless ($f_dev{$f} and $f_prod{$f}) {
      printf STDERR "skipping column %s: not in both files\n", $f;
      next;
    }

    my $out_dev = write_sorted($rows_dev, $f, "dev");
    my $out_prod = write_sorted($rows_prod, $f, "prod");

    printf "diff of field %s:\n", $f;
    system "diff $out_prod $out_dev";

  }
}

sub load_rows {
  my ($infile) = @_;
  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  my @rows;
  while (my $row = $df->get_hash()) {
    push @rows, $row;
  }
  return \@rows;
}

sub write_sorted {
  my ($rows, $field, $tag) = @_;
  my $outfile = sprintf 'diff_%s_%s', $field, $tag;
  open(OUT, "|sort > $outfile") || die;
  foreach my $row (@{$rows}) {
    printf OUT "%s\n", join "\t", map {$row->{$_}} qw(Chr WU_HG19_Pos ReferenceAllele MutantAllele), $field;
  }
  close OUT;
  return $outfile;
}
