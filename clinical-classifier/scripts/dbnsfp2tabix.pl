#!/bin/env perl
# convert dbNSFP files to tabix format
#
# TO DO:
# - some indication of target build/field index in output filename?

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use File::Basename;
use FileHandle;

use Cluster;
use CommandLineRebuilder;
use FileUtils qw(find_binary newer_than);
use MiscUtils qw(dump_die build_argv_list);
use TdtConfig;
use TabixPrep qw(tabix_concatenate);
use TabixFile;
use DelimitedFile;
use Counter;
use TARTANUtils qw(tartan_genome_put_helper);
use CommandLineRebuilder;

my $RAM = 4000;

#my $F_POS = "pos(1-coor)";
# for 2.x, the only coordinates (i.e. hg19)
# for 3.0a, hg38 coordinates and the new default sort order.
# hg19 coordinates are "hg19_pos(1-based)"

my %FLAGS;
my @clopts = (

	      "-tartan-index=s",
	      # output dir

	      "-dir=s",
	      "-genome=s",
	      "-max=i",
	      "-f-chr=s",
	      "-f-pos=s",
	      "-force",

	      "-sort-needed=i",
	      "-zip=s",
	      "-ping=i",

	      "-chr=s",
	      "-split",
	      "-join",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

if (my $out_dir = $FLAGS{"tartan-index"}) {
  tartan_genome_put_helper(
			   "-genome" => $FLAGS{genome},
			   "-subdir" => "SUPPORT/dbNSFP_tabix",
			   "-out" => $out_dir,
			  );
} elsif ($FLAGS{split}) {
  split_prep();
  # step 1: convert each chr individually
} elsif ($FLAGS{join}) {
  join_split_files();
  # step 2: join split files to single output file
} else {
  main_loop();
}

sub main_loop {
  my $F_CHR = $FLAGS{"f-chr"} || die "specify -f-chr";
  my $F_POS = $FLAGS{"f-pos"} || die "specify -f-pos";
  my $sort_needed = $FLAGS{"sort-needed"};
  die "-sort-needed [0|1]" unless defined $sort_needed;
  my $zipfile = $FLAGS{zip};

  my $infiles = get_infiles();
  my $outfile = get_outfile($infiles->[0]);

  my $df = get_df($infiles->[0], $zipfile);
  my $headers = $df->headers_raw();
  die "1st header entry not commented" unless $headers->[0] =~ /^#/;

  my $max_lines = $FLAGS{max};

  my $needed = -s $outfile ? 0 : 1;
  $needed = 1 if $FLAGS{force};
  foreach my $infile (@{$infiles}) {
    $needed = 1 if newer_than($infile, $outfile);
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
			     "-headers" => $headers,
			    );

      my $c = new Counter($infiles);
      my $total_rows = 0;
      my $ping = $FLAGS{ping};
      foreach my $infile (@{$infiles}) {
	my $df = get_df($infile, $zipfile);
	my $lines = 0;
	my $skipped = 0;
	while (my $row = $df->get_hash()) {
	  if ($row->{$F_POS} eq ".") {
	    # skip rows with no mapping in target genome
	    $skipped++;
	  } elsif ($row->{$F_CHR} eq ".") {
	    dump_die($row, "fix me: invalid $F_CHR");
	  } else {
	    $tp->add_row("-row" => $row);
	    last if $max_lines and ++$lines >= $max_lines;
	  }
	  if ($ping and ++$total_rows % $ping == 0) {
	    printf STDERR "%s: processed %d, %s\n", scalar(localtime), $total_rows, $df->last_line;
	  }
	}
	$c->next($infile);
      }
      $tp->finish();

    } else {
      #
      #  use primary coordinates, which are already sorted.
      #  Much faster to build.
      #
      die "zipfile mode not implemented yet" if $zipfile;

      my $first = 1;
      # my $c = new Counter(\@infiles);
      my $tmp = $outfile . ".tmp";
      open(OUT, sprintf "|bgzip > %s", $tmp) || die;
      my $c = new Counter($infiles);
      foreach my $infile (@{$infiles}) {
	open(IN, $infile) || die;
	if ($first) {
	  $first = 0;
	} else {
	  # skip 2nd instance of header line in subsequent files
	  my $header = <IN>;
	}
	my $lines = 0;
	while (<IN>) {
	  chomp;
	  s/\r$//;
	  # remove DOS-style carriage returns
	  print OUT $_ . "\n";
	  last if $max_lines and ++$lines >= $max_lines;
	}
	close IN;
	$c->next($infile);
      }
      close OUT;
      die if $?;
      rename($tmp, $outfile) || die;
    }
  }

  my $tf = new TabixFile(
			 "-file" => $outfile,
			 "-index" => 1,
			 "-f_chr" => $F_CHR,
			 "-f_start" => $F_POS,
			);
}


sub get_zip_fh {
  my ($zipfile, $fn) = @_;
  my $cmd = sprintf 'unzip -p %s %s|', $zipfile, $fn;
  my $fh = new FileHandle();
  $fh->open($cmd) || die;
  return $fh;
}

sub get_df {
  my ($file, $zipfile) = @_;
  my @dfo = ("-headers" => 1);
  if ($zipfile) {
    push @dfo, "-fh" => get_zip_fh($zipfile, $file);
  } else {
    push @dfo, "-file" => $file;
  }
  return new DelimitedFile(@dfo);
}

sub get_outfile {
  my ($file) = @_;
  my $F_POS = $FLAGS{"f-pos"} || die "specify -f-pos";
  my $chr = $FLAGS{chr};

  basename($file) =~ /^(.*)_/ || die;
  my $bn = $1;
  my $extra = $F_POS;
  $extra =~ s/\W/_/g;
  $extra =~ s/_$//;
  $extra .= "_" . $chr if $chr;

  my $outfile = sprintf '%s.%s.tabix.gz', $bn, $extra;
  return $outfile;
}

sub split_prep {
  my $clrb = new CommandLineRebuilder(
                                      "-parameters" => \@clopts,
                                      "-flags" => \%FLAGS,
                                     );
  $clrb->exclude_parameter("-split");
  my $infiles = get_infiles();

  foreach my $chr (1..22, qw(X Y M)) {
    $FLAGS{chr} = $chr;
    my $cmd = $clrb->get_command_line("-chr" => $chr);
    my $outfile = get_outfile($infiles->[0]);
#    printf "%s => %s\n", $cmd, $outfile;

    my $c = new Cluster(
			"-outfile" => $outfile,
			"-project" => "G4K",
			#                     "-tracking_dir" => $CURRENT_DIR,
		       );

    $c->node_class("");
    $c->memory_reserve_mb($RAM);
    $c->memory_limit_mb($RAM);
    $c->command($cmd);
    $c->run();
  }
}

sub get_infiles {
  my $zipfile = $FLAGS{zip};

  my @infiles;
  if ($zipfile) {
    open(IN, "unzip -l $zipfile|") || die;
    # hack: Archive::Zip not on clinical
    my $in_list;
    while (<IN>) {
      chomp;
      if (/^\---/) {
	if ($in_list) {
	  last;
	} else {
	  $in_list = 1;
	}
      } elsif ($in_list) {
	s/^\s+//;
	my @f = split /\s+/, $_;
	my $file = $f[3];
#	printf STDERR "in list: %s\n", $file;
	push @infiles, $file if $file and $file =~ /variant\./;
      }
    }
    @infiles = sort @infiles;
  } else {
    my $nsfp_dir = $FLAGS{dir} || die "specify -dir DIRECTORY (unpacked tabix files)\n";
    @infiles = sort glob($nsfp_dir . "/*variant.*");
  }
  die unless @infiles;
  my $chr = $FLAGS{chr};
  if ($chr) {
    my @hits = grep {/\.chr${chr}$/} @infiles;
    die "ambiguous hits for $chr in file list" if @hits > 1;
    die "can't find $chr in file list: " . join ",", @infiles unless @hits == 1;
    @infiles = @hits;
  }
  die "ERROR: no infiles" unless @infiles;

  return \@infiles;
}

sub join_split_files {
  my $F_POS = $FLAGS{"f-pos"} || die "specify -f-pos";

  my @outfiles;
  my $error;
  my $infiles = get_infiles();
  foreach my $chr (1..22, qw(X Y M)) {
    $FLAGS{chr} = $chr;
    my $outfile = get_outfile($infiles->[0]);
    if (-s $outfile) {
      push @outfiles, $outfile;
    } else {
      printf STDERR "ERROR: missing %s\n";
      $error = 1;
    }
  }

  unless ($error) {
    delete $FLAGS{chr};
    my $outfile = get_outfile($infiles->[0]);

    tabix_concatenate(
		      "-in" => \@outfiles,
		      "-out" => $outfile,
		      "-sorted" => 1,
		      "-ping" => 100000,
		      "-header_chr" => "#chr",
		      "-header_start" => $F_POS,
#		      "-max" => 100
		     );
  }

}
