#!/bin/env perl
# wrapper for basic somatic and germline classification

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

#use lib "/nfs_exports/apps/internal/script_bin/packages/edmonson_perl_lib/";
#use DevelopmentPath;
# optionally replaces production path with development version

use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die log_message);
use FileUtils qw(newer_than read_simple_file);
use CommandLineRebuilder;
use Cluster;
use TemporaryFileWrangler;

my $CLUSTER_RAM = 10000;

my %FLAGS;

my $CONFIG = "/nfs_exports/genomes/1/projects/ClinicalSeq/medal_config/2013_09_20/medal_config_2013_09_20.txt";

my $GENOME = "GRCh37-lite";
# only one supported for now

my @clopts = (
	      "-file=s",
	      "-no-sort",
	      "-list=s",
	      "-force",

	      "-cluster=s",
	      # OBSOLETE

	      "-debug",

	      "-outfile-fq",
	      "-no-germline",
	      "-bsub",

	      "-gl-custom-genes=s",
	      "-gl-no-reviewable",
	      "-gl-custom-allow-interval-fail",
	      # germline passthrough options (hack)
	     );

GetOptions(
	   \%FLAGS,
	   @clopts
	  );

if ($FLAGS{bsub}) {
  my $clrb = new CommandLineRebuilder(
				      "-parameters" => \@clopts,
				      "-flags" => \%FLAGS
				     );
  $clrb->exclude_parameter("-bsub");
  my $cmd = $clrb->get_command_line();

  my $c = new Cluster(
		      "-outfile" => "bogus_no_outfile.$$",
		      "-project" => "PCGP",
		     );

  $c->node_class("");
  $c->memory_reserve_mb($CLUSTER_RAM);
  $c->memory_limit_mb($CLUSTER_RAM);
  $c->command($cmd);
  $c->run();
  exit(0);
}

my $SOMATIC_LIST_FILE;
my $NEED_SOMATIC;
my $NEED_GERMLINE;

my $tfw = new TemporaryFileWrangler();

if (my $lf = $FLAGS{list}) {
  $SOMATIC_LIST_FILE = $lf;
} elsif (my $f = $FLAGS{file}) {
  $SOMATIC_LIST_FILE = $tfw->get_tempfile("-append" => ".files_to_medal.txt");
  open(TMP, ">" . $SOMATIC_LIST_FILE) || die;
  printf TMP "%s\n", $f;
  close TMP;
} else {
  die "specify -file [file] or -list [listfile]\n";
}

my $in_files = read_simple_file($SOMATIC_LIST_FILE);

my $need_somatic;
$need_somatic = 1 if $FLAGS{force};
my @somatic_out;
foreach my $infile (@{$in_files}) {
  my $outfile = get_outfile($infile);
  push @somatic_out, $outfile;
  if (-s $outfile) {
    $need_somatic = 1 if newer_than($infile, $outfile);
  } else {
    $need_somatic = 1;
  }
}

my @mc_args;

if ($FLAGS{"no-sort"}) {
  push @mc_args, "-no-sort";
} else {
  printf STDERR "NOTE: using default behavior (sorted output)\n";
}
push @mc_args, "-outfile-fq" if $FLAGS{"outfile-fq"};
my $mc_string = join " ", @mc_args;


if ($need_somatic) {
  my $log_base = sprintf 'somatic.%d', $$;
  my $cmd = sprintf '%s -no-pvr -in-snv-indel %s %s > %s.out 2>%s.err',
    get_cmd_base(),
      $SOMATIC_LIST_FILE,
	$mc_string,
	  $log_base, $log_base;
  run_command($cmd);
}

my $need_germline;
$need_germline = 1 if $FLAGS{force};
$need_germline = 0 if $FLAGS{"no-germline"};
my $GERMLINE_LIST_FILE = $tfw->get_tempfile("-append" => ".germline.lst");
open(TMP, ">" . $GERMLINE_LIST_FILE) || die;
foreach my $infile (@somatic_out) {
  printf TMP "%s\n", $infile;
#  my $outfile = $infile . ".medals.tab";
  my $outfile = get_outfile($infile);
  if (-s $outfile) {
    $need_germline = 1 if newer_than($infile, $outfile);
  } else {
    $need_germline = 1;
  }
}
close TMP;

if ($need_germline) {
  my $log_base = sprintf 'germline.%d', $$;
  my $cmd = sprintf '%s -no-pvr -in-gl %s %s ',
    get_cmd_base(),
      $GERMLINE_LIST_FILE,
	$mc_string;

  my @extra;
  push @extra, "-gl-no-reviewable" if $FLAGS{"gl-no-reviewable"};
  push @extra, "-gl-custom-allow-interval-fail" if $FLAGS{"gl-custom-allow-interval-fail"};
  push @extra, "-gl-custom-genes", $FLAGS{"gl-custom-genes"} if $FLAGS{"gl-custom-genes"};
  $cmd .= sprintf " %s ", join " ", @extra if @extra;
  $cmd .= sprintf ' > %s.out 2>%s.err', $log_base, $log_base;

  run_command($cmd);
}

sub run_command {
  my ($cmd) = @_;
  if ($FLAGS{debug}) {
    printf STDERR "DEBUG, NOT running: %s\n", $cmd;
  } else {
    log_message("running $cmd");
    system $cmd;
    die "error executing $cmd" if $?;
  }
}

sub get_cmd_base {
  my $cmd;
  if (my $cluster = $FLAGS{cluster}) {
    # now obsolete: debugging of old behavior only
    if ($cluster eq "research-old") {
      $cmd = sprintf 'medal_ceremony.pl -config %s', $CONFIG;
    } elsif ($cluster eq "clinical-old") {
      $cmd = "medal_ceremony_using_configs.sh GRCh37-lite";
    } else {
      die "-cluster is obsolete, use -genome";
    }
  } else {
    $cmd = "medal_ceremony.pl -genome $GENOME";
  }
  return $cmd;
}

sub get_outfile {
  my ($infile) = @_;
  my $outfile = sprintf '%s.medals.tab', $FLAGS{"outfile-fq"} ? $infile : basename($infile);
  return $outfile;
}
