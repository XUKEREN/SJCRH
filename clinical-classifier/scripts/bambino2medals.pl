#!/bin/env perl
# run Bambino calls through Annovar+ and medal ceremony
# MNE 8/2016

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use FileUtils qw(find_binary);
use MiscUtils qw(dump_die build_argv_list get_hash_option);
use DelimitedFile;
use Reporter;
use OutfileStack;
use CommandLineRebuilder;
use Cluster;

my $RAM = 10000;
# TO DO: tune annovar+ cache size down for smaller memory request?
# latest medal ceremony code no longer requires 10g

my $ANNOVAR_MAX_TRANSCRIPT_CACHE = 100;
# tweak to lower memory usage

my %FLAGS;
my @clopts = (
	      "-file=s",
	      "-files=s",
	      "-glob=s",
	      "-single=s",

	      "-genome=s",

	      "-ram=i" => \$RAM,
	      "-annovar-no-scratch",
	      "-annovar-max-transcript-cache=i",
	      "-now",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $genome = $FLAGS{genome} || die "-genome";
my @binaries = qw(
		   bambino2annovar.pl
		   annotate_variation.pl
		   tabix
		   medal_ceremony.pl
		);

foreach my $bin (@binaries) {
  find_binary($bin, "-die" => 1);
}


if (my $one = $FLAGS{single}) {
  run_one($one);
} else {
  my $files = build_argv_list("-flags" => \%FLAGS,
			      "-single" => "file",
			      "-set" => "files",
			      "-glob" => "glob",
			     );

  foreach my $infile (@{$files}) {
    my $clrb = new CommandLineRebuilder(
					"-parameters" => \@clopts,
					"-flags" => \%FLAGS,
				       );
    $clrb->exclude_parameter("-file");
    $clrb->exclude_parameter("-files");
    $clrb->exclude_parameter("-glob");
    my $cmd = $clrb->get_command_line("-single" => $infile);

    my $outfile = sprintf "%s.annovar_merged.tab.mp.tab.medals.tab.medals.tab", basename($infile);

    if ($FLAGS{now}) {
      printf STDERR "running: %s\n", $cmd;
      system $cmd;
    } else {
      my $c = new Cluster(
			  "-outfile" => $outfile,
			  "-project" => "PCGP",
			 );

      $c->memory_reserve_mb($RAM);
      $c->memory_limit_mb($RAM);
      $c->command($cmd);
      $c->run();
    }
  }
}

sub run_one {
  #
  # process on a single node
  #
  my ($infile) = @_;

  my $os = new OutfileStack(
			    "-use_basename" => 1,
			    "-start_file" => $infile,
			   );

  #
  # Annovar+:
  #
  my $template = "bambino2annovar.pl -genome $genome -file %s";
  $template .= sprintf " -max-transcript-cache %d", $ANNOVAR_MAX_TRANSCRIPT_CACHE;
  $template .= " -scratch" unless $FLAGS{"annovar-no-scratch"};
  $os->add_level_and_run(
			 "-suffix" => ".annovar_merged.tab",
			 "-template" => $template,
			);

  #
  # medal prep:
  #
  $os->add_level_and_run(
			 "-suffix" => ".mp.tab",
			 "-template" => "annovar2medals.pl -genome $genome -file %s"
			);


  #
  # somatic medals:
  #
  $os->add_level_and_run(
			 "-suffix" => ".medals.tab",
			 "-template" => "medal_ceremony.pl -genome $genome -no-pvr -single-snv %s"
			);

  #
  # germline medals:
  #
  $os->add_level_and_run(
			 "-suffix" => ".medals.tab",
			 "-template" => "medal_ceremony.pl -genome $genome -no-pvr -single-gl %s"
			);

  exit(0);
  # success

}

