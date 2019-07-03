#!/bin/env perl
# describe me

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list get_hash_option);
use DelimitedFile;
use Reporter;
use TdtConfig;
use FileUtils qw(read_simple_file);

my %FLAGS;
my @clopts = (
	      "-reportable",
	      "-genome=s",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $gene_list;
if ($FLAGS{reportable}) {
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  $gene_list = read_simple_file($config_genome->{CLINCLS_GERMLINE_REPORTABLE_GENES}) || die;
}
die "-reportable" unless $gene_list;

die scalar @{$gene_list};

# - map query genes to symbols used in refflat
# - foreach refflat transcript, bucket into (a) query list and (b) other list
# - foreach refflat transcript in query transcript set, look for overlap
#   with transcripts in other set
# - various ways to count: txStart/end, CDS start/end, exons only


