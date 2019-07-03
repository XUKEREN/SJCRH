package ReadReport;
# caching wrapper to Bambino's ReadReport.java
# currently used to report on a single site: this is inefficient
# but simplifies embedding for ad-hoc programs
#
# see BAMReadReport for a version that works with a set of sites
#
# TO DO: add reference sanity check
# MNE 9/2014

use strict;
use Carp qw(confess);

use Digest::MD5 qw(md5_hex);
use File::Basename;

use ReferenceSanityCheck;
use MiscUtils qw(dump_die);
use Reporter;
use DelimitedFile;
use JavaRun;
use UniqueCacheFile;

use Configurable;
use Exporter;

@ReadReport::ISA = qw(Configurable Exporter);
@ReadReport::EXPORT_OK = qw();

use MethodMaker qw(
	cache
		  );

my $JAVA_RAM = 500;


sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  return $self;
}

sub get_info {
  my ($self, %options) = @_;
  my ($ref_name, $base_number, $bam);
  if (my $r = $options{"-row"}) {
    $ref_name = $r->{Chr} || die "chr";
    $bam = $r->{bam} || die "bam";
    $base_number = $r->{WU_HG19_Pos} || die "pos";
  }

  die "no ref name" unless $ref_name;
  die "no base num" unless $base_number;
  die "no bam" unless $bam;
  confess "where is $bam" unless -s $bam;

  my $ucf = new UniqueCacheFile("-prefix" => "RR");
  $ucf->add($bam, $ref_name, $base_number);
  my $outfile = $ucf->get_file();

  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
bam
reference
base_num
					)
				      ],
			 "-auto_qc" => 1,
			);
  my %r;
  $r{bam} = $bam;
  $r{reference} = $ref_name;
  $r{base_num} = $base_number;
#  dump_die(\%r, "debug", 1);
  $rpt->end_row(\%r);
  $rpt->finish();

  #
  #  call read report code:
  #
  my $outfile2 = $outfile . ".read_report.tab";
  my $needed = $self->cache ? not(-s $outfile2) : 1;
  if ($needed) {
    my $jr = new JavaRun("-classname" => "Ace2.ReadReport");
    $jr->ram($JAVA_RAM);
    my $cmd = sprintf '-min-mapq 0 -bam %s -passthrough %s', $bam, $outfile;
    my $cl = $jr->run("-command" => $cmd, "-return" => 1);

    system($cl);
    die "error running $cl: $?" if $?;
  }
  die unless -s $outfile2;

  print STDERR "reading $outfile2\n";
  my $df = new DelimitedFile(
			     "-file" => $outfile2,
			     "-headers" => 1,
			    );

  my @rows;
  while (my $r = $df->get_hash()) {
    push @rows, $r;
  }
  die unless @rows == 1;

  return $rows[0];
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
