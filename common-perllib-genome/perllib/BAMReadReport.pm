package BAMReadReport;
# wrapper to Bambino Ace2.ReadReport utility
# MNE 2/2015

use strict;
use Configurable;
use Exporter;

use File::Basename;
use MiscUtils qw(dump_die);
use List::Util qw(sum);

@BAMReadReport::ISA = qw(Configurable Exporter);
@BAMReadReport::EXPORT_OK = qw();

use MethodMaker qw(
bam
sites
cache
min_quality
min_mapping_quality
results
cache_basename
extra
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->sites({});
  $self->configure(%options);
  return $self;
}

sub add_site {
  # add a location to check
  my ($self, %options) = @_;
  my $ref = $options{"-ref"} || die;
  my $pos = $options{"-base"} || die;
  $self->sites()->{$ref}{$pos} = 1;
}

sub get_report {
  #
  # generate and parse report
  #
  my ($self, %options) = @_;

  my $cache = $self->cache();
  my $bam = $self->bam() || die "-bam";

  my $base = $self->cache_basename || "";
  $base .= "." if $base;

  my $infile = $base . basename($bam) . ".rr_config.tab";
  my $outfile = $infile . ".read_report.tab";
  my $in_needed = $cache ? not(-s $infile) : 1;
  my $out_needed = $cache ? not(-s $outfile) : 1;

  if ($in_needed) {
    my $rpt = new Reporter(
			   "-file" => $infile,
			   "-delimiter" => "\t",
			   "-labels" => [
					 qw(
					     reference
					     base_num
					  )
					],
			   "-auto_qc" => 1,
			  );
    my $sites = $self->sites() || die;
    foreach my $ref (sort keys %{$sites}) {
      # generate infile by ref name and position for efficiency
      foreach my $pos (sort {$a <=> $b} keys %{$sites->{$ref}}) {
	my %r;
	$r{reference} = $ref;
	$r{base_num} = $pos;
	$rpt->end_row(\%r);
      }
    }
    $rpt->finish();
  }

  if ($out_needed) {
    my $cmd = sprintf 'java -Xmx500m Ace2.ReadReport -bam %s -passthrough %s',
      $bam, $infile;
    if (my $mq = $self->min_quality) {
      $cmd .= " -min-quality " . $mq;
    }
    my $min_mapq = $self->min_mapping_quality;
    if (defined $min_mapq) {
      $cmd .= " -min-mapq " . $min_mapq;
    }

    $cmd .= " " . $self->extra if $self->extra;

#    die $cmd;
    printf STDERR "running: %s\n", $cmd;
    system $cmd;
    die "ERROR running $cmd" if $?;
  }
  die "where is $outfile" unless -s $outfile;

  my $df = new DelimitedFile("-file" => $outfile,
			     "-headers" => 1,
			     );
  my %results;
  while (my $row = $df->get_hash()) {
    my $ref = $row->{reference} || die;
    my $pos = $row->{base_num} || die;

    my $coverage = sum(map {$row->{$_}} qw(A C G T));
    $row->{coverage} = $coverage;

    $results{$ref}{$pos} = $row;
  }
  $self->results(\%results);
  return \%results;
}

sub get_coverage {
  my ($self, %options) = @_;
  my $ref = $options{"-ref"} || die;
  my $pos = $options{"-base"} || die;
  $self->get_report() unless $self->results();
  my $results = $self->results() || die;
  return $results->{$ref}{$pos} || die "WTF: no coverage entry for $ref $pos";
}

sub is_suspicious_xt_m_ratio {
  # DEVELOPMENT
  my ($self, %options) = @_;
  my $row = $self->get_coverage(%options);
  my $var_base = $options{"-nt"} || die;
  my $SUSPICIOUS_RATIO = 0.30;
  # maybe low especially if just a few supporting reads?
  # OTOH raising it to 0.75 misses some known paralogous sites, e.g.
  # chr2    73677421        A       G       3       suspicious_variant_tumor_XT_M_ratio=0.33        http://bamviewer-rt:8080/BAMViewer/aceview/splash?tumorname=/rgs01/resgen/prod/tartan/index/data/SCMC/SCMC/SJALL018375_O2/VALCAP/bam/SJALL018375_O2.bam&ref=hg19&region=chr2&center=73677421&fullPath=true
  # this has many M reads that are tossed for other reasons (e.g. clips),
  # only a few reads left, 1 M.
  # => should we count XT status before read rejection?

  my $suspicious = 0;
  if (my $all = $row->{"XT_" . $var_base}) {
    my %counts;
    foreach my $entry (split /,/, $all) {
      my @f = split /=/, $entry;
      die unless @f == 2;
      $counts{$f[0]} = $f[1];
    }
    my $total = sum values %counts;
    my $xt_m = $counts{M} || 0;
    my $freq = $xt_m / $total;
    $suspicious = $freq if $freq > $SUSPICIOUS_RATIO;
  }
  return $suspicious;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
