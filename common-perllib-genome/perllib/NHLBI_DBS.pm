package NHLBI_DBS;
# binary search of NHLBI records

use strict;
use Configurable;
use Exporter;

use NHLBIParser;
use GenomeUtils qw(cook_chromosome_name);
use MiscUtils qw(dump_die);
use DelimitedBinarySearch;

@NHLBI_DBS::ISA = qw(Configurable Exporter);
@NHLBI_DBS::EXPORT_OK = qw();

use MethodMaker qw(
	directory
chr2file
np
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self, %options) = @_;
  my $dir = $self->directory() || die;
  my @files = glob(sprintf '%s/ESP6500.*.snps.txt', $dir);
  die "no NHLBI files in $dir" unless @files;
  my %chr2file;
  foreach my $fn (@files) {
    $fn =~ /(chr\w+)\./ || die;
    my $chr = cook_chromosome_name($1) || die;
    die "duplicate" if $chr2file{$chr};
    $chr2file{$chr} = $fn;
  }
  my $np = new NHLBIParser("-file" => $chr2file{1} || die);
  my $row = $np->next() || die;
  # read one row to initialize headers
  $self->np($np);
  $self->chr2file(\%chr2file);
}

sub find {
  my ($self, %options) = @_;
  if (my $r = $options{"-sj"}) {
    $options{"-chrom"} = $r->{Chr} || die;
    $options{"-base-number"} = $r->{WU_HG19_Pos} || $r->{WU_HG18_Pos} || die;
    $options{"-base-reference"} = $r->{ReferenceAllele} || die;
    $options{"-base-variant"} = $r->{MutantAllele} || die;
  }

  my $chr = cook_chromosome_name($options{"-chrom"} || die "-chrom");
  my $base_number = $options{"-base-number"} || die "-base-number";
  my $base_reference = $options{"-base-reference"};
  my $base_variant = $options{"-base-variant"};

  foreach ($base_reference, $base_variant) {
    my $bad;
    $bad = 1 if $_ and $_ eq "-";
    # query is an indel: won't match
    $bad = 1 if length($_) > 1;
    # MNV: won't match
    return [] if $bad;
  }

  my $nhlbi_file = $self->chr2file->{$chr};
  my @results;
  if ($nhlbi_file) {
    my $parse_callback = sub {
      my ($line) = @_;
      return $self->np->parse_line($line);
    };

    my $dbs = new DelimitedBinarySearch(
					"-file" => $nhlbi_file,
					"-is_headered" => 0,
					"-skip_comment_lines" => 1,
					"-line_parser_callback" => $parse_callback
				       );

    my $hits_raw = $dbs->find(
		       "-comparators" => [
					  {
					   "column_name" => "reference_position",
					   "type" => "number",
					   "value" => $base_number
					   # chromosome
					  }
					 ],
#		       "-verbose" => 1,
		      );
    foreach my $hit (@{$hits_raw}) {
      my $usable = 1;

      if ($base_reference) {
	# reference base given: ensure consistent
#	dump_die($hit, "debug vs user $base_reference", 1);
	if ($hit->{reference_base} ne $base_reference) {
	  dump_die($hit, "ERROR: reference sanity fail vs. user $base_reference", 1);
	  $usable = 0;
	}
      }

      if ($base_variant) {
	# variant base given: report only records containing this base
	my $found;
	foreach my $vb (@{$hit->{variant_bases}}) {
	  $found = 1 if uc($vb) eq uc($base_variant);
	}
#	dump_die($hit);
#	die join ",", @{$hit->{variant_bases}};
	$usable = 0 unless $found;
      }

      push @results, $hit if $usable;
#      dump_die($hit, "hit debug", 1) if $usable;
    }
  } else {
    printf STDERR "no NHLBI data for chr %s\n", $chr;
  }
  return \@results;
}

sub get_parser {
  my ($self) = @_;
  return $self->np;
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
