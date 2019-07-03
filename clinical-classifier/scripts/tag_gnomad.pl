#!/bin/env perl
# add annotations from gnomad
# MNE 5/2017

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list get_hash_option log_message);
use DelimitedFile;
use Reporter;
use TdtConfig;
use TabixFile;
use Variant;
use TabixBatchAnnotation;
use VariantParsePolicy;

my $TABIX_INDEL_EQUIVALENCE_DISTANCE = 25;
#my $TABIX_BATCH_SIZE_DBNSFP = 200;
my $FILE_BATCH_SIZE = 10000;
# rows to read from file before running batch job
my $TABIX_BATCH_SIZE = 500;

my %FLAGS;
my @clopts = (
	      "-file=s",
	      "-files=s",

	      "-genome=s",
	      "-gnomad=s",
	      "-gnomad-cooked=s",

	      VariantParsePolicy::get_command_line_parameters(),

	      "-require-af-lt=f",
	      # drop rows where global AF is < this frequency

	      "-require-af-lt-pop=s",
	      # same for a specific sub-population,
	      # format POP,freq, e.g. AFR,.01

	      "-batch-size=i" => \$TABIX_BATCH_SIZE,

	      "-no-annotation",

	      "-compress",

	      "-out=s",
	      "-verbose",

	      "-equivalence-distance=i" => \$TABIX_INDEL_EQUIVALENCE_DISTANCE,

	     );

GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $VPP = new VariantParsePolicy("-flags" => \%FLAGS);
my $infiles = build_argv_list(
			      "-flags" => \%FLAGS,
			      "-single" => "file",
			      "-set" => "files"
			     );

my $genome = $FLAGS{genome} || die "-genome";
my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
my $TWOBIT = $config_genome->{TWOBIT} || die "no TWOBIT";

unless ($FLAGS{gnomad}) {
  $FLAGS{gnomad} = config2tabix($config_genome, "GNOMAD_VCF_DIR");
}

my $require_af_lt = $FLAGS{"require-af-lt"};

my $require_af_lt_pop = $FLAGS{"require-af-lt-pop"};
my $require_af_lt_pop_name;
my $require_af_lt_pop_af;
if ($require_af_lt_pop) {
  my @f = split /,/, $require_af_lt_pop;
  die "specify POPULATION_CODE,FREQ: e.g. AFR,.01\n" unless @f == 2;
  ($require_af_lt_pop_name, $require_af_lt_pop_af) = @f;
  die "needs testing after gnomad format change";
}

die "can't use multiple infiles with -out" if @{$infiles} > 1;

my $GNOMAD_FIELDS = get_gnomad_fields();

foreach my $infile (@{$infiles}) {
  my $outfile = $FLAGS{out} || basename($infile) . ".gnomad.tab";

  my @extra = map {"gnomad_" . $_} @{$GNOMAD_FIELDS};
  push @extra, "gnomad_hit";
  push @extra, "gnomad_equiv";
  @extra = () if $FLAGS{"no-annotation"};

  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			    );
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => \@extra,
			      "-compress" => $FLAGS{compress},
			     );
  my @queue;
  my $flush = sub {
    printf STDERR "flushing batch of %d...\n", scalar @queue;
    add_batch_gnomad(
		     "-rows" => \@queue,
		    );
    foreach my $row (@queue) {
      foreach my $f (@{$GNOMAD_FIELDS}) {
	$row->{"gnomad_" . $f} = $row->{$f};
      }

      my $usable = 1;

      if (my $af = $row->{AF}) {
	# AF only defined if a match
#	dump_die($row, "WTF: . for allele freq") if $af eq ".";

	my $reject_reason;

	if ($require_af_lt and $af ne "." and not($af < $require_af_lt)) {
	  $usable = 0;
	  $reject_reason = "AF_cutoff $af";
	}
	# can be "." if user reference allele is wrong

	if ($require_af_lt_pop_name) {
	  my $field = "gnomad_AF_" . $require_af_lt_pop_name;
	  dump_die($row, "no column $field" unless exists $row->{$field};
	  my $paf = $row->{$field} || 0;
	  if (not($paf < $require_af_lt_pop_af)) {
	    $usable = 0;
	    $reject_reason = "pop_cutoff $paf";
	  }
	}

	if ($FLAGS{verbose}) {
	  printf STDERR "reject %s: %s\n", $VPP->get_variant_from_row($row)->get_snv4(), $reject_reason unless $usable;
	}
      }

#      printf "usable=%s\n", $usable;

      $rpt->end_row($row) if $usable;
    }
    @queue = ();
  };

  while (my $row = $df->get_hash()) {
    push @queue, $row;
    &$flush() if @queue >= $FILE_BATCH_SIZE;
  }
  &$flush();

  $rpt->finish();

}

sub add_batch_gnomad {
  my (%options) = @_;
  my $rows = $options{"-rows"} || die "-rows";
  my $start_time = time;
  log_message("batch gnomad start");
  my $tabix = get_tabix_gnomad();
  my $f_row_key = "_user_row";

  my @query_variants;
  foreach my $r (@{$rows}) {
    my $v = $VPP->get_variant_from_row($r);
    $v->{$f_row_key} = $r;
    push @query_variants, $v;
  }

  my $f_tabix_hits = "_tabix_results";

  my @opts;
  if ($FLAGS{"gnomad-cooked"}) {
    @opts = (
	     "-f_tabix_chr" => "Chr",
	     "-f_tabix_pos" => "Pos",
	     "-f_tabix_ref_allele" => "Ref",
	     "-f_tabix_var_allele" => "Alt",
	    );
    die "needs testing for update";
  } else {
    @opts = ("-vcf2tab" => "-gnomad");
    # "-no-insertion-adjustment" not needed because TabixFile already uses it
  }

  my $tba = new TabixBatchAnnotation(
				     "-tabix" => $tabix,
				     "-twobit" => $TWOBIT,
				     "-split_count" => $TABIX_BATCH_SIZE,

				     "-user_row_key" => $f_row_key,
#				     "-annotation_map" => \%map,
				     "-store_hits" => $f_tabix_hits,
				     @opts
				    );

  $tba->query(
	      "-query" => \@query_variants,
	     );

  foreach my $row (@{$rows}) {

    foreach my $f (@{$GNOMAD_FIELDS}) {
      $row->{$f} = "";
    }

    my @tabix_hits;
    my @tabix_equiv;
    if (my $hits = $row->{$f_tabix_hits}) {
      my $snv4_q = $VPP->get_variant_from_row($row)->get_snv4();
      foreach my $hit (@{$hits}) {
	my $snv4_t = $tba->get_snv4($hit);
	push @tabix_hits, $snv4_t;
	push @tabix_equiv, $snv4_t if $snv4_t ne $snv4_q;
      }


      foreach my $f (@{$GNOMAD_FIELDS}) {
	my @v;
	my %saw;
	foreach my $hit (@{$hits}) {
	  my $v = $hit->{$f};
	  die "no data for $f" unless defined $v;
	  unless ($saw{$v}) {
	    push @v, $v;
	    $saw{$v} = 1;
	  }
	}
#	dump_die($row, "FIX ME: multiple values found for $f") if @v > 1;
	$row->{$f} = join "|", @v;
      }
    }
    $row->{gnomad_hit} = join ",", @tabix_hits;
    $row->{gnomad_equiv} = join ",", @tabix_equiv;
  }

  log_message(sprintf "batch gnomad annotation for %d rows took %d",
	  scalar(@{$rows}), time - $start_time);

}


sub config2tabix {
  # find tabix file from a genome config directory (i.e. tartan output)
  my ($config_genome, $cv) = @_;
  my $tf;
  if (my $dir = $config_genome->{$cv}) {
    if (-d $dir) {
      my @gz = glob($dir . "/*.gz");
      if (@gz == 1) {
	($tf) = @gz;
	my $tbi = $tf . ".tbi";
	die "where is $tbi" unless -s $tbi;
      } else {
	confess "ERROR: not exactly one .gz file in $dir for $cv";
      }
    }
  }
  return $tf;
}

sub get_gnomad_file {
  return $FLAGS{gnomad} || $FLAGS{"gnomad-cooked"} || die "-gnomad/-gnomad-cooked";
}

sub get_tabix_gnomad {
  my $f_tabix = get_gnomad_file();
  printf STDERR "gnomad: %s\n", $f_tabix;
  return get_batch_tabix($f_tabix);
}

sub get_batch_tabix {
  my ($f_tabix) = @_;
  return new TabixFile(
		       "-file" => $f_tabix,
		       "-indel_wiggle_bases" => $TABIX_INDEL_EQUIVALENCE_DISTANCE
		      );
}

sub get_gnomad_fields {
  my $f_gnomad = get_gnomad_file();

  my %exclude = map {$_, 1} qw(
				variant_id
				variant_length
				Chr
				Pos
				Chr_Allele
				Alternative_Allele
				additional_shift_bases
				variant_note
				gnomad_variant_type2
			     );
  # keep:
  #   - variant_type: a gnomad field requested by Zhaoming,
  #     which clobbers the vcf2tab value.
  #   - variant_raw: for raw gnomAD lookup
  #

  open(H, "vcf2tab.pl -file $f_gnomad -gnomad -stdout|") || die;
  my $first = <H>;
  chomp $first;
  my @f = split /\t/, $first;
  my @h;
  my %saw;

  foreach (@f) {
    next if $saw{$_};
    next if $exclude{$_};
    push @h, $_;
    $saw{$_} = 1;
  }

  return \@h;
}
