#!/bin/env perl
# add annotations from dbNSFP
# - TO DO: tag equivalence if used

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
use Variant qw(new_variant_from_row);
use TabixBatchAnnotation qw(TABIX_VARIANT);
use VariantParsePolicy;

my $TABIX_INDEL_EQUIVALENCE_DISTANCE = 25;
#my $TABIX_BATCH_SIZE_DBNSFP = 200;
my $TABIX_BATCH_SIZE_EXAC = 500;

my $FREQ_FILTER = 0;
my $FREQ_MAX = .001;
# as in PCGP/medal ceremony
my $FREQ_USE_ALT = 0;
# if set, use AF_Adj rather than AF

my @TABIX_FIELDS = qw(
AC
AN
AC_Adj
AN_Adj
AF
AC_AFR
AN_AFR
AC_AMR
AN_AMR
AC_EAS
AN_EAS
AC_FIN
AN_FIN
AC_NFE
AN_NFE
AC_SAS
AN_SAS
AC_OTH
AN_OTH
		   );

my $F_EXAC_EQUIV = "ExAC_equivalent";
my $F_EXAC_VCF_VARIANT = "ExAC_VCF_variant";
my $F_EXAC_HG19 = "ExAC_original_hg19";

my @POPS = grep {/^AC_/ and not(/Adj/)} @TABIX_FIELDS;
foreach (@POPS) {
  s/^AC_//;
}

my %FLAGS;
my @clopts = (
	      "-file=s",
	      "-files=s",
	      "-genome=s",
	      "-exac-vcf2tab=s",

	      "-require-af-le=f",
	      "-adjusted",
	      "-no-annotation",

	      VariantParsePolicy::get_command_line_parameters(),
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

if (my $f = $FLAGS{"require-af-le"}) {
  $FREQ_FILTER = 1;
  $FREQ_MAX = $f;
}
$FREQ_USE_ALT = 1 if $FLAGS{adjusted};

my $VPP = new VariantParsePolicy("-flags" => \%FLAGS);

my $infiles = build_argv_list(
			      "-flags" => \%FLAGS,
			      "-single" => "file",
			      "-set" => "files"
			     );

my $genome = $FLAGS{genome} || die "-genome";
my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
my $TWOBIT = $config_genome->{TWOBIT} || die "no TWOBIT";

$FLAGS{"exac-vcf2tab"} = config2tabix($config_genome, "EXAC_COOKED_DIR")
  unless $FLAGS{"exac-vcf2tab"};

my $df = new DelimitedFile(
			   "-file" => $FLAGS{"exac-vcf2tab"},
			   "-headers" => 1,
			  );
my $first = $df->get_hash();
my $IS_LIFTOVER;
if (exists $first->{"Chr.orig"}) {
  # if liftover'd (i.e. to hg38), also include original allele info from 37
  my @orig = grep {/\.orig$/} @{$df->headers_raw};
  unshift @TABIX_FIELDS, @orig;
  $IS_LIFTOVER=1;
}


foreach my $infile (@{$infiles}) {
  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			    );
  my $outfile = basename($infile) . ".exac.tab";

  my @extra;
  unless ($FLAGS{"no-annotation"}) {
    @extra = map {"ExAC_" . $_} @TABIX_FIELDS;
    push @extra, map {"ExAC_AF_" . $_,
			"ExAC_rare_" . $_
		      } @POPS;
    push @extra, "ExAC_AF_Adj";
    push @extra, $F_EXAC_EQUIV;
    push @extra, $F_EXAC_VCF_VARIANT;
    push @extra, $F_EXAC_HG19 if $IS_LIFTOVER;
  }

  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => \@extra
			     );

  my $QUEUE_SIZE = 10000;

  my @queue;
  my $total = 0;

  my $flush = sub {
    add_batch_exac(
		   "-rows" => \@queue,
		  );
    foreach my $row (@queue) {
      foreach my $pop (@POPS) {
	my $ac = $row->{"ExAC_AC_" . $pop};
	my $an = $row->{"ExAC_AN_" . $pop};
	my $af = "";
	my $rare = "";
	if (defined $ac and $an) {
	  $af = $ac / $an;
	  $rare = ($af > 0 and $af < 0.001) ? 1 : 0;
	}
	$row->{"ExAC_AF_" . $pop} = $af;
	$row->{"ExAC_rare_" . $pop} = $rare;
      }

      my $ac_adj = $row->{"ExAC_AC_Adj"};
      my $an_adj = $row->{"ExAC_AN_Adj"};
      my $af_adj = "";
      if (defined $ac_adj and $ac_adj =~ /\d/ and
	  defined $an_adj and $an_adj =~ /\d/) {
	$af_adj = $an_adj == 0 ? 0 : $ac_adj / $an_adj;
      }
      $row->{ExAC_AF_Adj} = $af_adj;

      my $usable = 1;
      if ($FREQ_FILTER) {
	if (my $freq = $FREQ_USE_ALT ? $af_adj : $row->{"ExAC_AF"}) {
	  $usable = 0 if $freq > $FREQ_MAX;
#	  dump_die($row) unless $usable;
	}
      }

      $rpt->end_row($row) if $usable;
    }
    $total += @queue;

    @queue = ();
  };

  while (my $row = $df->get_hash()) {
    push @queue, $row;
    &$flush() if @queue >= $QUEUE_SIZE;
  }
  &$flush();

  $rpt->finish();
}

sub add_batch_exac {
  my (%options) = @_;
  my $rows = $options{"-rows"} || die "-rows";

  my $fn = $FLAGS{"exac-vcf2tab"} || die "-exac-vcf2tab";
#  printf STDERR "ExAC: %s\n", $fn;

  my $tabix = get_batch_tabix($fn);

  my $f_row_key = "_user_row";

  my $f_query_snv4 = "_query_snv4";

  my @query_variants;
  foreach my $r (@{$rows}) {
    delete @{$r}{@TABIX_FIELDS};
    # wipe in case present for some reason (shouldn't be)
#    my $v = new_variant_from_row("-row" => $r, "-flags" => \%FLAGS);
    my $v = $VPP->get_variant_from_row($r);
    $v->{$f_row_key} = $r;
    push @query_variants, $v;
    $r->{$f_query_snv4} = $v->get_snv4();
  }

  my %map = map {"ExAC_" . $_ => $_} @TABIX_FIELDS;
  # to do: add leading _ etc. for key safety?

  #
  #  add raw annotation fields to each row:
  #
  my $tba = new TabixBatchAnnotation(
				     "-tabix" => $tabix,
				     "-twobit" => $TWOBIT,
				     "-split_count" => $TABIX_BATCH_SIZE_EXAC,
				     "-f_tabix_chr" => "Chr",
				     "-f_tabix_pos" => "WU_HG19_Pos",
				     "-f_tabix_ref_allele" => "ReferenceAllele",
				     "-f_tabix_var_allele" => "MutantAllele",

				     "-user_row_key" => $f_row_key,
				     "-annotation_map" => \%map,
				    );

  $tba->query(
	      "-query" => \@query_variants,
	      "-tag" => "ExAC",
	     );

  foreach my $r (@{$rows}) {
    my $equiv = "";
    my $hg19 = "";
    my $vcf = "";

    if (my $vt = $r->{TABIX_VARIANT()}) {
      my $snv4_q = $r->{$f_query_snv4} || die;
      my $snv4_t = $vt->get_snv4();
      my $r_tabix = $vt->{_tabix} || die;
      # hack

      if ($snv4_q ne $snv4_t) {
	# if matching was to an equivalent variant, report those details
	$equiv = $snv4_t;
      }

      $vcf = $r_tabix->{variant_raw} || die "no variant_raw";

      if ($IS_LIFTOVER) {
	my $chr = $r_tabix->{"Chr.orig"} || die;
	my $pos = $r_tabix->{"WU_HG19_Pos.orig"} || die;
	my $ra = $r_tabix->{"ReferenceAllele.orig"} || die;
	my $va = $r_tabix->{"MutantAllele.orig"} || die;
	$hg19 = join ".", $chr, $pos, $ra, $va;
      }
    }
    $r->{$F_EXAC_EQUIV} = $equiv;
    $r->{$F_EXAC_HG19} = $hg19;
    $r->{$F_EXAC_VCF_VARIANT} = $vcf;
  }

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

sub get_tabix_dbnsfp {
  my $f_tabix = $FLAGS{dbnsfp} || die;
  printf STDERR "dbNSFP: %s\n", $f_tabix;
  return get_batch_tabix($f_tabix);
}

sub get_batch_tabix {
  my ($f_tabix) = @_;
  return new TabixFile(
		       "-file" => $f_tabix,
		       "-indel_wiggle_bases" => $TABIX_INDEL_EQUIVALENCE_DISTANCE
		      );
}

sub get_variant_from_row {
  my ($row) = @_;
  my $v = new Variant();
  $v->import_bambino_row("-row" => $row, "-postprocessed" => 1);
  return $v;
}
