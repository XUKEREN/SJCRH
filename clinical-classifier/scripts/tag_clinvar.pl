#!/bin/env perl
# standalone clinvar tagging utility
# 2 versions:
# - cooked version used by medal ceremony
# - VCF version that runs vcf2tab.pl on the fly
#
# TO DO: update to use VariantParsePolicy.pm

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
use Variant;
use TabixBatchAnnotation;

use constant INSERTION_BASE_BEFORE => "before";
use constant INSERTION_BASE_AFTER => "after";

use constant TABIX_BATCH_SIZE_CLINVAR => 1000;

use constant CLNSIG_LP => 4;
use constant CLNSIG_P => 5;
# see clinvar vcf distribution

my $QUEUE_SIZE = 1000;
my $TABIX_INDEL_EQUIVALENCE_DISTANCE = 25;
my $FIELD_CLINVAR_TABIX = "__clinvar_tabix";

my $F_OUT_CLINVAR = "ClinVar";
my $F_OUT_VARIANT = "ClinVar_variant";
my $F_OUT_CLNSIG = "ClinVar_CLNSIG";
my $F_OUT_CLNSIG_DESC = "ClinVar_CLNSIG_desc";
my $F_OUT_DBSNP = "ClinVar_dbSNP";
my $F_OUT_PUBMED = "PubMed";
my $F_OUT_GOLD_STARS = "ClinVar_gold_stars";
my $F_OUT_VID = "ClinVar_Variation_ID";

my $F_CHR;
my $F_POS;
my $F_RA;
my $F_VA;

my $F_TABIX_CHR = "Chr";
my $F_TABIX_POS = "WU_HG19_Pos";
my $F_TABIX_RA = "ReferenceAllele";
my $F_TABIX_VA = "MutantAllele";

my %FLAGS;
my @clopts = (
	      "-genome=s",
	      "-file=s",
	      "-tabix-clinvar=s",
	      # tabix-converted version,
	      # see clinical-classifier/clinvar2tabix.pl
	      "-f-chr=s" => \$F_CHR,
	      "-f-pos=s" => \$F_POS,
	      "-f-ra=s" => \$F_RA,
	      "-f-va=s" => \$F_VA,

	      "-sj-post",
	      "-bambino",
	      "-insertion-base=s",

	      "-p-or-lp",
	      # only annotate if P or LP in ClinVar

	      "-vcf=s",
	     );

GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $genome = $FLAGS{genome} || die "-genome";
my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
my $TWOBIT = $config_genome->{TWOBIT} || die "no TWOBIT";

if ($FLAGS{"sj-post"}) {
  $F_CHR = "Chr";
  $F_POS = "WU_HG19_Pos";
  $F_RA = "ReferenceAllele";
  $F_VA = "MutantAllele";
} elsif ($FLAGS{"bambino"}) {
  $F_CHR = "Chr";
  $F_POS = "Pos";
  $F_RA = "Chr_Allele";
  $F_VA = "Alternative_Allele";
  $FLAGS{"insertion-base"} = "after";
}

my $insertion_base = $FLAGS{"insertion-base"} || die "-insertion-base [before|after]";
# die $insertion_base
die "-insertion-base must be 'before' or 'after'" unless $insertion_base eq INSERTION_BASE_BEFORE or $insertion_base eq INSERTION_BASE_AFTER;

die "fields not defined" unless $F_CHR and $F_POS and $F_RA and $F_VA;


my $infile = $FLAGS{file} || die "-file";

my $outfile = basename($infile) . ".clinvar.tab";

my $df = new DelimitedFile(
			   "-file" => $infile,
			   "-headers" => 1,
			  );

my @OF;
if ($FLAGS{vcf}) {
  @OF = qw(
	    variant_raw
	    CLNSIG
	    CLNSIG_desc
	    CLNSIG_full
	    CLNORIGIN
	    CLNDSDB
	    CLNDSDBID
	    CLNDBN
	    CLNREVSTAT
	    dbSNP
	 );
  push @OF, $F_OUT_PUBMED;
} else {
  @OF = (
	 $F_OUT_CLINVAR,
	 $F_OUT_VARIANT,
	 $F_OUT_CLNSIG,
	 $F_OUT_CLNSIG_DESC,
	 $F_OUT_DBSNP,
	 $F_OUT_PUBMED,
	 $F_OUT_GOLD_STARS,
	 $F_OUT_VID,
	);
}

my $rpt = $df->get_reporter(
			    "-file" => $outfile,
			    "-extra" => \@OF,
			    "-clobber" => 1,
			   );
my @queue;
my $flush = sub {
  if (@queue) {
    if ($FLAGS{vcf}) {
      add_batch_clinvar_vcf(
			"-rows" => \@queue,
		       );
    } else {
      add_batch_clinvar(
			"-rows" => \@queue,
		       );
    }
    foreach my $row (@queue) {
      $rpt->end_row($row);
    }
  }
  @queue = ();
};

while (my $row = $df->get_hash()) {
  foreach my $f ($F_RA, $F_VA) {
    my $v = $row->{$f};
    $row->{$f} = "-" if $v eq "" or $v eq " ";
  }

  push @queue, $row;
  &$flush() if @queue >= $QUEUE_SIZE;
}
&$flush();

$rpt->finish();

sub get_tabix_clinvar {
  my $f_tabix = $FLAGS{vcf} || $FLAGS{"tabix-clinvar"} || die "-vcf | -tabix-clinvar";
  printf STDERR "ClinVar: %s\n", $f_tabix;
  return get_batch_tabix($f_tabix);
}

sub get_batch_tabix {
  my ($f_tabix) = @_;
  return new TabixFile(
		       "-file" => $f_tabix,
		       "-indel_wiggle_bases" => $TABIX_INDEL_EQUIVALENCE_DISTANCE
		      );
}


sub add_batch_clinvar {
  my (%options) = @_;

  unless ($FLAGS{"tabix-clinvar"}) {
    my $dir = $config_genome->{CLINVAR_TABIX_DIR} || die;
    my ($fn) = glob($dir . "/*gz");
    die "can't get tabix clinvar from genome config CLINVAR_TABIX_DIR" unless $fn;
    $FLAGS{"tabix-clinvar"} = $fn;
  }

  my $rows = $options{"-rows"} || die "-rows";
  my $start_time = time;
  my $label = "ClinVar";

  log_message("batch $label start");
  my $tabix = get_tabix_clinvar();
  my $f_row_key = "_user_row";

  my @query_variants;
  foreach my $r (@{$rows}) {
    my $v = get_variant_from_row($r);
    $v->{$f_row_key} = $r;
    push @query_variants, $v;
  }

  my $tba = new TabixBatchAnnotation(
				     "-tabix" => $tabix,
				     "-twobit" => $TWOBIT,
				     "-split_count" => TABIX_BATCH_SIZE_CLINVAR,
				     "-f_tabix_chr" => $F_TABIX_CHR,
				     "-f_tabix_pos" => $F_TABIX_POS,
				     "-f_tabix_ref_allele" => $F_TABIX_RA,
				     "-f_tabix_var_allele" => $F_TABIX_VA,

				     "-user_row_key" => $f_row_key,

				     "-store_hits" => $FIELD_CLINVAR_TABIX,
#				     "-store_site" => $FIELD_CLINVAR_TABIX_SITE,
				    );

  $tba->query(
	      "-query" => \@query_variants,
	     );

  my $p_or_lp = $FLAGS{"p-or-lp"};

  foreach my $row (@{$rows}) {
    my $v_clinvar = 0;
    my $v_variant = "";
    my $v_clnsig = "";
    my $v_clnsig_desc = "";
    my $v_dbsnp = "";
    my $v_pmid = "";
    my $v_gold_stars = "";
    my $v_vid = "";

    if (my $hits = $row->{$FIELD_CLINVAR_TABIX}) {
      my @variants;
      my @clnsig;
      my @clnsig_desc;
      my @dbsnp;
      my @pmid;
      my @gold_stars;
      my @vid;
      foreach my $hit (@{$hits}) {
	if (my $vr = $hit->{variant_raw}) {
	  # hit came from VCF parser, provided by vcf2tab.pl
	  push @variants, $hit->{variant_raw};
	} else {
	  # reconstruct from tabix hit
	  push @variants, join ".", @{$hit}{$F_TABIX_CHR, $F_TABIX_POS, $F_TABIX_RA, $F_TABIX_VA};
	}
#	dump_die($hit);
	push @clnsig, $hit->{CLNSIG};
	push @clnsig_desc, $hit->{CLNSIG_desc};
	push @dbsnp, $hit->{dbSNP} || "";
	push @pmid, $hit->{PubMed} if $hit->{PubMed};
	push @gold_stars, $hit->{ClinVar_gold_stars};
	push @vid, $hit->{ClinVar_Variation_ID};
      }

      my $usable;
      if ($p_or_lp) {
	$usable = 1 if grep {$_ eq CLNSIG_LP or $_ eq CLNSIG_P} @clnsig;
      } else {
	$usable = 1;
      }

      if ($usable) {
	$v_clinvar = 1;
	$v_clnsig = join "|", @clnsig;
	$v_clnsig_desc = join "|", @clnsig_desc;
	$v_dbsnp = join "|", @dbsnp;
	$v_variant = join "|", @variants;
	$v_pmid = join ",", @pmid;
	$v_gold_stars = join ",", @gold_stars;
	$v_vid = join ",", @vid;
      }
    }

    $row->{$F_OUT_CLINVAR} = $v_clinvar;
    $row->{$F_OUT_VARIANT} = $v_variant;
    $row->{$F_OUT_CLNSIG} = $v_clnsig;
    $row->{$F_OUT_CLNSIG_DESC} = $v_clnsig_desc;
    $row->{$F_OUT_DBSNP} = $v_dbsnp;
    $row->{$F_OUT_PUBMED} = $v_pmid;
    $row->{$F_OUT_GOLD_STARS} = $v_gold_stars;
    $row->{$F_OUT_VID} = $v_vid;
  }

  log_message(sprintf "batch %s annotation for %d rows took %d",
	  $label, scalar(@{$rows}), time - $start_time);


}


sub get_variant_from_row {
  my ($row) = @_;
  my $chr = $row->{$F_CHR} || dump_die($row, "no $F_CHR");
  my $pos = $row->{$F_POS} || dump_die($row, "no $F_POS");
  my $ra = $row->{$F_RA};
  my $va = $row->{$F_VA};
  foreach ($ra, $va) {
    dump_die($row, "undef reference or variant allele") unless defined $_;
    $_ = "-" unless $_;
  }

  my $v = new Variant();
  $v->import_generic(
		     "-reference-name" => $chr,
		     "-base-number" => $pos,
		     "-reference-allele" => $ra,
		     "-variant-allele" => $va
		    );
  if ($v->is_insertion() and $insertion_base eq INSERTION_BASE_AFTER) {
    # adjust Bambino-style base numbering to base before insertion
    $v->start($v->start - 1);
    $v->end($v->end - 1);
  }
  return $v;
}


sub add_batch_clinvar_vcf {
  my (%options) = @_;
  my $rows = $options{"-rows"} || die "-rows";
  my $start_time = time;
  my $label = "ClinVar";

  log_message("batch $label start");
  my $tabix = get_tabix_clinvar();
  my $f_row_key = "_user_row";

  my @query_variants;
  foreach my $r (@{$rows}) {
    my $v = get_variant_from_row($r);
    $v->{$f_row_key} = $r;
    push @query_variants, $v;
  }

  my $tba = new TabixBatchAnnotation(
				     "-tabix" => $tabix,
				     "-twobit" => $TWOBIT,
				     "-split_count" => TABIX_BATCH_SIZE_CLINVAR,
				     "-vcf2tab" => "-clinvar 2",

				     "-user_row_key" => $f_row_key,

				     "-store_hits" => $FIELD_CLINVAR_TABIX,
#				     "-store_site" => $FIELD_CLINVAR_TABIX_SITE,
				    );

  # -no-insertion-adjustment isn't needed because TabixFile.pm already uses it

  $tba->query(
	      "-query" => \@query_variants,
	     );

  die "p-or-lp not implemented yet in VCF mode" if $FLAGS{"p-or-lp"};

  foreach my $row (@{$rows}) {
    my $v_variant = "";
    my $v_clnsig = "";
    my $v_dbsnp = "";

    my %info;
    if (my $hits = $row->{$FIELD_CLINVAR_TABIX}) {
      my @variants;
      my @clnsig;
      my @dbsnp;
      foreach my $hit (@{$hits}) {
#	dump_die($hit);
	foreach my $f (@OF) {
	  push @{$info{$f}}, $hit->{$f} if defined $hit->{$f};
	}
      }
    }

    foreach my $f (@OF) {
      if ($info{$f}) {
	$row->{$f} = join ",", @{$info{$f}};
      } else {
	$row->{$f} = "";
      }
    }

  }

  log_message(sprintf "batch %s annotation for %d rows took %d",
	  $label, scalar(@{$rows}), time - $start_time);


}
