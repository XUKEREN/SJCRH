#!/bin/env perl
# add annotations from dbNSFP

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
my $TABIX_BATCH_SIZE_DBNSFP = 500;
my $FILE_BATCH_SIZE = 2000;

#Eigen-raw_rankscore
# earlier version

my @COPY_FIELDS = qw(
		      CADD_raw
		      CADD_raw_rankscore
		      CADD_phred
		      DANN_score
		      DANN_rankscore

Eigen_coding_or_noncoding
Eigen-raw
Eigen-phred
Eigen-PC-raw
Eigen-PC-phred
Eigen-PC-raw_rankscore

phastCons100way_vertebrate
phastCons100way_vertebrate_rankscore
phastCons20way_mammalian
phastCons20way_mammalian_rankscore
phyloP100way_vertebrate
phyloP100way_vertebrate_rankscore
phyloP20way_mammalian
phyloP20way_mammalian_rankscore
SiPhy_29way_pi
SiPhy_29way_logOdds
SiPhy_29way_logOdds_rankscore
GERP++_NR
GERP++_RS
GERP++_RS_rankscore
GenoCanyon_score
GenoCanyon_score_rankscore

REVEL_score
REVEL_rankscore

MetaSVM_score
MetaSVM_rankscore
MetaSVM_pred

MetaLR_score
MetaLR_rankscore
MetaLR_pred

VEST3_score
VEST3_rankscore
Transcript_id_VEST3
Transcript_var_VEST3

M-CAP_score
M-CAP_rankscore
M-CAP_pred

MutPred_score
MutPred_rankscore
MutPred_protID
MutPred_AAchange
MutPred_Top5features

integrated_fitCons_score
integrated_fitCons_score_rankscore
integrated_confidence_value
GM12878_fitCons_score
GM12878_fitCons_score_rankscore
GM12878_confidence_value
H1-hESC_fitCons_score
H1-hESC_fitCons_score_rankscore
H1-hESC_confidence_value
HUVEC_fitCons_score
HUVEC_fitCons_score_rankscore
HUVEC_confidence_value

MutationTaster_score
MutationTaster_converted_rankscore
MutationTaster_pred
MutationTaster_model
MutationTaster_AAE

SIFT_score
SIFT_converted_rankscore
SIFT_pred

PROVEAN_score
PROVEAN_converted_rankscore
PROVEAN_pred

Uniprot_acc_Polyphen2
Uniprot_id_Polyphen2
Uniprot_aapos_Polyphen2

Polyphen2_HDIV_score
Polyphen2_HDIV_rankscore
Polyphen2_HDIV_pred
Polyphen2_HVAR_score
Polyphen2_HVAR_rankscore
Polyphen2_HVAR_pred

LRT_score
LRT_converted_rankscore
LRT_pred
LRT_Omega

MutationAssessor_UniprotID
MutationAssessor_variant
MutationAssessor_score
MutationAssessor_score_rankscore
MutationAssessor_pred

FATHMM_score
FATHMM_converted_rankscore
FATHMM_pred
fathmm-MKL_coding_score
fathmm-MKL_coding_rankscore
fathmm-MKL_coding_pred
fathmm-MKL_coding_group

Ensembl_geneid
Ensembl_transcriptid
Ensembl_proteinid


		   );

my %FLAGS;
my @clopts = (
	      "-file=s",
	      "-files=s",

	      "-genome=s",
	      "-dbnsfp=s",
	      "-file-batch-size=i" => \$FILE_BATCH_SIZE,
	      "-tabix-batch-size=i" => \$TABIX_BATCH_SIZE_DBNSFP,
	      VariantParsePolicy::get_command_line_parameters()
	     );

GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $infiles = build_argv_list(
			      "-flags" => \%FLAGS,
			      "-single" => "file",
			      "-set" => "files"
			     );

my $genome = $FLAGS{genome} || die "-genome";
my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
my $TWOBIT = $config_genome->{TWOBIT} || die "no TWOBIT";

my $F_TABIX_POS;
$F_TABIX_POS = $config_genome->{DBNSFP3_TABIX_POS_FIELD} || die;
unless ($FLAGS{dbnsfp}) {
  $FLAGS{dbnsfp} = config2tabix($config_genome, "DBNSFP3_TABIX_DIR");
}

foreach my $infile (@{$infiles}) {
  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			    );

  my $outfile = basename($infile) . ".dbnsfp.tab";

  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => \@COPY_FIELDS
			     );
  my @rows;
  my $flush = sub {
    add_batch_dbnsfp(
		     "-rows" => \@rows,
		    ) if @rows;
    foreach my $row (@rows) {
      $rpt->end_row($row);
    }
    @rows = ();
  };

  while (my $row = $df->get_hash()) {
    push @rows, $row;
    &$flush() if @rows >= $FILE_BATCH_SIZE;
  }
  &$flush();

  $rpt->finish();

}

sub add_batch_dbnsfp {
  my (%options) = @_;
  my $rows = $options{"-rows"} || die "-rows";
  my $start_time = time;
  my $VPP = new VariantParsePolicy("-flags" => \%FLAGS);
  log_message(sprintf "batch dbNSFP start of %d rows", scalar @{$rows});
  my $tabix = get_tabix_dbnsfp();
  my $f_row_key = "_user_row";

  my @query_variants;
  foreach my $r (@{$rows}) {
    my $v = $VPP->get_variant_from_row($r);
    $v->{$f_row_key} = $r;
    push @query_variants, $v;
  }

#  my %map = map {$_, $_} @COPY_FIELDS;

  my $f_tabix_hits = "_tabix_results";

  my $tba = new TabixBatchAnnotation(
				     "-tabix" => $tabix,
				     "-twobit" => $TWOBIT,
				     "-split_count" => $TABIX_BATCH_SIZE_DBNSFP,
				     "-f_tabix_chr" => "chr",
				     "-f_tabix_pos" => $F_TABIX_POS,
				     "-f_tabix_ref_allele" => "ref",
				     "-f_tabix_var_allele" => "alt",

				     "-user_row_key" => $f_row_key,
#				     "-annotation_map" => \%map,
				     "-store_hits" => $f_tabix_hits
				    );

  $tba->query(
	      "-query" => \@query_variants,
	     );

  foreach my $row (@{$rows}) {

    foreach my $f (@COPY_FIELDS) {
      $row->{$f} = "";
    }

    if (my $hits = $row->{$f_tabix_hits}) {
      foreach my $f (@COPY_FIELDS) {
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
  }

  log_message(sprintf "batch dbNSFP annotation for %d rows took %d",
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

