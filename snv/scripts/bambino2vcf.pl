#!/bin/env perl
# convert Bambino->VCF
# MNE 5/2015
# TO DO:
# - streaming option for very large files??
# - prune unused columns from rows during loading to save RAM?

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use File::Basename;
use FileHandle;

use TdtConfig;
use MiscUtils qw(dump_die);
use FileUtils qw(universal_open);
use DelimitedFile;
use Reporter;
use WorkingFile;
use Variant;
use FAI;
use LabelDigest;
use VariantParsePolicy;

dump_die(\%INC) if $ENV{DEBUG_INC};

my $OUT_QUEUE_FLUSH_SIZE = 1000;
my $OUT_QUEUE_FLUSH_WRITE_COUNT = 250;

# spec 1.4.1:
use constant VCF_MISSING_VALUE => ".";
use constant VCF_FILTER_PASS => "PASS";
use constant VCF_INFO_DELIMITER => ";";
# 1.4.1

use constant INFO_BAMBINO_KEY => "BKY";

use constant INFO_PUBMED_IDS => "PMIDS";
use constant INFO_SOMATIC_STATUS => "SOM";
use constant INFO_GENOME_WIDE_SCREEN => "GWS";
use constant INFO_SAMPLE_ID => "SID";
use constant INFO_VARIANT_NAME => "VN";
use constant INFO_VARIANT_CLASS => "VC";

use constant INFO_CANCER_TYPE => "CT";
use constant INFO_CANCER_TYPE_FULL => "CTF";
use constant INFO_CANCER_SUBTYPE => "ST";
use constant INFO_CANCER_SUBTYPE_FULL => "STF";
use constant INFO_CANCER_SUBGROUP => "SG";
use constant INFO_CANCER_SUBGROUP_FULL => "SGF";
use constant INFO_DATASET => "DS";
use constant INFO_DISEASE_PHASE => "PH";
# Protein Paint db sample table

my %INFO_STRING = (
		   INFO_SOMATIC_STATUS() => "somatic status",
		   INFO_GENOME_WIDE_SCREEN() => "genome-wide screen (y or n)",
		   INFO_CANCER_TYPE() => "cancer type",
		   INFO_CANCER_TYPE_FULL() => "cancer type (full)",
		   INFO_CANCER_SUBTYPE() => "cancer subtype",
		   INFO_CANCER_SUBTYPE_FULL() => "cancer subtype (full)",
		   INFO_CANCER_SUBGROUP() => "cancer subgroup",
		   INFO_CANCER_SUBGROUP_FULL() => "cancer subgroup (full)",
		   INFO_SAMPLE_ID() => "sample ID",
		   INFO_VARIANT_NAME() => "variant name",
		   INFO_VARIANT_CLASS() => "variant class",
		   INFO_DATASET() => "data set",
		   INFO_DISEASE_PHASE() => "disease phase",
		  );

my %FLAGS;

my $NEED_BAMBINO_SORT = 0;
# raw Bambino files are only "mostly" sorted

my $INCLUDE_SAMPLE_IDS = 1;
# not sure if sharing SJXXX IDs externally is OK?

my $REMOVE_DUPLICATES = 1;

my $SKIP_INCOMPLETE = 1;
# if required data not populated

my $PP_SNV_INDEL = "/research/dept/compbio/common/proteinpaint-dev/tp/anno/db/pg_dmp_portal_prod_snvindel.tsv.gz";
my $PP_SAMPLE = "/research/dept/compbio/common/proteinpaint-dev/tp/anno/db/pg_dmp_portal_prod_sample_master.tsv.gz";

my @EXCLUDE_HEADERS;

GetOptions(\%FLAGS,
	   "-file=s",
	   "-out=s",
	   "-stdout",

	   "-genome=s",
	   "-mt",
	   # rename M to MT for compatibility w/downstream processing

	   "-sync-ref-names",
	   "-remove-duplicates=i" => \$REMOVE_DUPLICATES,

	   "-tabix",

	   "-pp",
	   "-pmid=i",

	   "-presort",
	   "-bambino-sort=i" => \$NEED_BAMBINO_SORT,

	   "-annot",
	   # pass through all additional columns as supplemental annotations
	   "-annot-label-len=i",

	   VariantParsePolicy::get_command_line_parameters(),

	   "-vcf",
	   # file contains untranslated alleles and positions (i.e.
	   # indels have padding base[s]), pass through as is

	   "-complex-ok",
	   "-fasta=s",

	   "-exclude=s" => \@EXCLUDE_HEADERS,
	  );


my $VPP = new VariantParsePolicy("-flags" => \%FLAGS);

my $VCF_MODE = $FLAGS{vcf};
if ($VCF_MODE) {

  foreach my $p (qw(insertion-base bambino)) {
    die sprintf "ERROR: the \"-%s\" parameter cannot be used when using \"-vcf\"\n", $p if exists $FLAGS{$p};
  }

  $FLAGS{"insertion-base"} = "before";
}
$VPP->setup_check();

my @pp_snvindel_headers = qw(
			   donor_name
			   donor_id
			   sample_name
			   cell_sample_id
			   sample_disease_phase
			   analyte_sample_id
			   analyte
			   variant_origin
			   variant_origin_pp
			   validation_status
			   quality
			   analysis_run_id
			   analysis_run_name
			   result_id_case
			   sample_id_case
			   sample_id_control
			   genotype_id
			   genotype_assay_id_case
			   allele_1_id_case
			   alelle_2_id_case
			   allele_1_signal_case
			   allele_2_signal_case
			   gv_id
			   variant_name
			   variant_class
			   variant_type
			   amino_acid_change
			   locus_id
			   chromosome
			   chr_position
			   gene_symbol
			   isoform_accession
			   protein_gi
			   gv2gene_id
			   analysis_type_name
			   analysis_type_type
			   variant_analysis
			   seq_src
			   result_id_ctrl
			   genotype_assay_id_ctrl
			   allele_1_id_ctrl
			   allele_2_id_ctrl
			   allele_1_signal_ctrl
			   allele_2_signal_crtl
			   allele_1
			   allele_2
			   allele_1_is_reference
			   allele_2_is_reference
			   mutant_reads_in_case
			   total_reads_in_case
			   mutant_reads_in_control
			   total_reads_in_control
			   src
			   donor_name_orig
			   sample_name_orig
			   cdna_coordinate
			   pubmed_id_list
			   dataset_label
			   allele_1_signal_rna
			   allele_2_signal_rna
			   rnaseq_maf
			   sample_has_loh_results
			   loh_seg_mean
			   loh
			   committee_classification
			   );

my @pp_sample_headers = qw(
			   donor_name
			   sample_name
			   sample_type
			   diagnosis_group_short
			   diagnosis_group_full
			   diagnosis_short
			   diagnosis_full
			   diagnosis_subtype_short
			   diagnosis_subtype_full
			   diagnosis_subgroup_short
			   diagnosis_subgroup_full
			   disease_code_short
			   dataset_label
			);

bambino2vcf();

sub bambino2vcf {
  my $fasta = $FLAGS{fasta};
  my $genome;
  unless ($fasta) {
    $genome = $FLAGS{genome} || die "specify -fasta FILE (must have .fai index) or -genome GENOME\n";
    my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
    $fasta = $config_genome->{FASTA} || die;
  }
#  printf STDERR "FASTA: %s\n", $fasta;
  my $fai = new FAI("-fasta" => $fasta);

  my $df;
  my $pp_mode = $FLAGS{pp};
  die "Protein Paint mode needs work for VariantParsePolicy compatibility" if $pp_mode;
  # not sure if this will be run again

  $REMOVE_DUPLICATES = 0 if $pp_mode and $INCLUDE_SAMPLE_IDS;
  my %samples;
  my $annot_mode = $FLAGS{annot};
  my $pmid_filter = $FLAGS{pmid};
  my $infile;

  if ($pp_mode) {
    $infile = $PP_SNV_INDEL;
    $df = new DelimitedFile(
			    "-file" => $PP_SNV_INDEL,
			    "-headers" => 0
			   );
    $df->headers(1);
    $df->headers_raw(\@pp_snvindel_headers);
    # unheadered, add manually (DANGEROUS: headered file or db query safer)

    my $df_s = new DelimitedFile(
				 "-file" => $PP_SAMPLE,
				 "-headers" => 0
				);
    $df_s->headers(1);
    $df_s->headers_raw(\@pp_sample_headers);
    # unheadered, add manually (DANGEROUS: headered file or db query safer)
    while (my $r = $df_s->get_hash()) {
      $samples{$r->{sample_name}} = $r;
    }
  } else {
    $infile = $FLAGS{file} || die "-file";

    if ($NEED_BAMBINO_SORT) {
      # for raw Bambino-format files ONLY
      my $cmd = sprintf 'bambino_sort.pl -in %s -stdout', $infile;
      my $fh = new FileHandle();
      $fh->open($cmd . "|") || die;
      $df = new DelimitedFile("-fh" => $fh,
			      "-headers" => 1,
			     );
    } else {
      $df = new DelimitedFile("-file" => $infile,
			      "-headers" => 1,
			     );
    }
  }

  my @input_rows;
  my $sort_mode = $FLAGS{presort};
  if ($sort_mode) {
    if ($pp_mode) {
      die "presort needs fixing for VariantParsePolicy";
    }
    my @rows;
    printf STDERR "presort: loading...";
    my $f_chr = $VPP->f_chr() || die;
    my $f_pos = $VPP->f_pos() || die;
    while (my $row = $df->get_hash()) {
      $row->{$f_chr} =~ s/^chr//i;
      push @rows, $row;
    }


    printf STDERR "sorting...";
    @input_rows = sort {$a->{$f_chr} cmp $b->{$f_chr} ||
			  $a->{$f_pos} <=> $b->{$f_pos}} @rows;
    printf STDERR "done\n";
  }

  my $rnm;
  my $mt_mode = $FLAGS{mt};
  my $sync_names = $FLAGS{"sync-ref-names"};

  my $outfile = $FLAGS{out} || basename($infile) . ".vcf";
  my $fh;
  my $wf;
  if ($FLAGS{stdout} or $outfile eq "-") {
    $fh = *main::STDOUT;
  } else {
    $wf = new WorkingFile($outfile);
    if ($FLAGS{tabix}) {
      $wf->compress("bgzip");
      $wf->tabix("vcf");
    }
    $fh = $wf->output_filehandle();
  }

  my @annot_info_h;
  my %annot_info2label;

  print $fh "##fileformat=VCFv4.1\n";
  printf $fh "##genome=%s\n", $genome if $genome;
  printf $fh '##INFO=<ID=%s,Number=1,Type=String,Description="Bambino lookup key">' . "\n", INFO_BAMBINO_KEY;
  if ($pp_mode) {
    printf $fh '##INFO=<ID=%s,Number=.,Type=String,Description="comma-delimited list of PubMed IDs">' . "\n", INFO_PUBMED_IDS;

    foreach my $key (sort keys %INFO_STRING) {
      next if $key eq INFO_SAMPLE_ID and not($INCLUDE_SAMPLE_IDS);
      printf $fh '##INFO=<ID=%s,Number=1,Type=String,Description="%s">' . "\n", $key, $INFO_STRING{$key};
    }
  } elsif ($annot_mode) {
    my %ignore = map {$_, 1} (
			      $VPP->f_chr,
			      $VPP->f_pos,
			      $VPP->f_ra,
			      $VPP->f_va
			     );
    foreach (@EXCLUDE_HEADERS) {
      $ignore{$_} = 1;
    }

    my $ld;
    if (my $alen = $FLAGS{"annot-label-len"}) {
      $ld = new LabelDigest("-max_length" => $alen);
    }

    foreach my $h (@{$df->headers_raw}) {
      next if $ignore{$h};
      push @annot_info_h, $h;
      my $label;
      if ($ld) {
	$label = $ld->get_brief_label("-label" => $h);
      } else {
	$label = $h;
      }
      $annot_info2label{$h} = $label;
      # TO DO: option for shorter labels;
      # however, will this save much space if output compressed?
    }

    foreach my $h (@annot_info_h) {
      printf $fh '##INFO=<ID=%s,Number=1,Type=String,Description="%s">' . "\n",$annot_info2label{$h}, $h;
    }

  }

  my $rpt = new Reporter(
			 "-fh" => $fh,
			 "-delimiter" => "\t",
			 "-labels" => [
				       "#CHROM",
				       "POS",
				       "ID",
				       "REF",
				       "ALT",
				       "QUAL",
				       "FILTER",
				       "INFO",
				      ]
			    );

  my %saw;
  my $duplicates = 0;
  my $malformed = 0;

  my @out_queue;

  my $flush = sub {
    my (%options) = @_;
    my $force = $options{"-force"};
    my $verbose = 0;
    printf STDERR "flush start: force=%d\n", $force if $verbose;
    die "-force [0|1]" unless defined $force;
    my @sorted = sort {$a->{POS} <=> $b->{POS}} @out_queue;

    my $flush_count = $force ? scalar @sorted : $OUT_QUEUE_FLUSH_WRITE_COUNT;
    $flush_count = @sorted if $flush_count > @sorted;
    for (my $i = 0; $i < $flush_count; $i++) {
      # flush the lowest sorted values
      my $r = shift @sorted;
      printf STDERR "write %s %d\n", $r->{"#CHROM"}, $r->{POS} if $verbose;
      $rpt->end_row($r);
    }
    @out_queue = @sorted;
    # leave remainder for next pass
    print STDERR "flush done\n" if $verbose;
  };

  my $last_rn = "";

  while (1) {
    my $row;
    if ($sort_mode) {
      last unless @input_rows;
      $row = shift @input_rows;
    } else {
      $row = $df->get_hash() || last;
    }

    my ($chr, $pos, $ra, $va);
    my %info;

    if ($pp_mode) {
      if ($pmid_filter) {
	my %pmid = map {$_, 1} split /,/, $row->{pubmed_id_list};
	next unless $pmid{$pmid_filter};
      }

      $chr = $row->{chromosome} || die;
      $pos = $row->{chr_position} || die;
      $ra = $row->{allele_1} || die;
      $va = $row->{allele_2} || die;
      foreach ($ra, $va) {
	if (/\-/) {
	  die if /[acgt]/i;
	  $_ = "";
	  # convert to bambino style: blank for indels
	}
      }
#      $row->{$F_CHR} = $chr;
#      $row->{$F_POS} = $pos;
#      $row->{$F_RA} = $ra;
#      $row->{$F_VA} = $va;
      die "VariantParsePolicy update";

      $info{INFO_PUBMED_IDS()} = $row->{pubmed_id_list};

      my $som = "unknown";
      my $vstat = $row->{validation_status} || die;
      my $vorg = $row->{variant_origin} || die;
      if (lc($vstat) eq "valid") {
	if (lc($vorg) eq "somatic") {
	  $som = "Confirmed somatic";
	  # COSMIC style
	}
      } else {
	dump_die($row, "not validated variant");
      }
      $info{INFO_SOMATIC_STATUS()} = $som;
      # Mutation somatic status: Confirmed somatic

      my $gws;
      my $src = $row->{seq_src} || die;

      if ($src eq "{NEXT_GEN_WGS}" or $src eq "{NEXT_GEN_EXCAP}") {
	$info{INFO_GENOME_WIDE_SCREEN()} = "y";
	# all variants are either WGS or WES
      } else {
	dump_die($row, "unhandled seq_src");
      }

      my $sid = $row->{sample_name} || die;

      $info{INFO_SAMPLE_ID()} = $sid;
      $info{INFO_VARIANT_NAME()} = $row->{variant_name};
      $info{INFO_VARIANT_CLASS()} = $row->{variant_class};

      my $si = $samples{$sid} || die;
      $info{INFO_CANCER_TYPE()} = $si->{diagnosis_short};
      $info{INFO_CANCER_TYPE_FULL()} = $si->{diagnosis_full};
      $info{INFO_CANCER_SUBTYPE()} = $si->{diagnosis_subtype_short};
      $info{INFO_CANCER_SUBTYPE_FULL()} = $si->{diagnosis_subtype_full};
      $info{INFO_CANCER_SUBGROUP()} = $si->{diagnosis_subgroup_short};
      $info{INFO_CANCER_SUBGROUP_FULL()} = $si->{diagnosis_subgroup_full};
      $info{INFO_DATASET()} = $si->{dataset_label};
      $info{INFO_DISEASE_PHASE()} = $si->{sample_type};

#      dump_die($si, "debug");
#      dump_die($row, "Debug");

    } else {
      $chr = $row->{$VPP->f_chr()} || "";
      $pos = $row->{$VPP->f_pos()} || "";
      $ra = $row->{$VPP->f_ra()} || "";
      $va = $row->{$VPP->f_va()} || "";

      if ($annot_mode) {
	foreach my $h (@annot_info_h) {
	  my $label = $annot_info2label{$h};
	  my $v = $row->{$h};
	  $info{$label} = $v;
	}
      }
    }

    my $is_bad;
    $is_bad = 1 unless $chr;
    $is_bad = 1 unless $pos;
    $is_bad = 1 unless $ra or $va;
    # need at least one allele

    my $raw_key = join ".", $chr, $pos, $ra, $va;

    if ($is_bad) {
      if ($SKIP_INCOMPLETE) {
	printf STDERR "skipping malformed row %s\n", $raw_key;
	$malformed++;
	next;
      } else {
	dump_die($row, "allele data incomplete");
      }
    }

    if ($REMOVE_DUPLICATES) {
      # don't reset tracking cache by chromosome, as input file
      # might not be sorted
      if ($saw{$raw_key}) {
	printf STDERR "skipping duplicate %s\n", $raw_key;
	$duplicates++;
	next;
      }
      $saw{$raw_key} = 1;
    }

    my $v = $VPP->get_variant_from_row($row);

    my %r;

#    my $rn = $v->reference_name();
    # cooking can break compatibility with reference FASTA (grch38)
    my $rn = $chr;
    if ($sync_names) {
      $rn = $fai->find_name($rn) || $rn;
    } elsif ($mt_mode) {
      die "mt mode obsolete?";
      if ($rn eq "M") {
	$rn = "MT";
      } elsif ($rn eq "chrM") {
	$rn = "chrMT";
      }
    }

    if ($rn ne $last_rn) {
      &$flush("-force" => 1);
      $last_rn = $rn;
    }

    $r{"#CHROM"} = $rn;

    foreach my $f (qw(ID QUAL)) {
      $r{$f} = VCF_MISSING_VALUE;
    }

    my ($ref, $alt, $pos_vcf);
    if ($v->is_substitution() or $VCF_MODE) {
      # SNV, MNV, or raw VCF alleles
      $pos_vcf = $v->start;
      $ref = $v->reference_allele;
      $alt = $v->variant_allele;
    } elsif ($v->is_deletion()) {
      $pos_vcf = $v->start - 1;
      my $before = $fai->get_chunk(
				   "-id" => $v->reference_name,
				   "-start" => $pos_vcf,
				   "-length" => 1,
				  );
      $ref = $before . $v->reference_allele;
      $alt = $before;
    } elsif ($v->is_insertion()) {
      $pos_vcf = $v->start();
      # already adjusted to base before insertion
      my $before = $fai->get_chunk(
				   "-id" => $v->reference_name,
				   "-start" => $pos_vcf,
				   "-length" => 1,
				  );
      $ref = $before;
      $alt = $before . $v->variant_allele;
    } else {
      my $vtype = $v->get_type();
      if ($vtype eq "complex" and $FLAGS{"complex-ok"}) {
	# VCF 4.2 spec:
	# "...this padding base is not required (although it is permitted)
	# for e.g. complex substitutions or other events where all alleles
	# have at least one base represented in their Strings."
	$pos_vcf = $v->start();
	$ref = $v->reference_allele();
	$alt = $v->variant_allele();
	# so, we can just pass through as-is
      } else {
	printf STDERR "ERROR: unhandled variant type: %s\n", $vtype;
	printf STDERR " - If this variant is a simple indel that uses unparsed VCF-style alleles, you might need to use the -vcf parameter.  vcf2tab.pl is strongly recommended to convert VCF files to tab-delimited.\n";
	printf STDERR " - If this is a complex substitution and that's expected, specify -complex-ok.\n";
	dump_die($row, "ERROR");
      }
    }

    $info{INFO_BAMBINO_KEY()} = $raw_key;

    $r{POS} = $pos_vcf || die;
    $r{REF} = $ref || die;
    $r{ALT} = $alt || die;
    $r{FILTER} = VCF_FILTER_PASS;
    $r{INFO} = join VCF_INFO_DELIMITER, map {$_ . "=" . $info{$_}} sort keys %info;

#    $rpt->end_row(\%r);
    push @out_queue, \%r;
    &$flush("-force" => 0) if @out_queue >= $OUT_QUEUE_FLUSH_SIZE;

  }
#  $rpt->finish();

  &$flush("-force" => 1);
  $wf->finish() if $wf;

  printf STDERR "skipped %d duplicate variants\n", $duplicates if $duplicates;
  printf STDERR "skipped %d malformed variants\n", $malformed if $malformed;

}
