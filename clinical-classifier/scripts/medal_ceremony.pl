#!/usr/bin/env perl
# assign gold/silver/bronze mutation status to variants (SNV/indel/CNV/SV)
# various modes for somatic SNV/indel, germline SNV/indel, CNV, SV
# MNE 9/2013 -
# TO DO:
# - shared/generic $vm importing from a $nsp!
# - standard/universal list of TS/oncogenes

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use FileHandle;
use File::Basename;
use File::Copy;
use Carp qw(confess cluck);
use Sys::Hostname;
use Digest::MD5;
use List::Util qw(max sum min);
use Data::Dumper;
use Bio::SeqIO;

use FileUtils qw(read_simple_file line_delimiter_qc md5_file);
use DelimitedFile;
use Reporter;
use GenomeUtils qw(cook_chromosome_name);
use DBTools qw(selectall_hashref get_dbi_gedi get_dbi_hgmd export_query_to_flatfile);
use AAParser;
use NucleotideSubstitutionParser;
use PostprocessedVariantReport;
use GeneAnnotation;
use VariantMatcher;
use MiscUtils qw(unique_ordered_list unique_list dump_die trim_flanking_whitespace split_list);
use HTMLUtils qw(parse_html_tables);
use NHLBIParser;
use WorkingFile;
use GenomeUtils qw(reverse_complement);
use dbNSFP;
use GeneSymbolStandardizer;
use MapSNVToGenome;
use SJPreferredIsoform;
use NHLBI_tabix;
use GenomicRangeFinder;
use HGNCParser;
use GeneSymbolMapper qw(new_gsm_lite);
use CommandLineRebuilder;
use TabixFile;
use TabixBatchAnnotation qw(TABIX_VARIANT);
use VCFUtils;
use TdtConfig;
use TemporaryFileWrangler;
use RefFlatFile;
use ChrBucketMap;
use VariantOverlap;
use DelimitedFileHP;
use TSOncoDB;
use OMIM_mim2gene;
use SimpleVariantDB qw(
F_VDB_CHR
F_VDB_POS
F_VDB_RA
F_VDB_VA
F_VDB_GENE
F_VDB_AA
F_VDB_CDS
F_VDB_TRANSCRIPT
F_VDB_PMID
);

use constant CLASS_UNKNOWN => "Unknown";
use constant CLASS_GOLD => "Gold";
use constant CLASS_SILVER => "Silver";
use constant CLASS_BRONZE => "Bronze";

use constant CLASS5_BENIGN => "B";
use constant CLASS5_LIKELY_BENIGN => "LB";
use constant CLASS5_UNKNOWN => "U";
use constant CLASS5_LIKELY_PATHOLOGIC => "LP";
use constant CLASS5_PATHOLOGIC => "P";

use constant MAX_CODON_DISTANCE_FOR_CLUSTERED_SEARCH => 5;
#use constant MAX_CODON_DISTANCE_FOR_CLUSTERED_SEARCH => 100;

use constant SV_MAX_BREAKPOINT_GENES_FOR_GOLD => 5;
# JZ email 9/23/2013: "If there are more than 5, make it silver."

my $CNV_MAX_GENES_MEDAL = 10;
# don't call medals (except in HLA) if more than this many genes
#my $CNV_MAX_GENES_ANNOTATE = 10;
#my $CNV_MAX_GENES_ANNOTATE = 100;
my $CNV_MAX_GENES_ANNOTATE = 50000;
# don't annotate if CNV contains more than this many genes
# 10/3/2017: disable annotation limit

use constant CNV_MAX_GENES_FOCAL => 5;
# max genes to consider a focal event
use constant CNV_MIN_LOG2_FOR_HLA => 2;
# minimum log ratio to qualify as highly amplified
use constant CNV_MIN_COPY_NUMBER_GAIN_FOR_GOLD => 10;
# minimum copy number to flag CNV amplifications as gold
use constant CNV_HIGH_LEVEL_ENABLE => 1;
# 7/2014: enable new logic for high-level amplification/deletion
use constant CNV_HIGH_LEVEL_DELETION_COPIES => -1.5;

use constant INDEL_NUCLEOTIDE_NT_FUZZY_MATCH => 3;
# number of bases to tolerate mismatch
# +/- 3 per 10/23/2013 meeting

#use constant GERMLINE_CLINVAR_CLNSIG_TYPES => (4, 5, 6, 7);
#use constant GERMLINE_CLINVAR_CLNSIG_TYPES => (4, 5);
use constant GERMLINE_CLINVAR_CLNSIG_TYPES => (0, 4, 5, 6, 7);
# 12/2013: 0 = uncertain significance

# clinvar: clinical significance codes to use
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Variant Clinical Significance, 0 - unknown, 1 - untested, 2 - non-pathogenic, 3 - probable-non-pathogenic, 4 - probable-pathogenic, 5 - pathogenic, 6 - drug-response, 7 - histocompatibility, 255 - other">
use constant CLINVAR_CLNSIG_U => 0;
use constant CLINVAR_CLNSIG_B => 2;
use constant CLINVAR_CLNSIG_LB => 3;
use constant CLINVAR_CLNSIG_LP => 4;
use constant CLINVAR_CLNSIG_P => 5;

my %CLINVAR_CLNREVSTAT_TO_GOLD_STARS = (
		"no_assertion_criteria_provided" => 0,
		"no_assertion_provided" => 0,
		# 0 gold stars

		"no_interpretation_for_the_single_variant" => 0,
		# not documented (rare)

		"criteria_provided,single_submitter" => 1,
		"criteria_provided,conflicting_interpretations" => 1,
		# 1 gold star

		"criteria_provided,multiple_submitters,no_conflicts" => 2,
                # 2 gold stars

                "reviewed_by_expert_panel" => 3,
		# 3 stars

		"practice_guideline" => 4,
		# 4 stars
				       );
# https://www.ncbi.nlm.nih.gov/clinvar/docs/details/
# known codes appearing in cooked version of CLNREVSTAT as generated
# by vcf2tab.pl and whether these are considered trustable for medal purposes.
# These codes have undergone major changes at least once though;
# maybe better to sanity-check them when that file is built??
my $CLINVAR_MIN_GOLD_STARS_TO_TRUST = 2;

#use constant GERMLINE_MAX_NHLBI => 3;
use constant GERMLINE_MAX_NHLBI => 6;
my $GERMLINE_MAX_NHLBI_FREQ = 0.001;
# maximum allowable records in NHLBI for silver deleterious check
# 12/18/2013: changed to 6 to be consistent with PCGP, or .1% (0.001)
# Gang Wu:
# "Yes, let.s stick to <0.1% (i.e. <7 patients out of from 6500) from now on."
# 7 / 6500 = 0.00107692307692308 (> .001)
my $GERMLINE_ACMG_PM2_MAX_POPULATION_FREQ = 0.0001;
my $GERMLINE_ACMG_BA1_MAX_POPULATION_FREQ = 0.05;

use constant TAG_CLASSIFY => "AutoClass";
use constant ACMG_2015 => "ACMG_2015";

use constant TABIX_BATCH_SIZE_NHLBI => 2000;
use constant TABIX_BATCH_SIZE_DBSNP => 2000;
use constant TABIX_BATCH_SIZE_COSMIC => 2000;
use constant TABIX_BATCH_SIZE_DBNSFP => 500;
use constant TABIX_BATCH_SIZE_THOUSAND_GENOMES => 500;
use constant TABIX_BATCH_SIZE_EXAC => 500;
use constant TABIX_BATCH_SIZE_EXAC_COVERAGE => 2000;
use constant TABIX_BATCH_SIZE_GEDI => 1000;
use constant TABIX_BATCH_SIZE_CLINVAR => 2000;
# how many rows to batch-tabix at once.
# NOT the number of rows tabix itself calls at once, see TabixFile.pm

use constant GERMLINE_ANNOTATION_BOTH_TS_AND_ONCO => 2;
use constant ACMG_PVS1_PENULTIMATE_EXON_WARN_BASES => 50;

my $GERMLINE_USE_MUTATIONASSESSOR_OVER_SIFT = 1;

my $GERMLINE_GSB_LOCKED_TO_COMMITTEE = 1;
my $SOMATIC_GSB_LOCKED_TO_COMMITTEE = 1;

my $IARC_TP53_NM = "NM_000546.5";
# http://p53.iarc.fr/p53Sequences.aspx
# - page doesn't contain .5
# - current version at this time
# - want subversion for collapsed reports
# HACK: should be parsed from HTML!

my $FIELD_GENE = "GeneName";
my $FIELD_REFSEQ = "mRNA_acc";
my $FIELD_AACHANGE = "AAChange";
my $FIELD_DBSNP = "dbSNP";
my $FIELD_NHLBI_FREQ = "NHLBI_frequency";
my $FIELD_PCGP_COUNT_SOMATIC = "PCGP_count_somatic";
my $FIELD_PCGP_DENOMINATOR_SOMATIC = "PCGP_count_somatic_denominator";
my $FIELD_PCGP_COUNT_GERMLINE = "PCGP_count_germline";
my $FIELD_PCGP_DENOMINATOR_GERMLINE = "PCGP_count_germline_denominator";
my $FIELD_FOLDX = "FoldX";
my $FIELD_THOUSAND_GENOMES = "thousand_genomes_frequency";
my $FIELD_OMIM_ID = "OMIM_ID";
my $FIELD_COSMIC = "COSMIC";
my $FIELD_PANEL_DECISION = "paneldecision";
my $FIELD_CHR = "Chr";
my $FIELD_VPOS = "WU_HG19_Pos";
my $FIELD_VPOS_ALT = "WU_HG18_Pos";
my $FIELD_REF_ALLELE = "ReferenceAllele";
my $FIELD_MUT_ALLELE = "MutantAllele";
my $FIELD_CLASS = "Class";
my $FIELD_EXAC_COVERAGE = "__exac_coverage";
my $FIELD_PUBMED_IDS = "PubMed_IDs";

my $FIELD_GEDI_TABIX = "__gedi_tabix";
my $FIELD_COSMIC_TABIX = "__cosmic_tabix";
my $FIELD_COSMIC_TRUNCATING_FLAG = "_cosmic_truncating";
my $FIELD_COSMIC_RECURRENT = "__cosmic_recurrent";
my $FIELD_COSMIC_PUBMED = "__cosmic_pubmed";

my $FIELD_CLINVAR_TABIX = "__clinvar_tabix";
my $FIELD_CLINVAR_TABIX_SITE = "__clinvar_tabix_site";
my $FIELD_TAYLOR_TABIX = "__taylor_tabix";

my $INTERNAL_FIELD_QUEUE_MEDALS = "__medals";
my $INTERNAL_FIELD_QUEUE_REASONS = "__reasons";
my $INTERNAL_FIELD_QUEUE_EVIDENCE = "__evidence";
my $INTERNAL_FIELD_NO_POPULATION_FILTER = "__no_popfilt";

my @COMMITTEE_GL_MEDALS;
# flatfiles representing committee medal decisions.
# hopefully eventually migrated to a database so can be held in one file.

my $CSITH_MEDALS_OVERRIDE_LEGACY_MEDALS = 1;
# csith has the latest medal decisions, which may differ from those
# in old legacy spreadsheets.  If set, load csith variants first
# and ignore old spreadsheet entries for this variants.
# e.g. 4.41747904.C.T: U in csith, LB in an old spreadsheets.


my $ARUP_MEN2_RET_NM = "NM_020630.4";
# HACK: add to parser!

#my $NHGRI_BRCA1_NM = "NM_007300";
# THIS IS A TEMPORARY HACK, the isoform is actually NOT 100% identical!
# JZ 2/19/14
my $NHGRI_BRCA1_NM = "NM_007294.3";
# MNE 2/28/14:
# http://www.ncbi.nlm.nih.gov/nuccore/555931 = U14680.1, protein=AAA73985.1
# https://www.ncbi.nlm.nih.gov/gene/672
# NM_007294.3 -> NP_009225.1
#    - SJ new isoform:
#      - NM_007294 -> NP_009225.1
#      - blastp: 100% match!
#    - SJ old:
#      - NM_007300 -> NP_009231.2
#      - blastp: GAP!

my $NHGRI_BRCA2_NM = "NM_000059.3";
# JZ 2/19/14: 100% identity, safe
# MNE 2/28/14:
# http://www.ncbi.nlm.nih.gov/nuccore/1161383 = U43746.1 -> AAB07223.1
# http://www.ncbi.nlm.nih.gov/gene/675
# NM_000059.3 -> NP_000050.2
#  - blastp: no gaps, 2 mismatches, very close

my @GEDI_RECURRENT_RESTRICT;

# lists of files to process:
my @INPUT_SNV_INDEL_FILES;
my @INPUT_CNV_FILES;
my @INPUT_SV_FILES;
my @INPUT_GL_FILES;

my $DEBUG_PCGP_INDEL = 0;
my $REVIEWABLE_FLANKING_BUFFER_NT = 1000;
my $TRIM_COSMIC_INFO = 1;
my $TRIM_HGMD_INFO = 1;
my $TRIM_TP53_INFO = 1;

my $DOWNGRADE_SOMATIC_MEDAL_IF_HIGH_NHLBI = 1;
# JZ 3/8/2016: downgrade somatic SNV/indel medals if high NHLBI frequency

my $SOMATIC_GENE_GSM = 1;
my $SOMATIC_TS_ONCO_DB = 1;

my $GERMLINE_GENE_GSM = 1;
my $GERMLINE_TS_ONCO_DB = 1;

my $SV_GSM = 1;
my $CNV_GSM = 1;

my $ENABLE_SJPI_GSM = 1;
# resolve gene symbol ambiguity w/SJ preferred isoform list

my @NSFP_FIELDS = qw(
		      Interpro_domain
		      UniSNP_ids
		      SIFT_pred
		      SIFT_score
		      SIFT_score_converted
		      Polyphen2_HDIV_pred
		      Polyphen2_HDIV_score
		      Polyphen2_HVAR_pred
		      Polyphen2_HVAR_score
		      LRT_Omega
		      LRT_pred
		      LRT_score
		      LRT_score_converted
		      MutationAssessor_pred
		      MutationAssessor_score
		      MutationAssessor_score_converted
		      MutationTaster_pred
		      MutationTaster_score
		      MutationTaster_score_converted
		      FATHMM_pred
		      FATHMM_score
		      FATHMM_score_converted
		      CADD_raw
		      CADD_phred
		      REVEL_score
		   );
# passthrough



my %MEDAL_SORT_WEIGHT = (
			 CLASS_GOLD() => 300,
			 CLASS_SILVER() => 200,
			 CLASS_BRONZE() => 100,
			 CLASS_UNKNOWN() => 0.01,  # for hash presence check

			 CLASS5_BENIGN() => 50,
			 CLASS5_LIKELY_BENIGN() => 100,
			 CLASS5_UNKNOWN() => 0.01,
			 CLASS5_LIKELY_PATHOLOGIC() => 200,
			 CLASS5_PATHOLOGIC() => 300,

			);

my %FLAGS;
load_config_file();
# pre-load @ARGV

my $FUSION_BUILDER_COLUMN = "Fusion Gene";
my $INDEL_MATCH_BASIC_TYPE = 1;
#my $INDEL_MATCH_BASIC_TYPE = 0;
# if set, basic type (insertion/deletion) must match reference database.
# counterargument: indels can both cause frameshifts at same position, so ???
my $INDEL_MATCH_SIZE = 1;

my $REFERENCE_SANITY_CHECK = 0;
# sanity-check reference base when loading variant databases

my $GERMLINE_MIN_PLATFORMS_FOR_DAMAGING = 2;

my $COSMIC_HOTSPOT_MIN_CONFIRMED_PATIENTS = 10;


my $TABIX_INDEL_EQUIVALENCE_DISTANCE = 25;

my $ENABLE_TABIX_COSMIC = 1;
my $ENABLE_TAYLOR_HOTSPOTS = 0;
my $ENABLE_SOMATIC_EXAC = 1;

my $ENABLE_KNOWN_PROMOTER_SITES = 0;
my $ENABLE_KNOWN_PROMOTER_REGIONS = 0;

my $VERBOSE_COSMIC = 0;
my $VERBOSE_SOMATIC_PROGRESS = 0;
my $VERBOSE_GERMLINE_PROGRESS = 0;
my $VERBOSE_AA_PARSE_WARNING = 0;
my $VERBOSE_SILENT_AND_SUPPRESSED = 0;

my $COMMITTEE_CONFLICT_WARN = 1;
# if true, warn rather than exiting with an error if
# inconsistencies detected in committee medal decisions.

my $COSMIC_MIN_VALIDATED_SAMPLES_FOR_RECURRENT = 3;
# COSMIC variants must have at least this many validated
# samples to be considered recurrent

my $CNV_GERMLINE_MEDAL_ANY_SIZE = 0;
my $CNV_GERMLINE_MIN_COPY_DELTA;

my @EXAC_POPULATION_TAGS = qw(
			  AFR
			  AMR
			  EAS
			  FIN
			  NFE
			  SAS
			  OTH
		       );
my @EXAC_HEADERS;

my @EXAC_COVERAGE_H_TABIX = qw(
			 mean
			 median
			 1
			 5
			 10
			 15
			 20
			 25
			 30
			 50
			 100
			     );
my @EXAC_COVERAGE_H = map {"ExAC_cvg_" . $_} @EXAC_COVERAGE_H_TABIX;

my @COMMAND_OPTIONS = (
	   "-config=s",
	   # alternative for specifying many of the below parameters,
	   # e.g.
	   # /nfs_exports/genomes/1/projects/ClinicalSeq/medal_config/2013_09_20/medal_config_2013_09_20.txt
           # DEPRECATED!

           "-genome=s",
           # new/preferred method: use genome config

	   "-gold=s",
	   "-silver=s",
	   # gold and silver gene lists, e.g.
	   # /nfs_exports/genomes/1/projects/ClinicalSeq/Gene4Report/GoldGene.lst
	   # /nfs_exports/genomes/1/projects/ClinicalSeq/Gene4Report/SilverGene.lst

	   "-no-pvr",
	   "-no-nhlbi",
		       "-no-cosmic",

	   "-in=s",
	   # single variant report to process (type autodetected)
	   "-in-list=s",
	   # list of variant reports to process (type autodetect)
	   # "-in-snv-indel=s" => \@INPUT_SNV_INDEL_FILES,
	   # "-in-cnv=s" => \@INPUT_CNV_FILES,
	   # "-in-sv=s" => \@INPUT_SV_FILES,
	   "-in-snv-indel=s",
	   "-in-cnv=s",
	   "-single-cnv=s" => \@INPUT_CNV_FILES,
	   "-in-sv=s",
	   "-in-gl=s",

	   "-single-gl=s" => \@INPUT_GL_FILES,
	   "-single-snv-indel=s" => \@INPUT_SNV_INDEL_FILES,
	   "-single-sv=s" => \@INPUT_SV_FILES,

#	   "-out=s",
	   # outfile name (optional)

	   "-outfile-fq",
	   # if specified outfile is written to same dir as infile
	   # (default: cwd)

	   "-genes-manual=s",
	   # manual gene classification (gold/silver, TS/oncogene)
	   # e.g. SMC3, which was not classified as either TS/oncogene
	   # and was field promoted by Matt Parker

	   "-cosmic-cleaned=s",
	   "-cosmic=s",
	   "-cosmic-gold=s",
	   "-cosmic-silver=s",

		       "-cosmic-pubmed-summary=s",
		       # new cosmic code for hg38+:
		       # replaces -cosmic, -cosmic-gold, -cosmic-silver

	   "-gedi-recurrent=s",
	   # GeDI dump of recurrent sites, e.g.
	   # /nfs_exports/genomes/1/projects/ClinicalSeq/PCGPRecurrentMutation/GeDI_mutation_20130425.txt
	   "-gedi-recurrent-restrict=s" => \@GEDI_RECURRENT_RESTRICT,

	   "-cnv=s",
	   # CNV configuration data, e.g.
	   # /nfs_exports/genomes/1/projects/ClinicalSeq/FusionGenes/CNVCheck.txt
	   "-cnv-annotations=s",
	   # supplementary annotations, e.g.
	   # /nfs_exports/genomes/1/projects/ClinicalSeq/FusionGenes/CancerGeneCensus/Deletion.txt

	   "-gene-exon-region-dir=s",
	   # in CNV processing, used to find genes associated with intervals
	   # (we need edge genes and they are alpha-sorted rather than genomic)
	   # e.g.
	   # /nfs_exports/genomes/1/projects/ClinicalSeq/medal_config/2013_09_20/medal_config_2013_09_20.txt

	   "-sv=s",
	   # SV configuration data, e.g.
	   # /nfs_exports/genomes/1/projects/ClinicalSeq/FusionGenes/SVCheck.txt
	   "-sv-manual=s",
	   # manually-curated/promoted entries, e.g.
	   # /nfs_exports/genomes/1/projects/ClinicalSeq/FusionGenes/SVCheck_manual.txt
	   # keep these in a separate file: when -sv file is refreshed,
	   # any manual changes made to it will be lost.

	   "-dump-gedi",
	   # generate GeDI SNV/indel dump files
	   "-dump-hgmd",
	   # export HGMD flatfile from GeDI
		       "-hgmd-table=s",

	   "-dump-asu-tert",
	   # parse ASU TERT HTML and export
	   "-generate-nhlbi-blacklist=s",
	   # e.g. /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/EVS_SNPs/ESP6500/SNP_txt/
	   "-nhlbi-x",
	   "-dump-clinvar",
	   # export ClinVar query from GeDI

	   "-verbose-nhlbi",
	   "-verbose-cosmic",

	   "-cnv-test=i",

	   "-fb-column=s" => \$FUSION_BUILDER_COLUMN,
	   # manual override of fusion column name

	   #
	   # GERMLINE classification options:
	   #
	   "-iarc-tp53-germline=s",
	   "-iarc-tp53-somatic=s",
	   # these parse the raw database files, however
	   # this is now deprecated due to manual curation
	   # (e.g. P72R).
#	   "-iarc-tp53-curated=s",
	   # manually curated file, currently just AA changes

#	   "-hgmd-dm=s",
	   # old format, obsolete
#	   "-hgmd-dm-mp=s",
	   # SNVs only, obsolete
#	   "-hgmd-dm-mp-indels=s",
	   # 1/2014, obsolete
#	   "-hgmd-dm-mp-indels-and-splices=s",
	   "-hgmd-dm-mp-all=s",
	   # 3/2014: everything
		       "-no-hgmd",

	   "-asu-tert=s",
#	   "-clinvar=s",
	       # direct VCF parsing: obsolete

	   "-gl-test=i",
	   "-gl-gold-cancer-ranges=s",
	   # raw file

	   "-gl-grf-debug",

#	   "-gl-gold-non-cancer=s",
#	   "-gl-silver=s",
#	   "-gl-bronze=s",
#	   "-gl-pcgp-gold-db=s",
#	   "-gl-umd-dir=s",
#	   "-germline-report-gold-only=i" => \$GERMLINE_REPORT_GOLD_ONLY,

	   "-gl-arup=s",
	   "-gl-umd-flatfile=s",
	   # end germline classification options

	   "-fasta-dir=s",
	   # directory of chromosomes
	   # (used to sanity-check reference bases in SNVs)

	   "-rsc=i" => \$REFERENCE_SANITY_CHECK,

		       "-hack-1kg",
	   "-hack",
	   "-hack2",
	   "-hack3",
	   "-hack4",
	   "-hack5",
	   "-hack6",
	   "-hack-nhlbi-tabix",
	   "-hack-nsfp",
		       "-hack-refgene-sv-symbols",
	   "-hack8",
	   "-hack9",
	   "-hack10",
	   "-hack12",
	   "-hack13",
	   "-hack14",
	   "-hack15",
	   "-hack16",
	   "-hack17",
		       "-debug-genes-to-ger=s",
		       "-hack18",
	   "-hack-gsm",
	   "-hack-gsm-bare",

		       "-hack-committee=s",

           "-hack-cosmic-orig",
           "-hack-cosmic-tabix",
           "-hack-cosmic-tabix-gl",
           "-hack-nsfp-tabix",
           "-hack-clinvar-tabix",
		       "-hack-taylor",

		       "-hack-mc-db",

	   "-aa-verbose",

	   "-one-list-to-rule-them-all=s",
	   # argument = committee-approved list of germline genes
	   # - ACMG genes + committee additions
	   # e.g.
	   # /nfs_exports/genomes/1/projects/ClinicalSeq/medal_config/2013_09_20/gold_genes_committee_2014_01_02.txt
	   "-debug-one",
	   "-debug-no-rsc",

	   "-rb1-flatfile=s",
	   # /nfs_exports/genomes/1/projects/ClinicalSeq/germline/RB1/scrape_2013_12_20/export_RB1.tab.map.tab
#	   "-clinvar-gedi-flatfile=s",
#	   "-asu-tert-flatfile=s",
	   # /nfs_exports/genomes/1/projects/ClinicalSeq/germline/ASU_TERT/export_ASU_TERT.tab
	   "-collapse-meta-list=s",
	   # remove duplicates, etc. from

	   "-germline-gene-survey=s",
	   "-clinvar-manual=s",
	   "-gene-info=s",

	   "-parse-rb1=s",
	   # /nfs_exports/genomes/1/projects/ClinicalSeq/germline/RB1/scrape_2013_12_20
	   "-rb1-type=s",

	   "-map-rb1=s",
	   "-rb1-fasta=s",
	   "-rb1-hack=s",

	   "-uniprot-idmapping=s",
	   # ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz

	   "-preferred-isoforms=s",
	   # /nfs_exports/apps/gnu-apps/NextGen/nextgensupport2/NCBI/RefSeq/gene_transcript_matrix.withNM.mod
#	   "-mini-gl-prep=s",
	   "-meta-reannotate=s",
	   # reannotate a meta gene report to SJ accessions

	   "-mini-gl-run=s",
	   # annotate a SJ variants file with annotations
	   # from a meta report (typically reannotated with -mini-gl-prep
	   # - requires gl file spec options

	   "-gl-bad-snv=s",
	   # manual blacklist of bad variants, e.g. TP53 P72R
	   # will probably be eventually replaced with committee review
	   # feedback system TBD

	   "-generate-cosmic-hotspot-list",
	   "-cosmic-hotspots=s",

#	   "-gl-committee-medals=s",
	   # will be replaced by -gl-committee-medals2
#	   "-gl-committee-medals2=s" => \@COMMITTEE_GL_MEDALS,
#	   "-gl-current-committee-medals=s",
	   "-gl-current-committee-medals=s" => \@COMMITTEE_GL_MEDALS,
		       "-committee-debug",

		       "-tabix-dbsnp=s",
		       "-tabix-nhlbi=s",
		       "-tabix-cosmic=s",
		       "-tabix-gedi=s",
		       "-tabix-clinvar=s",

		       "-tabix-dbnsfp=s",
		       "-tabix-dbnsfp-pos=s",
		       # override column name for position field

	   "-gl-pcgp-population=s",
	   "-gl-apc-flatfile=s",
	   "-gl-msh2-flatfile=s",

	   "-debug-start",

	   "-debug-db-hit",

	   "-seatbelts-are-for-sissies",

	   "-nhgri-brca1=s",
	   "-nhgri-brca2=s",

	   "-no-sort",

	   "-hack-reannotate=s",

	   "-show-isoforms-used=s",
	   # input: meta report

	   "-committee-conflict-warn",
	   # if a conflict detected in committee calls,
	   # warn instead of dying

	   "-indel-match-size=i" => \$INDEL_MATCH_SIZE,

	   "-sv-gene-check",
	   "-hgnc=s",
   # ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc_complete_set.txt.gz
	   # /nfs_exports/genomes/1/Homo_sapiens/HGNC/hgnc_complete_set.txt

	   "-cnv-gene-patch",
	   "-cnv-sv-gene-patch",

	   "-dump-sv-genes",
	   "-dump-cnv-genes",
	   # export list of gene symbols used in CNV analysis,
	   # for consistency check with HUGO etc.

	   "-verify-cnv-genes",
	   # verify that CNV genes are findable via gene lookup
	   # interface used in CNV classification

	   "-map-snv-gold-genes-to-sv",

	   "-refflat-fb=s",
	   # refFlat file used by FusionBuilder ("sharp" version)

#	   "-gold-genes-mapped-to-fb=s",

	   "-erin-genes-cnv",
	   # report lists of genes used in CNV analysis
	   "-erin-genes-sv",
	   # report lists of genes used in SV analysis

	   "-sv-silver-genes-fb=s",
	   # SV analysis: gene list to assign silver.
	   # gene symbols must be available in the set used
	   # by FusionBuilder (i.e. "sharp" refFlat).


	   "-generate-germline-cnv-config",
	   # generate germline CNV analysis config file from
	   # (a) committee-approved gene list and
	   # (b) tumor suppressor status lists
	   "-gl-reportable-genes=s",
	   # committee-approved reportable germline gene list

#	   "-gl-reportable-genes-ger=s",
#	   "-gl-reviewable-genes-ger=s",
           # ...mapped to GENE_EXON_REGION

#	   "-gl-reportable-genes-refflat=s",
#	   "-gl-reviewable-genes-refflat=s",
           # ...mapped to refFlat (SV / FusionBuilder)

#	   "-gl-cnv=s",
	   # germline CNV analyis config file

	   "-cnv-germline",
	   # flag indicating CNV analyis should be run in germline mode

	   "-map-gl-extended-to-ger",
	   # translate/verify gene symbols in "gene ranges" matches
	   # GENE_EXON_REGION symbols

	   "-sanity-somatic-cosmic-genes-in-main-lists",

	   "-meta-classifier-gene-list=s",

	   "-sheila-sv-to-manual=s",
	   # append sheila's SV additions to manual SV list

	   "-meta-qc",

	   "-dump-paneldecision",
	   "-dump-germline-reportable-genes",

	   "-gedi-user=s",
	   "-gedi-password=s",
	   "-gedi-dev",

	   "-hgmd-parse-test",
	   "-committee-parse-test",
		       "-committee-export",

		       "-config-checksum=i",
		       # 0 = don't report basename of config value in output.
		       #     What we want most of the time, will make for
		       #     simpler diffs between research and clinical.
		       # 1 = include basename (diagnostic if differences)
		       # 2 = include fully-qualified
		       "-config-dump",

		       "-export-for-clinical=s",
		       # export data for a single param (may be a list)

		       "-exac=s",
		       # ExAC database, TABIX-indexed .vcf
		       "-exac-vcf2tab=s",
		       # cooked and tabix'd version
		       "-exac-coverage=s",
		       # coverage (merged tabix)

		       "-hack-exac",
		       "-hack-exac-vcf2tab",
		       "-hack-exac-coverage",

		       "-check-gold-isoforms=s",

		       "-private-variants=s",
		       "-gl-custom-genes=s",
		       "-gl-custom-genes-lite=s",
		       "-gl-custom-allow-interval-fail",
		       "-gl-no-reviewable",

		       "-hack-exac2",

		       "-tabix-thousand-genomes=s",


		       "-known-promoter-intervals=s",
		       "-known-promoter-sites=s",

		       "-enable-tabix-cosmic=i" => \$ENABLE_TABIX_COSMIC,
		       "-enable-taylor-hotspots=i" => \$ENABLE_TAYLOR_HOTSPOTS,

		       "-enable-known-promoter-regions=i" => \$ENABLE_KNOWN_PROMOTER_REGIONS,

		       "-taylor-hotspots=s",

		       "-2bit=s",
		       # currently not used as config system is queried
		       # directly, however add so will be part of
		       # config checksum

		       "-debug-ram",

		       "-cnv-max-genes-annotate=i" => \$CNV_MAX_GENES_ANNOTATE,

		       "-cnv-germline-any-size=i" => \$CNV_GERMLINE_MEDAL_ANY_SIZE,
		       "-cnv-germline-min-copy-delta=f" => \$CNV_GERMLINE_MIN_COPY_DELTA,

		       "-config-somatic-cnv-add-manual-genes",


		       "-max-population-frequency=f" => \$GERMLINE_MAX_NHLBI_FREQ,
		       "-count-medalable-nt=s",

		       "-hack-sjpi",

		       "-dump-modules",

		       "-enable-sv-gsm=i" => \$SV_GSM,
		       "-enable-cnv-gsm=i" => \$CNV_GSM,

		       "-clinvar-min-gold-stars=i" => \$CLINVAR_MIN_GOLD_STARS_TO_TRUST,
		       "-sqlite=s",

		       "-user-variants=s",
		       "-user-blacklist=s",
		       "-hack-user-variants",
		      );

GetOptions(\%FLAGS, @COMMAND_OPTIONS);
load_genome_config();

my $CLASS_TRUNCATION = 1;
my $CLASS_SPLICE = 2;
my $CLASS_SPLICE_INTRON = 3;

my %variant_class_truncating = (
				"silent" => 0,
				"missense" => 0,
				"utr_5" => 0,
				"utr_3" => 0,
				"exon" => 0,
				"intron" => 0,

				"proteindel" => 0,
				"proteinins" => 0,
				# in-frame

				"nonsense" => $CLASS_TRUNCATION,
				"frameshift" => $CLASS_TRUNCATION,

				"splice" => $CLASS_SPLICE,

				"splice_region" => $CLASS_SPLICE_INTRON,
				#				"splice_region" => 0,
				# Gang:
				# spce_region is an intronic mutation not a truncating mutation.
				# JZ: SILVER, not gold

				"promoter_region" => 0,
			       );
# store everything lowercase:
# - standard postprocessed reports use lowercase
# - SITH reports use uppercase  :/

my %variant_class_inframe_indel = (
				   "proteinins" => 1,
				   "proteindel" => 1,
				  );

if ($FLAGS{"erin-genes-sv"}) {
  erin_genes_sv();
  exit(0);
} elsif ($FLAGS{"dump-hgmd"}) {
  dump_hgmd();
  exit(0);
} elsif (my $fn = $FLAGS{"one-list-to-rule-them-all"}) {
  dump_raw_lists($fn);
  exit(0);
} elsif (my $fn3 = $FLAGS{"collapse-meta-list"}) {
  collapse_meta_list($fn3);
  exit(0);
} elsif (my $fn2 = $FLAGS{"germline-gene-survey"}) {
  germline_gene_survey($fn2);
  exit(0);
} elsif ($FLAGS{"dump-clinvar"}) {
  dump_clinvar();
  exit(0);
} elsif ($FLAGS{"dump-asu-tert"}) {
  my $rows = load_asu_tert();
  my $outfile = "export_ASU_TERT.tab";

  my @fields = grep {not(/References/)} sort keys %{$rows->[0]};
  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => \@fields
			);

  foreach my $row (@{$rows}) {
    $rpt->end_row($row);
  }
  $rpt->finish();
  exit(0);
} elsif ($FLAGS{hack}) {
  $FLAGS{"asu-tert"} = "/nfs_exports/genomes/1/projects/ClinicalSeq/germline/ASU_TERT/diseases_excerpt_TERT.html";
  my $rows = load_asu_tert();
  #  my $rows = load_asu_tert_new();

  foreach my $row (@{$rows}) {
    foreach (sort keys %{$row}) {
      printf "%s: %s\n", $_, $row->{$_};
    }
  }
  exit(0);
} elsif ($FLAGS{hack2}) {
  condense_output("germline_variants_list.tab");
  exit(0);
} elsif ($FLAGS{hack3}) {
  sanity_check_ordering("export_RB1.tab.map.tab");
  exit(0);
} elsif ($FLAGS{hack4}) {
  $FLAGS{"debug-ram"} = 1;
  ram_debug("HGMD start");
  my $vm = get_vm_hgmd();
  ram_debug("HGMD end");
  exit(0);
} elsif ($FLAGS{hack6}) {
  tp53_hack();
  exit(0);
} elsif ($FLAGS{"hack-sjpi"}) {
  hack_sjpi();
  exit(0);
} elsif ($FLAGS{"dump-modules"}) {
  print Dumper \%INC;
  exit(0);
} elsif ($FLAGS{"hack-gsm"}) {
  hack_gsm();
  exit(0);
} elsif ($FLAGS{"hack-gsm-bare"}) {
  hack_gsm_bare();
  exit(0);
} elsif ($FLAGS{"hack-refgene-sv-symbols"}) {
  hack_refgene_sv_symbols();
  exit(0);
} elsif ($FLAGS{"hack-nsfp"}) {
  init_infiles();
  hack_nsfp_tabix();
  exit(0);
} elsif ($FLAGS{"hack-cosmic-orig"}) {
  init_infiles();
  hack_cosmic_orig();
  exit(0);
} elsif ($FLAGS{"hack-cosmic-tabix"}) {
  init_infiles();
  hack_cosmic_tabix();
  exit(0);
} elsif ($FLAGS{"hack-cosmic-tabix-gl"}) {
  init_infiles();
  hack_cosmic_tabix_germline();
  exit(0);
} elsif ($FLAGS{"hack-nsfp-tabix"}) {
  init_infiles();
  hack_nsfp_tabix();
  exit(0);
} elsif ($FLAGS{"hack-clinvar-tabix"}) {
  init_infiles();
  hack_clinvar_tabix();
  exit(0);
} elsif ($FLAGS{"hack-clinvar2"}) {
  init_infiles();
  hack_clinvar2();
  exit(0);
} elsif ($FLAGS{"hack-taylor"}) {
  init_infiles();
  hack_taylor();
  exit(0);
} elsif ($FLAGS{"hack-mc-db"}) {
  init_infiles();
  my $vm = get_vm_mc_db("-name" => "ALSoD");

  die $vm->find_aa_specific_codon("-gene" => "FUS", "-aa" => "K510N");

  exit(0);
} elsif ($FLAGS{"hack-1kg"}) {
  init_infiles();
  thousand_genomes_hack();
  exit(0);
} elsif ($FLAGS{"hack-exac-vcf2tab"}) {
  init_infiles();
  ram_debug("ExAC start");
  hack_exac_vcf2tab();
  ram_debug("end");
  exit(0);
} elsif ($FLAGS{"hack-exac-coverage"}) {
  init_infiles();
  hack_exac_coverage();
  exit(0);
} elsif ($FLAGS{"hack-nhlbi-tabix"}) {
  init_infiles();
  nhlbi_hack_tabix();
  exit(0);
} elsif ($FLAGS{"sv-gene-check"}) {
  sv_gene_check();
  exit(0);
} elsif ($FLAGS{"cnv-gene-patch"}) {
  cnv_gene_patch();
  exit(0);
} elsif ($FLAGS{"cnv-sv-gene-patch"}) {
  cnv_sv_gene_patch();
  exit(0);
} elsif ($FLAGS{hack10}) {
  init_infiles();
  grf_hack();
  exit(0);
} elsif ($FLAGS{hack12}) {
  umd_hack();
  exit(0);
} elsif ($FLAGS{hack13}) {
  pcgp_somatic_hack();
  exit(0);
} elsif ($FLAGS{hack14}) {
  column_qc_hack();
  exit(0);
} elsif ($FLAGS{hack17}) {
  pcgp_indel_hack();
  exit(0);
} elsif ($FLAGS{"debug-genes-to-ger"}) {
  debug_genes2ger($FLAGS{"debug-genes-to-ger"});
  exit(0);
} elsif ($FLAGS{hack18}) {
  hack18();
  exit(0);
} elsif ($FLAGS{"hack-user-variants"}) {
  init_infiles();
  hack_user_variants();
  exit(0);
} elsif (my $mf = $FLAGS{"show-isoforms-used"}) {
  show_isoforms_used($mf);
  exit(0);
} elsif (my $hf = $FLAGS{"hack-reannotate"}) {
  hack_reannotate($hf);
  exit(0);
} elsif ($FLAGS{hack5}) {
  base_match_test();
  exit(0);
} elsif (my $r1d = $FLAGS{"parse-rb1"}) {
  parse_rb1($r1d);
  exit(0);
} elsif (my $ff = $FLAGS{"map-rb1"}) {
  map_rb1($ff);
  exit(0);
} elsif ($FLAGS{"meta-reannotate"}) {
  meta_reannotate();
  exit(0);
} elsif ($FLAGS{"mini-gl-run"}) {
  run_mini_germline();
  exit(0);
} elsif ($FLAGS{"generate-cosmic-hotspot-list"}) {
  generate_cosmic_hotspot_list();
  exit(0);
} elsif ($FLAGS{"dump-cnv-genes"}) {
  dump_cnv_genes();
  exit(0);
} elsif ($FLAGS{"dump-sv-genes"}) {
  dump_sv_genes();
  exit(0);
} elsif ($FLAGS{"verify-cnv-genes"}) {
  verify_cnv_genes();
  exit(0);
} elsif ($FLAGS{"map-snv-gold-genes-to-sv"}) {
  map_snv_gold_genes_to_sv();
  exit(0);
} elsif ($FLAGS{"generate-germline-cnv-config"}) {
  generate_germline_cnv_config();
  exit(0);
} elsif ($FLAGS{"map-gl-extended-to-ger"}) {
  map_gl_extended_to_ger();
  exit(0);
} elsif ($FLAGS{"sanity-somatic-cosmic-genes-in-main-lists"}) {
  sanity_check_cosmic_somatic_genes_in_main_lists();
  exit(0);
} elsif ($FLAGS{"meta-classifier-gene-list"}) {
  meta_classifier_gene_list();
  exit(0);
} elsif (my $shl = $FLAGS{"sheila-sv-to-manual"}) {
  sheila_sv_to_manual($shl);
  exit(0);
} elsif ($FLAGS{"meta-qc"}) {
  meta_qc_checks();
  exit(0);
} elsif ($FLAGS{"dump-paneldecision"}) {
  export_panel_decisions();
  exit(0);
} elsif ($FLAGS{"dump-germline-reportable-genes"}) {
  export_germline_reportable_genes();
  exit(0);
} elsif ($FLAGS{"hgmd-parse-test"}) {
  my $vm = get_vm_hgmd();
  printf STDERR "parsing finished\n";
  exit(0);
} elsif ($FLAGS{"committee-parse-test"}) {
  committee_hack();
  exit(0);
} elsif ($FLAGS{"committee-export"}) {
  committee_export();
  exit(0);
} elsif (exists $FLAGS{"config-checksum"}) {
  load_genome_config("-checksum" => $FLAGS{"config-checksum"});
  exit(0);
} elsif (exists $FLAGS{"config-dump"}) {
  load_genome_config();
  exit(0);
} elsif ($FLAGS{"export-for-clinical"}) {
  export_for_clinical();
  exit(0);
} elsif (exists $FLAGS{"hack-exac"}) {
  # exac testing
  hack_exac();
  exit(0);
} elsif (exists $FLAGS{"hack-exac2"}) {
  # exac-only annotation
  hack_exac2();
  exit(0);
} elsif (exists $FLAGS{"check-gold-isoforms"}) {
  check_gold_isoforms();
  exit(0);
} elsif ($FLAGS{"config-somatic-cnv-add-manual-genes"}) {
  config_somatic_cnv_add_manual_genes();
  exit(0);
} elsif ($FLAGS{"count-medalable-nt"}) {
  count_medalable_nt();
  exit(0);
}

init_infiles();

foreach my $ref (
		 \@INPUT_CNV_FILES,
		 \@INPUT_SV_FILES,
		 \@INPUT_GL_FILES,
		 \@INPUT_SNV_INDEL_FILES
		) {
  if (@{$ref}) {
    outfile_safety_check($ref);
    file_sanity_check($ref);
  }
}

if (@INPUT_CNV_FILES) {
  my %options;
  if ($FLAGS{"cnv-germline"}) {
    %options = ("-germline" => 1);
  } else {
    %options = ("-somatic" => 1);
  }
  if ($CNV_GSM) {
    printf STDERR "running GSM CNV classifier\n";
    run_cnvs_integrated_new(%options);
  } else {
    run_cnvs_integrated_old(%options);
  }
}

if (@INPUT_SV_FILES) {
  if ($SV_GSM) {
    printf STDERR "running GSM SV classifier\n";
    run_svs_new();
  } else {
    run_svs_old();
  }
}

run_germline() if @INPUT_GL_FILES;
run_somatic_snv_indel() if @INPUT_SNV_INDEL_FILES;

sub run_somatic_snv_indel {
  #
  #  primary gold/silver gene lists:
  #
  my $gold_genes = flatfile_to_hash($FLAGS{"gold"} || die "-gold");
  my $silver_genes = flatfile_to_hash($FLAGS{"silver"} || die "-silver");
  my $omim = new OMIM_mim2gene();

  my $ts_onco_db;
  if ($SOMATIC_TS_ONCO_DB) {
    $ts_onco_db = new TSOncoDB(
			       "-gsm" => new_gsm_lite()
			      );
  }

  my $gsm_gold_genes = new_gsm_lite();
  foreach my $g (keys %{$gold_genes}) {
    $gsm_gold_genes->add_gene("-gene" => $g);
  }
  my $gsm_silver_genes = new_gsm_lite();
  foreach my $g (keys %{$silver_genes}) {
    $gsm_silver_genes->add_gene("-gene" => $g);
  }

  my $gene2class;
  # retire in future release
  if ($SOMATIC_TS_ONCO_DB) {
    $gene2class = {};
  } else {
    $gene2class = load_gene_class_files();
  }

  #
  # additional genes later promoted to gold/silver:
  #
  my $f_manual = $FLAGS{"genes-manual"} || die "-genes-manual";
  my $df = new DelimitedFile(
			     "-file" => $f_manual,
			     "-headers" => 1
			    );
  while (my $row = $df->get_hash()) {
    my $gene = $row->{Gene} || die;
    my $type = lc($row->{Type} || die);
    if ($type eq "gold") {
      $gold_genes->{$gene} = 1;
      $gsm_gold_genes->add_gene("-gene" => $gene);
    } elsif ($type eq "silver") {
      $silver_genes->{$gene} = 1;
      $gsm_silver_genes->add_gene("-gene" => $gene);
    } else {
      die;
    }
    $gene2class->{$gene} = $row->{Class} || die;
  }

  my $germline_reportable = get_germline_reportable_genes();
  # no need for context-corrected version as using symbol mapper
  foreach my $g (@{$germline_reportable}) {
    $gsm_gold_genes->add_gene("-gene" => $g);
  }

  ram_debug("somatic_snv/start");

  #
  #  load COSMIC sites:
  #
  my $vm_cosmic_somatic;
  my $cosmic_indels;
  my $cosmic;
  my $new_cosmic = $FLAGS{"cosmic-pubmed-summary"};

  if ($FLAGS{"no-cosmic"}) {
    printf STDERR "DEBUG: COSMIC disabled\n";
    $vm_cosmic_somatic = new VariantMatcher();
    $cosmic_indels = {};
    $cosmic = {};
  } else {

    if ($new_cosmic) {
      # required info will be in batch tabix results
      $cosmic = $cosmic_indels = "do_not_use";
    } else {
      $cosmic_indels = load_cosmic_indels();
      ram_debug("somatic_snv/after cosmic indels");
      $cosmic = load_cosmic_data();
      ram_debug("somatic_snv/after cosmic data");
      # - early parsers (summer 2013)
      # - replace with VariantMatcher model?
    }

    unless ($FLAGS{"tabix-cosmic"}) {
      # old style
      $vm_cosmic_somatic = parse_cosmic(
					"-file" => ($FLAGS{"cosmic-cleaned"} || die "-cosmic-cleaned"),
					"-all-genes" => 1,
					"-all-variants" => 1
				       );
      # 2/2014: global reporting
    }
  }
  ram_debug("somatic_snv/after all cosmic");

  my $GL_REPORTABLE_GENES_GER;
  # TO DO: prune out in next release
  if ($SOMATIC_GENE_GSM) {
    $GL_REPORTABLE_GENES_GER = "do_not_use";
    # won't be accessed
  } else {
    $GL_REPORTABLE_GENES_GER = read_simple_file(($FLAGS{"gl-reportable-genes-ger"} || die "-gl-reportable-genes-ger"), "-hash1" => 1);
    # this list scheduled for demolition
  }

  my $GL_TS_INFO;
  if ($SOMATIC_TS_ONCO_DB) {
    $GL_TS_INFO = {};
  } else {
    # retire in future
    $GL_TS_INFO = get_gl_ts_list();
  }

  add_custom_ts(
		"-ts-genes" => $GL_TS_INFO,
		"-ts-onco-db" => $ts_onco_db
	       );
  # TO DO: remove $GL_TS_INFO in favor of ts_onco_db

  #
  #  process SNVs/indels:
  #
  my $somatic_start_time = time;

  ram_debug("somatic_snv/before file loop");

  foreach my $infile_raw (@INPUT_SNV_INDEL_FILES) {
    my ($infile, $outfile) = get_infile_and_outfile($infile_raw);
    printf STDERR "processing %s...\n", $infile;
    ram_debug("somatic_snv/infile $infile");

    line_delimiter_qc("-file" => $infile, "-delimiter" => "\t");

    if ((-s $infile || 0) <= 1) {
      # some files contain only a newline
      printf STDERR "WARNING: skipping empty file %s\n", $infile;
      next;
    }

    my @headers;
    my $rows;

    if ($FLAGS{"no-pvr"}) {
      # - not a standard postprocessed variant report
      # - standard header labels are still required
      # - example: Gang's PCGP_germline_indels_in35reportable_genes.txt
      $rows = [];
      my $df = new DelimitedFile(
				 "-file" => $infile,
				 "-headers" => 1,
				);
      while (my $row = $df->get_hash()) {
	push @{$rows}, $row;
      }
      @headers = @{$df->headers_raw};
    } else {
      my $pvr = new PostprocessedVariantReport();
      $rows = $pvr->parse("-file" => $infile);
      #  my @headers = @{$pvr->headers()};
      @headers = @{$pvr->get_headers()};
      # use raw user-supplied headers if present, otherwise autodetected
    }
    check_required_snv_indel_headers("-file" => $infile,
				     "-headers" => \@headers);

    push @headers, $FIELD_DBSNP;
    push @headers, $FIELD_NHLBI_FREQ;
    push @headers, $FIELD_COSMIC;
    push @headers, $FIELD_PCGP_COUNT_SOMATIC;
    push @headers, $FIELD_PCGP_DENOMINATOR_SOMATIC;
    push @headers, $FIELD_FOLDX;
    push @headers, @NSFP_FIELDS;
    push @headers, $FIELD_THOUSAND_GENOMES;
    push @headers, $FIELD_OMIM_ID;

    if ($ENABLE_SOMATIC_EXAC) {
      # 1/2018
      my @h;
      foreach my $f (qw(AF AC AN AC_Adj AN_Adj)) {
	push @h, "ExAC_" . $f;
      }
      foreach my $pop (@EXAC_POPULATION_TAGS) {
	foreach my $f (qw(AC AN)) {
	  push @h, sprintf "ExAC_%s_%s", $f, $pop;
	}
      }
      @EXAC_HEADERS = @h;
      push @headers, @h;
      push @headers, @EXAC_COVERAGE_H;
      # old Reason example field example:
      # ExAC_AF=9.415e-06;ExAC_AC=1;ExAC_AN=106210;ExAC_AC_Adj=1;ExAC_AN_Adj=105794;ExAC_pop=AFR:0:8942:0,AMR:0:11216:0,EAS:0:7866:0,FIN:0:6614:0,NFE:0:54070:0,SAS:1:16392:6.101e-05,OTH:0:694:0
    }

    push @headers, $FIELD_PUBMED_IDS;

    push @headers, qw(
		       GSBClass
		       Reason
		       Evidence
		    );

    my @out_rows;

    add_batch_dbsnp(
		    "-rows" => $rows
		   );

    my $dbs_nhlbi;
    if (my $tf = $FLAGS{"tabix-nhlbi"}) {
      add_batch_nhlbi_frequency(
				"-rows" => $rows,
			       );
    } else {
      die "-tabix-nhlbi"
    }

    if (not($FLAGS{"no-cosmic"}) and my $tf = $FLAGS{"tabix-cosmic"}) {
      add_batch_cosmic(
		       "-rows" => $rows,
		      );
    }

    my $sjpi = get_sjpi();

    my $dbnsfp;
    if ($FLAGS{"tabix-dbnsfp"}) {
      add_batch_nsfp(
		     "-rows" => $rows,
		     "-sjpi" => $sjpi,
		    );
    } else {
      $dbnsfp = get_nsfp();
    }

    add_batch_thousand_genomes(
			       "-rows" => $rows,
			      );



    if ($FLAGS{"tabix-gedi"}) {
      add_batch_gedi(
		     "-rows" => $rows
		    );
    }

    if ($ENABLE_SOMATIC_EXAC) {
      # 1/2018: need to move ExAC annotation to somatic since
      # it is replacing NHLBI ESP.  May ultimately use gnomAD
      # if TCGA-subtracted version is available.
      add_batch_exac_columnar(
			      "-rows" => $rows
			     );
      add_batch_exac_coverage_columnar(
				       "-rows" => $rows
				      );
    }

    my $vm_pcgp_somatic = get_vm_pcgp_somatic();
    ram_debug("somatic_snv/vm_pcgp_somatic");
    # FIX ME: move these out of file loop!

    my $vm_committee = get_vm_committee();

    foreach my $row (@{$rows}) {
      dump_die($row, "start variant lookup", 1) if $FLAGS{"debug-start"};

      my $ref_allele = get_field($row, "ReferenceAllele");
      my $mut_allele = get_field($row, "MutantAllele");

      foreach ($ref_allele, $mut_allele) {
	# 9/19/2018: sometimes Excel parsing damage converts "-" to "0" (fail)
	dump_die($row, "numbers not allowed in allele fields, found value \"$_\"") if /\d/;
      }

      # obsolete start: use only if batch tabix versions not available
      add_nhlbi_frequency($row, $dbs_nhlbi) if $dbs_nhlbi;
      add_cosmic($row, $vm_cosmic_somatic) if $vm_cosmic_somatic;
      add_nsfp(
	       "-row" => $row,
	       "-dbnsfp" => $dbnsfp,
	       "-sj" => $row,
	       "-aa" => ($row->{$FIELD_AACHANGE} || ""),
	       # may not be present, e.g. variant in promoter region
	       "-nm" => $row->{$FIELD_REFSEQ},
	       "-preferred" => $sjpi,
	      ) if $dbnsfp;
      # obsolete end

      add_pcgp_somatic_population($row, $vm_pcgp_somatic);
      add_foldx($row);

      my $gene = get_field($row, "GeneName");

      my $omim_ids = "";
      if (my $omim_hits = $omim->find_gene($gene)) {
	$omim_ids = join ",", @{$omim_hits};
      }
      $row->{$FIELD_OMIM_ID} = $omim_ids;

      my ($is_gold, $is_silver);

      if ($SOMATIC_GENE_GSM) {
	$is_gold = $gsm_gold_genes->find($gene);
	# includes primary gold + germline reportable
	$is_silver = $gsm_silver_genes->find($gene);
      } else {
	$is_gold = $gold_genes->{$gene};
	$is_gold = 1 if $GL_REPORTABLE_GENES_GER->{$gene};
	# also include any germline-reportable gene
	$is_silver = $silver_genes->{$gene};
      }


      my $is_indel;
      foreach ($ref_allele, $mut_allele) {
	#      $is_indel = 1 unless /^[ACGT]$/;
	$is_indel = 1 if $_ eq "-";
      }

      my $is_multibase_substitution;
      unless ($is_indel) {
	foreach ($ref_allele, $mut_allele) {
	  $is_multibase_substitution = 1 if length($_) > 1;
	}
      }

      #    confess sprintf "indels not implemented: ref=%s mut=%s", $ref_allele, $mut_allele if $is_indel;

      #    my $aa = sprintf "p.%s", $row->{AAChange} || die "no AAChange";
      my $aa = sprintf "p.%s", $row->{AAChange} || "unknown";
      # may not be present (e.g. promoter)
      my $cosmic_key = join "_", $gene, $aa;
      my $cosmic_snv_hits = $new_cosmic ? undef : $cosmic->{$cosmic_key};
      my $cosmic_indel_hits;
      my $cosmic_indel_hits_cluster;
      if ($is_indel and !$new_cosmic) {
	$cosmic_indel_hits = find_cosmic_indel_hits($cosmic_indels, $row);
	$cosmic_indel_hits_cluster = find_cosmic_indel_hits($cosmic_indels,
							    $row,
							    "-max-codon-distance" => MAX_CODON_DISTANCE_FOR_CLUSTERED_SEARCH,
							    "-gold" => $is_gold,
							    "-silver" => $is_silver
							   );
      }

      my $medal;
      my @reasons;
      if (my $list = $row->{$INTERNAL_FIELD_QUEUE_REASONS}) {
	# from batch tabix
	push @reasons, @{$list};
	delete $row->{$INTERNAL_FIELD_QUEUE_REASONS};
      }

      my $gedi_snv_hits;
      my $gedi_recurrent;

      #    my $variant_class = lc($row->{Class} || die "no Class field");
      my $variant_class = lc($row->{Class} || "");

      #    printf STDERR "gold:%d silver:%d\n", ($is_gold || 0), ($is_silver || 0);

      if ($is_gold or $is_silver) {
	my ($is_ts, $is_recurrent);

	if ($SOMATIC_TS_ONCO_DB) {
	  # new model: single source
	  $is_ts = $ts_onco_db->is_lof("-gene" => $gene);
	  $is_recurrent = $ts_onco_db->is_gof("-gene" => $gene);
	  # TO DO: fatal error if record not found for gold/silver genes
	} elsif (my $gene_class = $gene2class->{$gene}) {
	  # original annotations
	  ($is_ts, $is_recurrent) = parse_gene_class($gene_class);
	} elsif (exists $GL_TS_INFO->{$gene}) {
	  # germline annotations
	  my $info = $GL_TS_INFO->{$gene};
	  if ($info == 0) {
	    $is_recurrent = 1;
	    # gain of function
	  } elsif ($info == 1) {
	    $is_ts = 1;
	    # loss of function
	  } elsif ($info == GERMLINE_ANNOTATION_BOTH_TS_AND_ONCO) {
	    # both gain and loss
	    $is_ts = $is_recurrent = 1;
	  } else {
	    die "unknown ts/onco code $info for $gene";
	  }
	} else {
	  die "no gene2class entry for $gene";
	}

	my $is_truncating;
	if ($variant_class) {
	  # may not be defined (e.g. promoter)
	  $is_truncating = $variant_class_truncating{$variant_class};
	  die "unhandled variant class " . $variant_class unless defined $is_truncating;
	}

	#
	#  PCGP main and recurrent sites:
	#
	$gedi_snv_hits = $row->{$FIELD_GEDI_TABIX};
        $gedi_recurrent = $vm_pcgp_somatic->find_snv("-sj" => $row);
	$gedi_recurrent = $vm_pcgp_somatic->find_literal_variant("-sj" => $row) unless $gedi_recurrent;
	# for SJ-to-SJ lookup, require indel position annotations
	# to match perfectly

	my $cosmic_recurrent;
	if ($new_cosmic) {
	  # assigned in batch tabix
	  $cosmic_recurrent = $row->{$FIELD_COSMIC_RECURRENT};
	} elsif ($cosmic_snv_hits) {
	  my $wanted = $is_gold ? "is_gold" : "is_silver";
	  ($cosmic_recurrent) = grep {$_->{$wanted}} @{$cosmic_snv_hits};
	  # recurrent mutation (COSMIC sub-lists)
	  printf STDERR "checking COSMIC SNV %s, hit=%d\n", $cosmic_key, $cosmic_recurrent ? 1 : 0 if $VERBOSE_COSMIC;
	  #	printf STDERR "COSMIC miss for %s %s: avail=%s %s\n",
	  #	  $cosmic_key, $wanted,
	  #	  join(",", map {$_->{Gene}} @{$cosmic_snv_hits}),
	  #	    join(",", map {$_->{AAchange}} @{$cosmic_snv_hits});
	} elsif ($cosmic_indel_hits) {
	  my $wanted = $is_gold ? "is_gold" : "is_silver";
	  ($cosmic_recurrent) = grep {$_->{$wanted}} @{$cosmic_indel_hits};
	  # recurrent mutation (COSMIC sub-lists)
	  printf STDERR "checking COSMIC indel %s, hit=%d\n", $cosmic_key, $cosmic_recurrent ? 1 : 0 if $VERBOSE_COSMIC;
	}

	my $is_inframe_indel;
	if ($is_indel) {
	  my $aa = $row->{AAChange} || "";
	  if ($variant_class_inframe_indel{$variant_class}) {
	    $is_inframe_indel = 1;
	    confess "say what? fs AA change $aa in $variant_class" if $aa =~ /fs$/i;
	  } else {
	    printf STDERR "WARNING: non-fs AA change, AA=%s class=%s ref=%s var=%s\n", $aa, $variant_class, $ref_allele, $mut_allele unless $aa =~ /fs$/ or $variant_class eq "intron";
	  }
	}

	if ($is_truncating) {
	  #
	  #  Truncation mutation
	  #
	  if ($is_gold) {
	    my $label;
	    if ($is_ts) {
	      $medal = ts_truncating_medal($is_truncating);
	      $label = "tumor suppressor";
	    } elsif ($is_recurrent) {
	      $medal = CLASS_SILVER;
	      $label = "recurrent/oncogene";
	    } else {
	      die;
	    }
	    push @reasons, sprintf "event=truncation in gold gene (%s)", $label;
	  } elsif ($is_silver) {
	    $medal = CLASS_SILVER;
	    my $label;
	    if ($is_ts) {
	      $label = "tumor suppressor";
	    } elsif ($is_recurrent) {
	      $label = "recurrent/oncogene";
	    } else {
	      die;
	    }
	    push @reasons, sprintf "event=truncation in silver gene (%s)", $label;
	  } else {
	    die;
	  }
	} elsif ($is_inframe_indel) {
	  if ($is_gold) {
	    $medal = CLASS_GOLD;
	  } elsif ($is_silver) {
	    $medal = CLASS_SILVER;
	  } else {
	    die;
	  }
	  push @reasons, "event=in-frame indel";
	}

	#
	#  the lookups below can assign medals, but can also add
	#  supplemental information only.  e.g. an in-frame indel is
	#  assigned a medal above, but it may also have recurrency info
	#  which is appended below without assigning a medal.
	#
	if ($cosmic_recurrent or $gedi_recurrent) {
	  #
	  #  recurrent variant (COSMIC or GeDI/PCGP)
	  #
	  unless ($medal) {
	    if ($is_gold) {
	      $medal = CLASS_GOLD;
	    } elsif ($is_silver) {
	      $medal = CLASS_SILVER;
	    } else {
	      die;
	    }
	  }
	  my @where;
	  push @where, "COSMIC" if $cosmic_recurrent;
	  #	push @where, "GeDI" if $gedi_recurrent;
	  push @where, "PCGP" if $gedi_recurrent;
	  push @reasons, sprintf "event=recurrent mutation (%s)", join ",", @where;
	} elsif ($cosmic_snv_hits or $row->{$FIELD_COSMIC} or
		 $gedi_snv_hits) {
	  #
	  # Other known mutation (COSMIC or PCGP)
	  #
	  unless ($medal) {
	    if ($is_gold) {
	      $medal = CLASS_SILVER;
	    } elsif ($is_silver) {
	      $medal = CLASS_BRONZE;
	    } else {
	      die;
	    }
	  }
	  my @where;
	  push @where, "COSMIC" if $cosmic_snv_hits or $row->{$FIELD_COSMIC};
	  push @where, "PCGP" if $gedi_snv_hits;
	  push @reasons, sprintf "event=other known mutation (%s)", join ",", @where;
	} elsif (!$medal) {
	  #
	  # other variant in a gold/silver gene.
	  #
	  $medal = CLASS_BRONZE;
	  push @reasons, "event=other";
	}

	my $class = $is_gold ? "gold" : ($is_silver ? "silver" : "other");
	push @reasons, sprintf "gene_medal=%s", $class;

	my @class;
	push @class, "TS" if $is_ts;
	push @class, "recurrent/oncogene" if $is_recurrent;
	push @reasons, sprintf "gene_class=%s", join ",", @class;
      } else {
	# n/a
	$medal = CLASS_UNKNOWN;
	@reasons = "not a gold/silver gene";
      }

      my $is_silent_intron_utr = 0;
      if ($variant_class eq "silent" or
	  $variant_class eq "intron" or
	  $variant_class eq "utr_5" or
	  $variant_class eq "utr_3") {
	# JZ 9/30/2013:
	# We should NOT include UTR and silent mutation in
	# Gold/Silver/Bronze class
	$medal = CLASS_UNKNOWN;
	#      unshift @reasons, "silent/intron/utr";
	# totally reset reasons field?
	@reasons = ("silent/intron/utr");
	# reset field so other information is not reported, e.g.
	# PCGP.  May have to be revisited if gene class reporting is
	# desired in these cases.  Alternatively, undesirable classes
	# of variants can be filtered from PCGP database as they are
	# loaded.
	$is_silent_intron_utr = 1;
      }

      my ($rm, $irrelevant, $panel_decision, $committee_call_medal);
      ($rm, $panel_decision, $committee_call_medal) =
	committee_check(
			"-vm" => $vm_committee,
			"-row" => $row,
			"-is-silent" => $is_silent_intron_utr,
			"-protect-flag" => \$irrelevant,
			"-reasons" => \@reasons,
		       );
      $medal = request_medal($medal, $rm) if $rm;

      if (0 and $ENABLE_KNOWN_PROMOTER_SITES) {
	die "known promoter site check in somatic not implemented";
      }

      if ($DOWNGRADE_SOMATIC_MEDAL_IF_HIGH_NHLBI and
	  ($medal eq CLASS_GOLD or $medal eq CLASS_SILVER) and
	  not(get_is_rare($row))) {
	push @reasons, sprintf "somatic_popfreq_downgrade_from_%s", $medal;
	$medal = CLASS_BRONZE;
      }

      if ($SOMATIC_GSB_LOCKED_TO_COMMITTEE and $committee_call_medal) {
	# committee call always overrides other logic
	$medal = $committee_call_medal;
      }

      $row->{GSBClass} = $medal || die "no medal";
      #    die unless @reasons;

      #    $row->{Reason} = join ",", @reasons;
      #    $row->{Reason} = join ";", @reasons;
      $row->{Reason} = join_reason_list(\@reasons);

      #
      #  Evidence field: patch in COSMIC PubMed info, if available:
      #
      my @evidence;
      if ($new_cosmic) {
	if (my $ci = $row->{$FIELD_COSMIC_PUBMED}) {
	  my @f = split /,/, $ci;
	  @f = @f[0 .. 4] if @f > 5;
	  push @evidence, sprintf "COSMIC=%s", join ",", @f;
	}
      } elsif ($cosmic_snv_hits or $cosmic_indel_hits) {
	my @pmi;
	my %saw;
	my $set = $cosmic_snv_hits ? $cosmic_snv_hits : $cosmic_indel_hits;
	foreach (@{$set}) {
	  my $pmi = $_->{PubmedInfo} || next;
	  #	$pmi =~ s/\s+$//;
	  $pmi =~ s/\t.*$//;
	  # strip out "removed" section for now

	  next if $saw{$pmi};
	  $saw{$pmi} = 1;
	  push @pmi, $pmi;
	}
	# STILL NEEDS WORK if multiple entries!!
	if (@pmi) {
	  push @evidence, sprintf "COSMIC=%s", join ",", @pmi;
	}
      }

      if ($gedi_snv_hits) {
	my $set = $gedi_snv_hits;

	my %projects = map {$_, 1} grep {$_} map {$_->{official_project}} @{$set};
	die "no projects??" unless %projects;

	if (0) {
	  foreach my $ref (@{$set}) {
	    foreach (sort keys %{$ref}) {
	      printf STDERR "%s = %s\n", $_, $ref->{$_};
	    }
	    print STDERR "\n";
	  }
	  die;
	}
	push @evidence, sprintf "PCGP=%s", join ",", sort keys %projects;
      } elsif ($gedi_recurrent) {
	# recurrent hits use a database dump containing some putative variants,
	# which will not be present in the primary GeDI export
	my %diseases = map {$_, 1} grep {$_} map {$_->{Disease}} @{$gedi_recurrent};
	push @evidence, sprintf "PCGP=%s", join ",", sort keys %diseases;
      }

      my @pmid;
      if (my $pms = $row->{$INTERNAL_FIELD_QUEUE_EVIDENCE}) {
	# PubMed IDs from batch annotations, e.g. ExAC
	push @pmid, @{$pms};
	delete $row->{$INTERNAL_FIELD_QUEUE_EVIDENCE};
      }
      foreach my $ev (@evidence) {
	# extract PMIDs from COSMIC evidence tags
	if ($ev =~ /^COSMIC=(.*)/) {
	  my @things = split /,/, $1;
	  foreach my $thing (@things) {
	    my @f = split /;/, $thing;
	    push @pmid, $f[0];
	  }
	}
      }
      if (my $f = $row->{vep_pubmed}) {
	# extract PMIDs from VEP+ output, if present
	push @pmid, split /,/, $f;
      }

      $row->{$FIELD_PUBMED_IDS} = join ",", sort {$a <=> $b} @{unique_list(\@pmid)};
      # presumably PMIDs are assigned chronologically

      #    $row->{Evidence} = join ":", @evidence;
      $row->{Evidence} = join_evidence_field(\@evidence);
      # delimiter trainwreck
      # WARNING: this field parsed by SiTH! notify erin/jared before modifying

      if ($row->{Reason} =~ /recurrent mutation/i and
	  $row->{Evidence} !~ /\w/) {
	# QC check
	printf STDERR "ERROR: QC FAIL: no evidence entry for medal %s!:\n", $medal;
	foreach (sort keys %{$row}) {
	  printf STDERR "  %s: %s\n", $_, $row->{$_};
	}
      }

      my $sort_score;
      my $state = $row->{State} || "";
      # might not be present for some files
      # e.g. indels pre-review
      if ($state =~ /good/i) {
	$sort_score = 2000;
      } elsif ($state =~ /bad/i) {
	$sort_score = 1000;
      } else {
	$sort_score = 0;
      }

      if ($medal eq CLASS_GOLD) {
	$sort_score += 300;
      } elsif ($medal eq CLASS_SILVER) {
	$sort_score += 200;
      } elsif ($medal eq CLASS_BRONZE) {
	$sort_score += 100;
      }

      printf STDERR "sort score for %s.%s %s %s = %d\n",
	$row->{Chr},
	  get_sj_pos($row),
	    $state,
	      $medal, $sort_score if $VERBOSE_SOMATIC_PROGRESS;

      $row->{sort_score} = $sort_score;

      push @out_rows, $row;
    }

    #  my @sorted = sort {$b->{sort_score} <=> $a->{sort_score} ||
    #		       $a->{Chr} cmp $b->{Chr} ||
    #			 $a->{WU_HG19_Pos} || $b->{WU_HG19_Pos}} @out_rows;

    my @sorted;
    if ($FLAGS{"no-sort"}) {
      @sorted = @out_rows;
    } else {
      @sorted = sort {$b->{sort_score} <=> $a->{sort_score}} @out_rows;
    }

    #
    #  write sorted output:
    #
    my $rpt = new Reporter(
			   "-delimiter" => "\t",
			   "-file" => $outfile,
			   "-labels" => \@headers
			  );
    $rpt->write_headers();
    # 2/2014: force writing header line even if file contains
    # no data rows (per Matt Parker)

    foreach my $row (@sorted) {
      $rpt->end_row($row);
    }
    $rpt->finish();
    print STDERR "report done\n";
  }

  ram_debug("somatic_snv/done");
  log_msg(sprintf "somatic classification done, took %d", time - $somatic_start_time);
}

sub get_field {
  my ($row, $label) = @_;
  my $value = $row->{$label};
  die "no value for $label" unless defined $value;
  return $value;
}

sub flatfile_to_hash {
  return { map {$_, 1} @{read_simple_file($_[0])} };
}

sub parse_gene_class {
  my ($string) = @_;
  my @things = split /\+/, $string;
  my ($is_ts, $is_recurrent);
  foreach my $thing (@things) {
    if ($thing eq "Recur") {
      $is_recurrent = 1;
    } elsif ($thing eq "TS" or $thing eq "TS?") {
      $is_ts = 1;
    } else {
      confess "unknown category string $thing";
    }
  }
  return ($is_ts, $is_recurrent);
}

sub load_cosmic_data {
  my $cosmic_all = $FLAGS{cosmic} || die "-cosmic";
  # currently JZ's filtered version with special PubMed annotations.
  # however this is missing detailed coordinates for indels.
  my $cosmic_gold = $FLAGS{"cosmic-gold"} || die "-cosmic-gold";
  my $cosmic_silver = $FLAGS{"cosmic-silver"} || die "-cosmic-silver";

  printf STDERR "loading COSMIC AA (%s)...\n", $cosmic_all;
  my $df = new DelimitedFile(
			     "-file" => $cosmic_all,
			     "-headers" => 1,
			    );
  my $cosmic_headers = $df->headers_raw();

  my %cosmic2rows;
  # key might might be unique, e.g. EGFR_p.E746_A750delELREA

  my @trim_fields = ("", qw(AAPos Chr Pos1 Pos2 TotalSample TotalVerifiedSample));

  while (my $row = $df->get_hash()) {

    if ($TRIM_COSMIC_INFO) {
      delete @{$row}{@trim_fields};
#      dump_die($row);
    }
    my $key = get_cosmic_key($row);
    #    print STDERR "stash key $key\n";
    push @{$cosmic2rows{$key}}, $row;
  }

  assign_cosmic_gold_silver(\%cosmic2rows, $cosmic_gold, $cosmic_headers, "is_gold");
  assign_cosmic_gold_silver(\%cosmic2rows, $cosmic_silver, $cosmic_headers, "is_silver");

  return \%cosmic2rows;
}

sub assign_cosmic_gold_silver {
  my ($cosmic2rows, $file, $headers, $field, $missing_ok) = @_;
  my $df = new DelimitedFile(
			     "-file" => $file,
			     "-headers" => 0,
			    );

  while (my @row = $df->next()) {
    my %r;
    @r{@{$headers}} = @row;
    my $key = get_cosmic_key(\%r);
    my $list = $cosmic2rows->{$key};
    if ($list) {
      printf STDERR "assigning %s to %s\n", $field, $key if $VERBOSE_COSMIC;
      foreach my $r (@{$list}) {
	$r->{$field} = 1;
      }
    } elsif (!$missing_ok) {
      confess "can't find $key in master list" if not($list) and not($missing_ok);
      # indel mode won't have all entries
    }

  }
}

sub get_cosmic_key {
  my ($row) = @_;
  #  return join "_", @{$row}{qw(Gene AAchange)};
  my $gene = $row->{Gene} || $row->{"Gene name"};
  die "gene name undef" unless defined $gene;

  my $aa = $row->{AAchange} || $row->{"Mutation AA"};
  #  confess "can't identify AA: " . join ",", map {$_ . "=" . $row->{$_}} sort keys %{$row} unless defined $aa;
  # either JZ's version or raw COSMIC
  $aa = "_undef_" unless defined $aa;
  # some appear undef, e.g. FAM54B
  return join "_", $gene, $aa;
}

sub get_gedi_key {
  my ($row) = @_;
  # my $key = join ".", @{$row}{qw(Chr Pos RefAllele MutAllele)};
  my $chr = $row->{Chr};
  $chr =~ s/^chr//;
  my $pos = get_sj_pos($row) || die;
  my $ra = $row->{ReferenceAllele} || die;
  my $va = $row->{MutantAllele} || die;

  #  my $key = join ".", $chr, @{$row}{qw(WU_HG19_Pos ReferenceAllele MutantAllele)};
  my $key = join ".", $chr, $pos, $ra, $va;

  return $key;
}

sub get_gedi_recurrent_restrict {
  die "-gedi-recurrent-restrict" unless @GEDI_RECURRENT_RESTRICT;
  die "-gedi-recurrent-restrict uses multiple files" unless @GEDI_RECURRENT_RESTRICT > 1;
  # JZ 3/7/2014:
  # "Please only include the genes listed in the following two files:
  # cancer_gene_SJ_mutation.txt non_cancer_gene_SJ_mutation_mod.txt
  # for look up."
  my %wanted;
  foreach my $fr (@GEDI_RECURRENT_RESTRICT) {
    open(FRTMP, $fr) || die;
    while (<FRTMP>) {
      chomp;
      my @f = split /\t/, $_;
      die unless @f == 2;
      my ($gene) = @f;
      $wanted{$gene} = 1;
    }
    close FRTMP;
  }
  return \%wanted;
}

sub get_outfile {
  my ($infile) = @_;
  my $outfile;
  if ($FLAGS{"outfile-fq"}) {
    $outfile = $infile . ".medals.tab";
  } else {
    $outfile = basename($infile) . ".medals.tab";
  }
  return $outfile;
}

sub load_cosmic_indels {
  #
  #  load all indels observed in COSMIC
  #  2 use cases:
  #   - perfect match using genomic coordinates
  #   - fuzzy match using codon numbers for expanded "cluster" annotation
  #
  my $fn = $FLAGS{"cosmic-cleaned"} || die "-cosmic-cleaned";
  printf STDERR "loading COSMIC indels (%s)...\n", $fn;
  my $df = new DelimitedFile(
			     "-file" => $fn,
			     "-headers" => 1,
			    );
  my %indels;
  my %indel_descs = (
		     "Substitution - Missense" => 0,
		     "Substitution - Nonsense" => 0,
		     "Substitution - coding silent" => 0,
		     "Unknown" => 0,
		     "Complex" => 0,
		     "Nonstop extension" => 0,
		     "Complex - compound substitution" => 0,
		     "No detectable mRNA/protein" => 0,

		     "Whole gene deletion" => 1,
		     # regional? disqualify?

		     "Complex - insertion inframe" => 1,
		     "Complex - deletion inframe" => 1,
		     "Complex - frameshift" => 1,
		     # ???

		     "Insertion - Frameshift" => 1,
		     "Insertion - In frame" => 1,
		     "Deletion - Frameshift" => 1,
		     "Deletion - In frame" => 1
		    );

  my @cosmic_trim;

  while (my $row = $df->get_hash()) {
    my $desc = $row->{"Mutation Description"} || die;
    my $is_indel = $indel_descs{$desc};
    die "unhandled description $desc" unless defined $is_indel;

    if ($is_indel) {
      #      print STDERR "fix me!\n";
      #      die;

      if ($TRIM_COSMIC_INFO) {
	unless (@cosmic_trim) {
	  my %all = map {$_, 1} keys %{$row};
	  my %keep = map {$_, 1} ("Gene name", "Mutation AA");
	  foreach (keys %keep) {
	    delete $all{$_};
	  }
	  @cosmic_trim = keys %all;
	}
	delete @{$row}{@cosmic_trim};
      }

      my $gene = $row->{"Gene name"} || die;
      my $aa = $row->{"Mutation AA"} || die;
      my $aa_raw = $aa;
      my $type;
      my $nt_size;

      #      printf STDERR "indel desc:%s CDS:%s\n", $desc, $row->{"Mutation CDS"};

      my $aap = new AAParser();
      my $parsed_ok;
      if ($aap->parse($aa)) {
	$parsed_ok = 1;
      } else {
	printf STDERR "ERROR: can't identify codon in COSMIC AA $aa $desc, skipping\n" if $VERBOSE_COSMIC and $desc !~ /whole gene deletion/i;
      }

      if ($desc eq "Deletion - Frameshift" or
	  $desc eq "Deletion - In frame" or
	  $desc eq "Complex - deletion inframe",
	 ) {
	$type = "D";
      } elsif ($desc eq "Insertion - Frameshift" or
	       $desc eq "Insertion - In frame" or
	       $desc eq "Complex - insertion inframe",
	      ) {
	$type = "I";
      } else {
	printf STDERR "FIX ME: skipping indel %s\n", $desc if $VERBOSE_COSMIC;
      }

      if ($type and $parsed_ok) {
	#	printf STDERR "codon stash of %s\n", join " ", $type, $gene, $codon_number;
	my $ckey = get_cosmic_key($row);
	#	printf STDERR "saving ckey %s\n", $ckey;
	push @{$indels{protein_key}{$ckey}}, $row;
	# for mapping between raw COSMIC and JZ's cooked silver/gold lists,
	# used below

	my $cs = $aap->codon_start || die;
	my $ce = $aap->codon_end || die;

	foreach my $codon_number ($cs .. $ce) {
	  push @{$indels{gene2codon}{$type}{$gene}{$codon_number}}, $row;
	}
      }
    }
  }

  #
  #  patch in gold/silver annotation status from JZ's files:
  #
  printf STDERR "begin gold/silver assignment for COSMIC indels:\n" if $VERBOSE_COSMIC;
  my $cosmic_all = $FLAGS{cosmic} || die "-cosmic";
  # currently JZ's filtered version with special PubMed annotations.
  # however this is missing detailed coordinates for indels.
  my $cosmic_gold = $FLAGS{"cosmic-gold"} || die "-cosmic-gold";
  my $cosmic_silver = $FLAGS{"cosmic-silver"} || die "-cosmic-silver";
  $df = new DelimitedFile(
			  "-file" => $cosmic_all,
			  "-headers" => 1,
			 );
  my $cosmic_headers = $df->headers_raw();
  # just get the headers for parsing JZ's gold/silver excerpt files

  assign_cosmic_gold_silver($indels{protein_key}, $cosmic_gold, $cosmic_headers, "is_gold", 1);
  assign_cosmic_gold_silver($indels{protein_key}, $cosmic_silver, $cosmic_headers, "is_silver", 1);
  # this will assign gold/silver status from JZ's custom reports
  # to the main row
  printf STDERR "end gold/silver assignment for COSMIC indels\n" if $VERBOSE_COSMIC;

  #
  #  patch in PubmedInfo field from JZ's report:
  #
  patch_pubmed_info($indels{protein_key});

  return \%indels;
}

sub find_cosmic_indel_hits {
  my ($cosmic_indels, $row, %options) = @_;
  # $cosmic_indels{gene2codon}{$type}{$gene}{$codon_number}

  my $max_codon_distance = $options{"-max-codon-distance"} || 0;

#  foreach (sort keys %{$row}) {
#    printf STDERR "%s: %s\n", $_, $row->{$_};
#  }

  my $type;
  if ($row->{ReferenceAllele} eq "-") {
    # insertion
    $type = "I";
  } elsif ($row->{MutantAllele} eq "-") {
    # deletion
    $type = "D";
  } else {
    confess sprintf "WTF: not an indel ref=%s var=%s", @{$row}{qw(ReferenceAllele MutantAllele)};
  }
  my $gene = $row->{GeneName} || "";
  my $aa = $row->{AAChange} || "";
  my ($codon_start, $codon_end);

  my @rough_with_duplicates;

  my $aap = new AAParser();
  if ($aap->parse($aa)) {
    $codon_start = $aap->codon_start();
    $codon_end = $aap->codon_end();
    if ($VERBOSE_COSMIC) {
      printf STDERR "parsed SJ AA %s: %s %s\n", $aa, $codon_start, $codon_end;
      printf STDERR "codon TARGET lookup for %s/%s %s-%s\n", $type, $gene, $codon_start, $codon_end;
    }

    my @search;
    if ($max_codon_distance) {
      @search = $codon_start - $max_codon_distance .. $codon_end + $max_codon_distance;
    } else {
      @search = $codon_start .. $codon_end;
    }

    foreach my $i (@search) {
      my $set = $cosmic_indels->{gene2codon}{$type}{$gene}{$i};
      push @rough_with_duplicates, @{$set} if $set;
    }
  } else {
    printf STDERR "ERROR: can't parse SJ AA %s for lookup in %s\n", $aa, $gene if $VERBOSE_COSMIC;
  }

  return \@rough_with_duplicates;
  #
  #  attempt to remove duplicates??
  #
}

sub patch_pubmed_info {
  my ($indels) = @_;
  my $cosmic_all = $FLAGS{cosmic} || die "-cosmic";
  # JZ's filtered version with special PubMed annotations.
  printf STDERR "loading %s...\n", $cosmic_all;
  my $df = new DelimitedFile(
			     "-file" => $cosmic_all,
			     "-headers" => 1,
			    );

  while (my $row = $df->get_hash()) {
    my $key = get_cosmic_key($row);
    if (my $hits = $indels->{$key}) {
      foreach my $r2 (@{$hits}) {
	$r2->{PubmedInfo} = $row->{PubmedInfo};
      }
    }
  }
}

sub load_config_file {
  # pre-populate command line params from configuration file
  my @things;
  for (my $i = 0; $i < @ARGV; $i++) {
    if ($ARGV[$i] eq "-config") {
      my $config_file = $ARGV[$i + 1];
      die "where is $config_file" unless -s $config_file;
      open(CF, $config_file) || die;
      while (<CF>) {
	chomp;
	next if /^\#/;
	push @things, split /\s+/, $_;
      }
    }
  }
  #  push @ARGV, @things if @things;
  unshift @ARGV, @things if @things;
  # put config file entries first in @ARGV rather than appending.
  # This allows a user's manual specification to override the config
  # file specification.
}


sub init_infiles {
  # detect processing type from filename
  my @raw;
  if (my $lf = $FLAGS{"in-list"}) {
    @raw = @{read_simple_file($lf)};
    die "can't use -out with -in-list" if $FLAGS{out};
  } elsif (my $f = $FLAGS{in}) {
    @raw = $f;
  }

  #
  # TO DO: peek at headers/data to determine file type
  #
  foreach my $fn (@raw) {
    my $bn = basename($fn);
    if ($bn =~ /_indel/ or $bn =~ /tier\d+_putative/) {
      # somatic SNVs or indels
      push @INPUT_SNV_INDEL_FILES, $fn;
    } elsif ($bn =~ /CONSERTING/) {
      # CONSERTING (CNV)
      push @INPUT_CNV_FILES, $fn;
    } elsif ($bn =~ /\.predSV/) {
      # CREST (SV)
      push @INPUT_SV_FILES, $fn;
    } elsif ($bn =~ /event_fusion/) {
      # FusionBuilder (SV)
      push @INPUT_SV_FILES, $fn;
    } elsif ($bn =~ /germline/) {
      # germline SNVs/indels
      push @INPUT_GL_FILES, $fn;
    } else {
      die "don't recognize filename format for $fn, maybe try -in-snv-indel/-in-cnv/in-sv etc.?";
    }
  }

  foreach my $ref (
		   ["in-snv-indel", \@INPUT_SNV_INDEL_FILES],
		   ["in-cnv", \@INPUT_CNV_FILES],
		   ["in-sv", \@INPUT_SV_FILES],
		   ["in-gl", \@INPUT_GL_FILES],
		  ) {
    my ($flag, $array) = @{$ref};
    if (my $fn = $FLAGS{$flag}) {
      my $files = read_simple_file($fn);
      push @{$array}, @{$files};
    }
  }

  die "nothing to process" unless @INPUT_SNV_INDEL_FILES or
    @INPUT_CNV_FILES or
      @INPUT_SV_FILES or
	@INPUT_GL_FILES;
}

sub run_svs_old {
  # DELETE ME
  my $svc = $FLAGS{"sv"} || die "specify -sv";
  my $lines = read_simple_file($svc);
  # primary SVs

  my $gene2class = load_gene_class_files();
  # tumor suppressor/oncogene annotations
  # NOTE: may need revision, there is an additional list for germline

  # manually-curated SVs:
  my $svc_manual = $FLAGS{"sv-manual"} || die "-sv-manual";
  line_delimiter_qc("-file" => $svc_manual, "-delimiter" => "\t");
  printf STDERR "SV manual: %s\n", $svc_manual;
  my $df = new DelimitedFile(
			     "-file" => $svc_manual,
			     "-headers" => 1,
			    );
  while (my $row = $df->get_hash()) {
    my ($gene1, $gene2, $submitter) = @{$row}{qw(pair1 pair2 contact)};
    push @{$lines}, sprintf "%s\n", join "\t", $gene1, $gene2, $submitter;
    # reformat
  }

  my @reference_svs;
  #
  # load reference SV list to check against:
  #
  my %sv_genes;
  foreach my $line (@{$lines}) {
    my @f = split /\t/, $line;
    die scalar @f unless @f == 3;
    my ($gene1, $gene2, $stuff) = @f;
    die unless $gene1 or $gene2;
    my %sv;
    $sv{breakpoints} = [];
    foreach ($gene1, $gene2) {
      if (defined($_) and /\w/) {
	push @{$sv{breakpoints}}, $_;
	# each reference SV:
	# - may contain one or two gene breakpoints
	# - each breakpoint refers to exactly one gene
	$sv_genes{$_} = 1;
      }
    }
    $sv{breakpoint_count} = scalar @{$sv{breakpoints}};
    push @reference_svs, \%sv;
  }

  my $gold_genes_patched = flatfile_to_hash($FLAGS{"gold-genes-mapped-to-fb"} || die "-gold-genes-mapped-to-fb");
  # gold gene list for SNV/indel patched to FusionBuilder gene symbol space

  my $gl_reportable = read_simple_file($FLAGS{"gl-reportable-genes-refflat"} || die "-gl-reportable-genes-refflat");
  foreach my $g (@{$gl_reportable}) {
    # 3/2016: also include germline-reportable genes
    $gold_genes_patched->{$g} = 1;
  }

  $df = new DelimitedFile(
			  "-file" => ($FLAGS{"sv-silver-genes-fb"} || die "-sv-silver-genes-fb"),
			  "-headers" => 1,
			 );
  my %silver_genes;
  while (my $row = $df->get_hash()) {
    my $gene = $row->{gene} || die;
#    die sprintf "silver gene %s is already in gold gene list!", $gene if $gold_genes_patched->{$gene};
    # obsolete: can't assume this now given use of germline reportable list
    die sprintf "silver gene %s is already in SV list!", $gene if $sv_genes{$gene};
    # ensure genes aren't already mentioned in other lists

    $silver_genes{$row->{gene} || die} = 1;
  }

  my ($gsm_gold_genes, $gsm_silver_genes);
  if ($SV_GSM) {
    $gsm_gold_genes = new_gsm_lite();
    foreach (keys %{$gold_genes_patched}) {
      $gsm_gold_genes->add_gene("-gene" => $_);
    }

    $gsm_silver_genes = new_gsm_lite();
    foreach (keys %silver_genes) {
      $gsm_silver_genes->add_gene("-gene" => $_);
    }
  }


  #
  #  check to make sure all gene symbols used in config actually
  #  appear in annotation source
  #
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";

  my $rf = $FLAGS{"refflat-fb"};
  unless ($rf) {
    $rf = $config_genome->{REFSEQ_NM_REFFLAT} || die "no config REFSEQ_NM_REFFLAT";
  }
  open(RFTMP, $rf) || die;
  my %refflat;
  while (<RFTMP>) {
    chomp;
    my @f = split /\t/, $_;
    my $gene = clean_sharp_gene_symbol($f[0]);
    $refflat{$gene} = 1 if $gene;
  }

  my $broken_symbols;
  my %all = (%{$gold_genes_patched}, %silver_genes, %sv_genes);
  foreach my $sym (sort keys %all) {
    unless ($refflat{$sym}) {
      printf STDERR "ERROR: SV config uses symbol %s but not found in %s!\n", $sym, $rf;
      $broken_symbols = 1;
    }
  }
  printf STDERR "ERROR: problematic SV configuration!  Classification may not work for the genes above.\n" if $broken_symbols;

  my %FB_USAGE = map {$_, 1} qw(
				 GENIC
				 TRUNCATING
				 INTERGENIC
				 HALF_INTERGENIC
				 INTRONIC
				 DUPLICATE
				 UNUSABLE

				 INVERTED_REPEAT
			      );
  # known Usage fields from FusionBuilder
  # http://hc-wiki.stjude.org/display/compbio/Fusion+Builder

  foreach my $infile (@INPUT_SV_FILES) {
    printf STDERR "processing %s...\n", $infile;
    my $df = new DelimitedFile(
			       "-file" => $infile,
			       "-headers" => 1
			      );

    my @out_rows;
    while (my $row = $df->get_hash()) {
      my $is_truncating;
      my @breakpoints;
      # terrible name: actually "sides" of the event / fusion partner lists
      if (exists $row->{$FUSION_BUILDER_COLUMN}) {
	# FusionBuilder format:
	# - either one or two partner lists
	#   (1 indicates internal rearrangement?)
	# - each side of breakpoint may contain one or more genes
	my @sides = split /_/, $row->{$FUSION_BUILDER_COLUMN};
	my %saw;
	foreach my $thing (@sides) {
	  die unless $thing =~ /\w/;
	  next if $saw{$thing};
	  $saw{$thing} = 1;
	  # special case for internal rearrangements, e.g.
	  # PAX5_PAX5, see below in CREST example
	  my %genes;
	  foreach my $gene (split /,/, $thing) {
	    $gene = clean_sharp_gene_symbol($gene);
	    $genes{$gene} = 1;
	  }
	  push @breakpoints, \%genes;
	}
	die "no FusionBuilder Usage column present" unless exists $row->{Usage};
	my $usage = $row->{Usage} || "";
	# may be blank in very rare cases
	# e.g. /rgs01/resgen/prod/tartan/runs/crest-post/fnEMxF4l/intmd/SJRB012408_X2_G1/SJRB012408_X2_G1.counts
	#	die "unknown FusionBulder Usage field value $usage" unless $FB_USAGE{$usage};
	# probably too finicky
	$is_truncating = $usage eq "TRUNCATING";
      } else {
	# CREST format:
	# - either one or two gene columns (fusion partners)
	# - each column may contain one or more genes
	die "where is GeneA?" unless exists $row->{GeneA};
	my $gene_a = $row->{GeneA};
	my $gene_b = $row->{GeneB};
	my %saw;
	foreach my $thing ($gene_a, $gene_b) {
	  next unless $thing;
	  next if $saw{$thing};
	  # special case: if both columns refer to the same gene,
	  # treat as a single entry. This will enable perfect matches
	  # to reference SVs (would otherwise receive a silver match
	  # b/c breakpoint count doesn't match [putative:2 reference:1])
	  # e.g. EBF1 in
	  # /nfs_exports/genomes/1/projects/.allType/ClinicalPilot/.allTumor/SJPHALL/NextGen/WholeGenome/StructuralVariation/CREST/SJPHALL005_D_G.predSV
	  # contains examples of BOTH types:
	  # - EBF1 in GeneA/B, the other blank
	  # - EBF1 in BOTH GeneA and GeneB.
	  $saw{$thing} = 1;

	  my @g = split /,/, $thing;
	  my %bucket;
	  foreach (@g) {
	    next unless /\w/;
	    $_ = clean_sharp_gene_symbol($_);
	    $bucket{$_} = 1;
	  }
	  push @breakpoints, \%bucket if %bucket;
	}
	if (0) {
	  foreach (sort keys %{$row}) {
	    printf "%s: %s\n", $_, $row->{$_};
	  }
	}
	#	die "WTF: no breakpoints in CREST row" unless @breakpoints;
	# never mind: rows might not have gene annotations

	die "ERROR: don't know how to get truncation status from CREST";
	# may not be needed, I think CREST may have been just the
	# initial/test data
      }

      #
      # find all matching SVs:
      #
      my @sv_hits;
      foreach my $sv_ref (@reference_svs) {
	my @found_i;
	my @fuzzy_i;
	foreach my $end_gene (@{$sv_ref->{breakpoints}}) {
	  # one gene per reference breakpoint
	  for (my $i = 0; $i < @breakpoints; $i++) {
	    my $genes = $breakpoints[$i];
	    if ($genes->{$end_gene}) {
	      # reference breakpoint gene found in local breakpoint
	      push @{$found_i[$i]}, $end_gene;
	      $fuzzy_i[$i] = 1 if scalar keys %{$genes} > SV_MAX_BREAKPOINT_GENES_FOR_GOLD;
	    }
	  }
	}

	my $local_end_hit_count = scalar grep {$_} @found_i;
	# count of ENDS hit, rather than gene count
	my @reasons;
	my ($is_fuzzy) = grep {$_} @fuzzy_i;
	push @reasons, "large_gene_list" if $is_fuzzy;

	my @match_genes = map {@{$_}} grep {$_} @found_i;

	my %hit;
	$hit{sv} = $sv_ref;
	$hit{match_genes} = \@match_genes;
	$hit{reasons} = \@reasons;

	if (@breakpoints == $sv_ref->{breakpoint_count} and
	    # local SV has the same number of breakpoints as the reference SV
	    $local_end_hit_count == $sv_ref->{breakpoint_count}
	    # all genes in the reference SV found in local SV,
	    # and all local breakpoints were covered.

	    # note this case should NOT match the example below,
	    # because all the gene symbols are found in GeneA:

	    # reference SV:
	    # /nfs_exports/genomes/1/projects/ClinicalSeq/FusionGenes/SVCheck.txt contains this SV involving C11orf95 and MAML2:
	    # C11orf95        MAML2   EPD;COSMIC
	    #
	    # In the CREST file /nfs_exports/genomes/1/projects/.allType/ClinicalPilot/.allTumor/SJETV/NextGen/WholeGenome/StructuralVariation/CREST/SJETV092_D_G.predSV there is an entry where GeneA contains a long list of genes including both C11orf95 and MAML2, however GeneB is blank.
	   ) {
	  # - local SV has the same number of breakpoints as the reference SV
	  # - both genes in reference SV found

	  my $handled;
	  if ($sv_ref->{breakpoint_count} == 1 and $is_truncating) {
	    # JZ 5/29/2014:
	    # If any single-gene marked as Recurrent, if it is a
	    # truncation, mark them as silver. For those marked as
	    # TS+Recurr, treat them as non-Recur.
	    my $gene = $sv_ref->{breakpoints}->[0];
	    my $class = $gene2class->{$gene};
	    if ($class) {
	      my ($is_ts, $is_recurrent) = parse_gene_class($class);
	      if ($is_recurrent and not($is_ts)) {
		# perfect match, but a truncation in an oncogene
		# which is less interesting than one in a tumor suppressor
		unshift @reasons, ("perfect_match", "oncogene_truncation");
		$hit{medal} = CLASS_SILVER;
		$hit{sort_score} = 6;
		$handled = 1;
	      }
	    }
	  }

	  unless ($handled) {
	    unshift @reasons, "perfect_match";
	    if ($is_fuzzy) {
	      $hit{medal} = CLASS_SILVER;
	      $hit{sort_score} = 5;
	      # pretty good because all genes hit,
	      # but less compelling because there are a lot of genes
	    } else {
	      $hit{medal} = CLASS_GOLD;
	      $hit{sort_score} = 10;
	    }
	  }

	  push @sv_hits, \%hit;
	} elsif (@match_genes) {
	  # a single gene matches, others don't (putative or reference)
	  $hit{sv} = $sv_ref;
	  unshift @reasons, "imperfect_match";
	  $hit{medal} = CLASS_SILVER;
	  $hit{sort_score} = 4;
	  # minor match
	  push @sv_hits, \%hit;
	}

	# my @hit_genes = grep {$sv_ref->{$_}} keys %genes;
	# my $hit_count = scalar @hit_genes;

	# if ($hit_count) {
	#   my $sv_gene_count = scalar keys %{$sv};
	#   my $gene_count_matches = scalar(keys %genes) == $sv_gene_count;
	#   my %hit;
	#   $hit{sv} = $sv_ref;
	#   $hit{match_genes} = \@hit_genes;
	#   if ($sv_gene_count == $hit_count and $gene_count_matches) {
	#     $hit{medal} = CLASS_GOLD;
	#     $hit{sort_score} = 10;
	#   } else {
	#     $hit{medal} = CLASS_SILVER;
	#     $hit{sort_score} = $gene_count_matches ? 5 : 1;
	#   }
	#   push @sv_hits, \%hit;
	# }
      }

      #
      #  check for matches to gold gene list:
      #
      my @search = map {keys %{$_}} @breakpoints;
      # all
      #      my @gold_hits = grep {$gold_genes->{$_}} @search;
      # FAIL: some of these not compatible with genes reported by FusionBuilder

      my @gold_hits;
      if ($SV_GSM) {
#	@gold_hits = grep {$gsm_gold_genes->find($_)} @search;
	foreach my $g (@search) {
	  my $hit = $gsm_gold_genes->find($g);
	  push @gold_hits, $hit if $hit;
	  # report symbol we're officially looking for rather than
	  # user symbol, which may be previous/alias
	}
      } else {
	@gold_hits = grep {$gold_genes_patched->{$_}} @search;
      }
      if (@gold_hits) {
	push @sv_hits, {
			"medal" => CLASS_SILVER,
			"sort_score" => 3,
			"reasons" => [ "gold_gene" ],
			"match_genes" => unique_ordered_list(\@gold_hits)
		       };
	# create new record so as not to interact with %hits above,
	# which may have already been added
      }

      #
      #  silver gene list:
      #
      my @silver_hits;
      if ($SV_GSM) {
	foreach my $g (@search) {
	  my $hit = $gsm_silver_genes->find($g);
	  push @silver_hits, $hit if $hit;
	  # report symbol we're officially looking for rather than
	  # user symbol, which may be previous/alias
	}
      } else {
	@silver_hits = grep {$silver_genes{$_}} @search;
      }

      if (@silver_hits) {
	push @sv_hits, {
			"medal" => CLASS_SILVER,
			"sort_score" => 2,
			"reasons" => [ "silver_gene" ],
			"match_genes" => unique_ordered_list(\@silver_hits)
		       };
	# create new record so as not to interact with %hits above
      }

      my $medal;
      my @gene_hits;
      my $reason = "";

      if (@sv_hits) {
	@sv_hits = sort {$b->{sort_score} <=> $a->{sort_score}} @sv_hits;
	my $best = $sv_hits[0];
	@gene_hits = @{$best->{match_genes}};
	printf STDERR "dump of %d SV hits:", scalar(@sv_hits);
	foreach my $ref (@sv_hits) {
	  my @ref_sides;
	  if ($ref->{sv}) {
	    @ref_sides = @{$ref->{sv}{breakpoints}};
	  } else {
	    @ref_sides = "gold_gene";
	  }
	  printf STDERR " %s (%s)",
	    $ref->{medal}, join ",", @ref_sides;
	}
	print STDERR "\n";
	$medal = $best->{medal};
	$row->{sort_score} = $best->{sort_score};
	$reason = join ",", @{$best->{reasons}};
      } else {
	$medal = "Unknown";
	$row->{sort_score} = 0;
      }

      $row->{GSBClass} = $medal;
      $row->{Reason} = $reason;
      $row->{Evidence} = join ",", @gene_hits;
      push @out_rows, $row;
    }

    my $rpt = $df->get_reporter(
				"-file" => get_outfile($infile),
				"-extra" => [
					     qw(
						 GSBClass
						 Reason
						 Evidence
					      )
					    ]
			       );
    $rpt->write_headers();
    # force header line even if empty file
    foreach my $row (sort {$b->{sort_score} <=> $a->{sort_score}} @out_rows) {
      $rpt->end_row($row);
    }
    $rpt->finish();
  }
  printf STDERR "SVs done\n";
}

sub general_label {
  my ($is_focal, $is_deletion) = @_;
  my @things;
  push @things, $is_focal ? "focal" : "non_focal";
  push @things, $is_deletion ? "deletion" : "amplification";
  return @things;
}

sub run_germline {
  my $start_time = time;
  log_msg(sprintf "germline classification: starting on %s", hostname());
  ram_debug("germline_snv/start");

  ram_debug("gold gene GSM start");
  my $gsm_gold_genes = new_gsm_lite();

  my $ts_onco_db;
  ram_debug("TSOncoDB start");
  if ($GERMLINE_TS_ONCO_DB) {
    $ts_onco_db = new TSOncoDB(
			       "-gsm" => new_gsm_lite()
			      );
  }
  ram_debug("TSOncoDB end");

  ram_debug("ACMG PVS1/start");
  my $acmg_pvs1_caution = get_acmg_pvs1_caution_regions();
  ram_debug("ACMG PVS1/end");

  my $vm_private = get_vm_private();

  my $exac_vcf2tab = $FLAGS{"exac-vcf2tab"} || die "-exac-vcf2tab";
  # new version: cooked
  my $exac_coverage = $FLAGS{"exac-coverage"} || die "-exac-coverage";

  my $tf_exac;
  unless ($exac_vcf2tab) {
    # old version
    my $exac_fn = $FLAGS{"exac"} || die "-exac";
    $tf_exac = new TabixFile("-file" => $exac_fn);
    $tf_exac->wiggle_bases(INDEL_NUCLEOTIDE_NT_FUZZY_MATCH);
  }

  my $use_original_ranges = 1;

  #
  #  gene list setup:
  #
  #  my $gl_gold_genes = flatfile_to_hash($FLAGS{"gl-gold-cancer"} || die "-gl-gold-cancer");
#  my $gl_gold_ranges_file = $FLAGS{"gl-gold-cancer-ranges"} || die "-gl-gold-cancer-ranges";
  my $gl_gold_ranges_file;
  if ($use_original_ranges) {
    # if we can use the original file, we can obsolete -gold-cancer-ranges-ger
    printf STDERR "using raw ranges file\n";
    $gl_gold_ranges_file = $FLAGS{"gl-gold-cancer-ranges"} || die "-gl-gold-cancer-ranges";
  } else {
    # use version standardized to GENE_EXON_REGION (used in SNV annotation)
    printf STDERR "using GER ranges file\n";
    $gl_gold_ranges_file = $FLAGS{"gl-gold-cancer-ranges-ger"} || die "-gl-gold-cancer-ranges-ger";
  }

  my $grf_gold = new GenomicRangeFinder();
  my $gl_gold_genes = {};

  add_custom_gl_genes(
		      "-grf-gold" => $grf_gold,
		      "-gl-gold-genes" => $gl_gold_genes,
		      "-gsm-gold-genes" => $gsm_gold_genes
		     );

  #
  #  2/2014: for gold genes, use genomic regions as a backup
  #  in case there is ever an inconsistency in HUGO symbols
  #
  unless ($FLAGS{"gl-no-reviewable"}) {
    # disable for some projects, e.g. ALS
    open(GLTMP, $gl_gold_ranges_file) || die;
    while (<GLTMP>) {
      chomp;
      my @f = split /\t/, $_;
      die unless @f == 2;
      my ($gene, $loc) = @f;
      $gl_gold_genes->{$gene} = 1
	unless $gene eq "TP53" and $FLAGS{"gl-grf-debug"};
      $gsm_gold_genes->add_gene("-gene" => $gene);

      $grf_gold->add(
		     "-range" => $loc,
		     "-value" => {
				  "gene" => $gene,
				  "location" => $loc
				 },
		    );
    }

    unless ($use_original_ranges) {
      my $gl_reviewable = read_simple_file($FLAGS{"gl-reviewable-genes-ger"} || die "-gl-reviewable-genes-ger");
      foreach my $gene (@{$gl_reviewable}) {
	$gl_gold_genes->{$gene} = 1;
	# backup: also use official mapping of reviewable genes to
	# GENE_EXON_REGION annotation used in cooked COSMIC
      }
    }
  }

  ram_debug("germline_snv/before PCGP germline");
  my $vm_pcgp_germline = get_vm_pcgp_germline();
  ram_debug("germline_snv/after PCGP germline");

  my $vm_lovd_apc = get_vm_generic_aa(
				      "-label" => "APC",
				      "-flag" => "gl-apc-flatfile",
				     );

  my $vm_lovd_msh2 = get_vm_generic_aa(
				       "-label" => "MSH2",
				       "-flag" => "gl-msh2-flatfile",
				      );


  ram_debug("germline_snv/before NHGRI");
  my $vm_nhgri_brca1 = parse_nhgri("-gene" => "BRCA1");
  my $vm_nhgri_brca2 = parse_nhgri("-gene" => "BRCA2");
  ram_debug("germline_snv/after NHGRI");

  my $vm_committee = get_vm_committee();
  #
  #  variant-level calls made by committee
  #  (these trump all automated calls)
  #

  #
  #  oncogenic/hotspot sites in COSMIC:
  #
  my $vm_cosmic_hotspots = get_vm_cosmic_hotspots();
  ram_debug("germline_snv/after COSMIC hotspots");

  #
  #  RB1 database:
  #
  my $vm_rb1 = get_vm_rb1();
  ram_debug("germline_snv/after RB1");

  my $vm_bad = get_vm_bad();
  # blacklisted variants: do not call medals for these
  # - will probably be replaced by formal committee review db

  my $ts_genes;
  if ($GERMLINE_TS_ONCO_DB) {
    $ts_genes = {};
  } else {
    $ts_genes = get_gl_ts_list();
    # tumor suppressor genes (cancer analysis).
  }
  add_custom_ts(
		"-ts-genes" => $ts_genes,
		"-ts-onco-db" => $ts_onco_db
	       );
  # for custom projects, whether or not a truncating variant should recieve a gold medal.

  #
  #  HGMD DM database:
  #
  ram_debug("germline_snv/before HGMD");
  my $vm_hgmd = get_vm_hgmd();
  ram_debug("germline_snv/after HGMD");

  #
  #  COSMIC (gold truncating for now):
  #
  my $vm_cosmic;
  unless ($FLAGS{"tabix-cosmic"}) {
    $vm_cosmic = parse_cosmic(
			      "-file" => ($FLAGS{"cosmic-cleaned"} || die "-cosmic-cleaned"),
			      "-genes" => $gl_gold_genes
			      # NOTE: gene list must be localized to
			      # COSMIC annotation, i.e. GENE_EXON_REGION
			     );
  }


  #
  #  ARUP MEN2 (i.e. RET)
  #
  ram_debug("germline_snv/before ARUP");
  my $vm_arup = parse_arup();
  ram_debug("germline_snv/after ARUP");

  #
  #  umd.be locus-specific databases:
  #
  ram_debug("germline_snv/before UMD");
  my $vm_umd = get_vm_generic_aa(
				 "-label" => "UMD",
				 "-flag" => "gl-umd-flatfile",
				);
  ram_debug("germline_snv/before UMD");

  #
  #  IARC TP53 database ("the french database"):
  #     http://p53.iarc.fr/DownloadDataset.aspx
  #
  ram_debug("germline_snv/before TP53");
  my $vm_tp53_gl = get_vm_tp53_gl();
  my $vm_tp53_somatic = get_vm_tp53_somatic();
  ram_debug("germline_snv/after TP53");

  #
  #  ASU TERT database:
  #
  #  QUESTIONS:
  # - what to do if variant present in ASU but presentation is blank?
  #   consider "polymorphism"
  # - indels?
  # - does TERT require any special logic outside of ASU matches?
  #
  ram_debug("germline_snv/before ASU TERT");
  printf STDERR "loading ASU TERT...\n";
  my $asu_tert_raw = load_asu_tert();
  my $vm_tert = get_new_vm();
  foreach my $row (@{$asu_tert_raw}) {
    # TERT db provides AA annotations only, no genomic coordinates
    $vm_tert->add_aa(
		     "-aa" => ($row->{"AA substitution"} || die),
		     "-gene" => ($row->{Gene} || die),
		     "-row" => $row
		    );
  }
  ram_debug("germline_snv/after ASU TERT");

  ram_debug("germline_snv/before ALSoD");
  my $vm_alsod = get_vm_mc_db("-name" => "ALSoD");
  ram_debug("germline_snv/after ALSoD");

  my $vm_ikzf1 = get_vm_mc_db("-name" => "IKZF1");

  my $vm_tp53_functional = get_vm_mc_db("-name" => "TP53_functional");

  my $vm_user_variants = get_vm_user_variants("-param" => "user-variants");
  # optional: user variants to medal against
  my $vm_user_blacklist = get_vm_user_variants("-param" => "user-blacklist");
  # optional: user variants that should not receive a medal

  #
  #  process reports:
  #
  my %field_rename = (
		      "GSBClass" => "Somatic_GSBClass",
		      "Reason" => "Somatic_Reason",
		      "Evidence" => "Somatic_Evidence",
		     );

  my $sjpi = get_sjpi();

  my $classification_start_time = time;

  log_msg(sprintf "startup overhead time: %d", $classification_start_time - $start_time);

  foreach my $glf_raw (@INPUT_GL_FILES) {
    my ($glf, $outfile) = get_infile_and_outfile($glf_raw);
    next unless $glf =~ /\w/;
    # in case blank line in input list
    log_msg(sprintf "processing %s...", $glf);
    ram_debug("germline_snv/main loop for $glf");
    line_delimiter_qc("-file" => $glf, "-delimiter" => "\t");
    my $df = new DelimitedFile(
			       "-file" => $glf,
			       "-headers" => 1
			      );

    my $rename_required = grep {$_ eq "GSBClass"} @{$df->headers_raw()};
    unless ($rename_required) {
      my %have = map {$_, 1} @{$df->headers_raw()};
      foreach my $to (values %field_rename) {
	die "ERROR: somatic classification must be run first, can't find either GSBClass or $to in $glf_raw" unless $have{$to};
      }
    }

    my @in_rows;
    my @out_rows;

    while (my $row = $df->get_hash()) {
      push @in_rows, $row;
    }

    if ($exac_vcf2tab and !$ENABLE_SOMATIC_EXAC) {
      printf STDERR "ExAC: starting batch vcf2tab style\n";
      ram_debug("germline_snv/before batch ExAC");
      add_batch_exac(
		     "-rows" => \@in_rows,
		    );
      ram_debug("germline_snv/after batch ExAC");
    }

    unless ($ENABLE_SOMATIC_EXAC) {
      ram_debug("germline_snv/before ExAC coverage");
      add_batch_exac_coverage(
			      "-rows" => \@in_rows,
			     );
      ram_debug("germline_snv/after ExAC coverage");
    }

    if ($FLAGS{"tabix-cosmic"}) {
      ram_debug("germline_snv/before batch COSMIC");
      add_batch_cosmic(
		       "-rows" => \@in_rows,
		       "-germline" => 1
		      );
      ram_debug("germline_snv/after batch COSMIC");
    }

    add_batch_clinvar(
		      "-rows" => \@in_rows,
		     );

    if ($ENABLE_TAYLOR_HOTSPOTS) {
      add_batch_taylor_hotspots("-rows" => \@in_rows);
    }


    foreach my $row (@in_rows) {
      dump_die($row, "start variant lookup", 1) if $FLAGS{"debug-start"};

      init_germline_tests($row);

      if ($rename_required) {
	while (my ($from, $to) = each %field_rename) {
	  $row->{$to} = $row->{$from};
	}
      }

      foreach my $f ($FIELD_DBSNP,
		     $FIELD_NHLBI_FREQ,
		     $FIELD_FOLDX,
		     $NSFP_FIELDS[0]) {
	die "missing field $f, somatic classification is now a prerequisite" unless exists $row->{$f};
      }
      add_pcgp_germline_population($row, $vm_pcgp_germline);

      my $gene = $row->{GeneName} || die "no GeneName field";
#      my $variant_class = lc($row->{Class} || die "no Class field");
      my $variant_class = lc($row->{Class} || "");
      my $aa = $row->{AAChange};

      my $is_rare = get_is_rare($row);
      printf STDERR "is_rare=%d\n", $is_rare if $VERBOSE_GERMLINE_PROGRESS;

      my ($damaging_call_available, $is_damaging, $is_tolerated) = parse_damaging($row);
      die "ERROR: sift/polyphen fields not available" unless defined $damaging_call_available;

      my $medal;

      my @reasons;
      if (my $list = $row->{$INTERNAL_FIELD_QUEUE_REASONS}) {
	# NOTE: these will be applied even to rows with non-gold genes
	push @reasons, @{$list};
	# from batch tabix
	delete $row->{$INTERNAL_FIELD_QUEUE_REASONS};
      }

      my @pmid;
      if (my $pmids = $row->{$INTERNAL_FIELD_QUEUE_EVIDENCE}) {
	# from batch tabix, e.g. ExAC
	push @pmid, @{$pmids};
	delete $row->{$INTERNAL_FIELD_QUEUE_EVIDENCE};
      }

      my $hits;

      #      my $SNP_FIELD = "UniSNP_ids";
      # from dbNSFP: old version (snp129), plus only nonsynonymous coverage

      my $is_novel;
      if (exists $row->{$FIELD_DBSNP}) {
	$is_novel = $row->{$FIELD_DBSNP} ? 0 : 1;
	# "novel" = no dbSNP annotation
      } else {
	die "no SNP field $FIELD_DBSNP";
      }

      #      die "no dbSNP137 field" unless exists $row->{dbSNP137};
      #      printf STDERR "WARNING: no dbSNP137 field!\n" unless exists $row->{dbSNP137};
      # change to a warning: some test data may not have this field

      my $is_novel_these_db = 1;
      # obsolete: whether we see the variants in the various other
      # database lookups

      my $medal_call_protected_from_nhlbi;
      # set if medal call cannot be overridden by NHLBI frequency check
      if ($row->{$INTERNAL_FIELD_NO_POPULATION_FILTER}) {
	$medal_call_protected_from_nhlbi = 1;
	delete $row->{$INTERNAL_FIELD_NO_POPULATION_FILTER};
      }

      my $user_blacklisted;
      my $panel_decision = "";

      my $is_gold_gene = 0;
      if ($GERMLINE_GENE_GSM) {
	$is_gold_gene = 1 if $gsm_gold_genes->contains($gene) or
	  $gsm_gold_genes->resolve("-symbol" => $gene);
      } else {
	$is_gold_gene = $gl_gold_genes->{$gene};
      }

      if (0) {
	# disabled 3/26/2018: this check is problematic as these large
	# intervals may overlap multiple genes, not all of them
	# reviewable; e.g. in GRCh37, a variant at 13.48996929 overlaps RB1
	# (reviewable) but also LPAR6 (not reviewable).  Disable as we
	# are already using localized gene symbols.
	#
	# These intervals are still used by e.g. cloud medal ceremony
	# for VCF pre-fitlering to genes of interest.
	unless ($is_gold_gene) {
	  my $hits = $grf_gold->find(
				     "-chrom" => ($row->{Chr} || die),
				     "-base" => get_sj_pos($row)
				    );

	  if (@{$hits}) {
	    $is_gold_gene = 1;
	    #	foreach my $h (@{$hits}) {
	    #	  dump_die($h, "debug!");
	    #	}
	    #	die;
	    my %hit_genes = map {$_->{gene}, 1} @{$hits};

	    printf STDERR "WARNING: found gold gene genomically at %s.%s but sym lookup for %s failed!! hit=%s\n",
	      $row->{Chr},
		get_sj_pos($row),
		  $gene, join ",", sort keys %hit_genes;
	  }
	}
      }

      my $silver_lock;
      my $committee_call_medal;
      my $is_silent_intron_utr = 0;
      if ($variant_class eq "silent" or
	  $variant_class eq "intron" or
	  $variant_class eq "utr_5" or
	  $variant_class eq "utr_3") {
	$is_silent_intron_utr = 1;
      }

#      if (not($is_gold_gene) and not($in_promoter_region)) {
      if (not($is_gold_gene)) {
	#	      $gl_gold_noncancer_genes->{$gene} or
	#	      $gl_silver_genes->{$gene} or
	#	      $gl_bronze_genes->{$gene}
	push @reasons, "not_gold_gene";
	# might be helpful to highlight discrepancies with
	# SJ gold variant list
      } else {  # gold gene
	#
	#  start matching:
	#

	if (my $mlist = $row->{$INTERNAL_FIELD_QUEUE_MEDALS}) {
	  # from batch tabix; only apply if in a desired gene
	  foreach my $request (@{$mlist}) {
	    $medal = request_medal($medal, $request);
	  }
	  delete $row->{$INTERNAL_FIELD_QUEUE_MEDALS};
	}

	#
	#  reviewed committee calls:
	#
	my $rm;
	($rm, $panel_decision, $committee_call_medal) =
	  committee_check(
			  "-vm" => $vm_committee,
			  "-row" => $row,
			  "-is-silent" => $is_silent_intron_utr,
			  "-protect-flag" => \$medal_call_protected_from_nhlbi,
			  "-reasons" => \@reasons,
			 );
	$medal = request_medal($medal, $rm) if $rm;
	if ($panel_decision) {
	  # variant has been seen before (perfect match)
	  if ($panel_decision eq "P") {
	    push @reasons, generate_acmg_reason(ACMG_2015, "PS1", "SJ_committee");
	  }
	} elsif ($variant_class eq "missense" and $aa) {
	  # no perfect match
	  my $hits = $vm_committee->find_aa_specific_codon(
						    "-gene" => $gene,
						    "-aa" => $aa
						   );
	  my @sj_hits;
	  foreach my $h (@{$hits}) {
	    next unless lc($h->{$FIELD_CLASS}) eq "missense" and
	      $h->{paneldecision} eq "P";
	    # "A novel missense amino acid change occurring at the same
	    # position as another pathogenic missense change (e.g.,
	    # Trp38Ser and Trp38Leu) is considered moderate evidence"
	    #
	    # => want to match to missense specifically, ignoring e.g.
	    #    frameshifts which are inherently pathogenic
	    push @sj_hits, join "_", @{$h}{$FIELD_GENE, $FIELD_AACHANGE};
	  }
	  push @reasons, generate_acmg_reason(ACMG_2015, "PM5", "SJ_committee", @sj_hits) if @sj_hits;
	}

	if (0) {
	  # DISABLED 2/6/2018 per JZ: the only PCGP germline variants
	  # that should get gold are NEJM P/LP.  Since NEJM is a subset
	  # of committee decisions that logic should already cover it.

	  # PCGP germline variants:
	  $medal = assign_aa_match(
				   "-vm" => $vm_pcgp_germline,
				   "-row" => $row,
				   "-is-silent" => $is_silent_intron_utr,
				   "-silent-allowed" => 0,
				   "-protect" => 0,
				   "-transcript-handshake" => $FIELD_REFSEQ,
				   "-current-medal" => $medal,
				   "-reasons" => \@reasons,
				   "-label" => "PCGP_germline",
				   "-novel-flag" => \$is_novel_these_db,
				   "-medal-perfect" => CLASS_GOLD,
				   "-medal-codon" => CLASS_SILVER,
				   "-medal-site" => CLASS_SILVER,
				   "-genomic-lookup" => 1,
				  );
	}

	if ($is_damaging and
	    ($row->{"Somatic_GSBClass"} || "") eq CLASS_GOLD) {
	  # MW "Definitions and TP53a.pptx", 10/28/2013:
	  # "Silver: A variant predicted pathologic mutation which has
	  # been identified as a somatic gold mutation"
	  # Q: should this be restricted to the germline gold list?
	  # A: yes (12/5/2013 committee meeting)
	  # (this code now moved inside the gene check block)
	  $medal = request_medal($medal, CLASS_SILVER);
	  push @reasons, "somatic_gold_predicted_deleterious";
	}

	#
	#  truncating and in-frame indel checks:
	#
	# MW email 10/28/2013: "We should run for all genes."
	my $tclass = "";
	if ($variant_class) {
	  $tclass = $variant_class_truncating{$variant_class};
	  die "truncation class status not defined for $variant_class" unless defined $tclass;
	}

	if ($tclass) {
	  # style as used in somatic classification
	  my $is_ts;
	  if ($GERMLINE_TS_ONCO_DB) {
	    unless ($ts_onco_db->find_row("-gene" => $gene)) {
	      printf STDERR "WARNING: no LoF/GoF annotations for reviewable gene %s\n", $gene;
	    }
	    $is_ts = $ts_onco_db->is_lof("-gene" => $gene);
	  } else {
	    printf STDERR "WARNING: no tumor suppressor status for %s\n", $gene
	      unless exists $ts_genes->{$gene};
	    $is_ts = $ts_genes->{$gene};
	  }
	  if ($is_ts) {
	    #
	    # tumor suppressor/loss-of-function gene: truncations are gold
	    #
	    #	    printf STDERR "truncating TS medal for %s is %s\n", $tclass, ts_truncating_medal($tclass);
	    my $request = ts_truncating_medal($tclass);
	    if ($request eq CLASS_GOLD) {
	      # nonsense/frameshift/splice in loss-of-function gene
	      my @info;

	      my $chr = get_chr($row);
	      my $pos = get_sj_pos($row);
	      my $nm = $row->{$FIELD_REFSEQ} || "";
	      if (my $warn_regions = $acmg_pvs1_caution->{$nm}) {
		my $use_caution;
		foreach my $region (@{$warn_regions}) {
		  my ($w_chr, $w_start, $w_end) = @{$region};
		  if ($w_chr eq $chr and $pos >= $w_start and $pos <= $w_end) {
		    $use_caution = 1;
		  }
		}
		if ($use_caution) {
		  # "One must also be cautious when interpreting truncating
		  # variants downstream of the most 3′ truncating variant
		  # established as pathogenic in the literature. This is
		  # especially true if the predicted stop codon occurs in
		  # the last exon or in the last 50 base pairs of the
		  # penultimate exon..."
		  push @info, "CAUTION_NEAR_TRANSCRIPT_END";
		  # TO DO: maybe also downgrade to silver??
		}
	      } else {
		push @info, "CAUTION_NO_TRANSCRIPT_ANNOTATION";
	      }
	      push @reasons, generate_acmg_reason(ACMG_2015, "PVS1", @info);
	    }
	    $medal = request_medal($medal, $request);
	  } else {
	    $medal = request_medal($medal, CLASS_SILVER);
	  }
	  push @reasons, "truncation";
	} elsif ($variant_class_inframe_indel{$variant_class}) {
	  $medal = request_medal($medal, CLASS_SILVER);
	  # spec: "Indel: Novels (not in HGMD) are Silver"
	  # QUESTION: but activating mutations should be gold???
	  push @reasons, "in-frame indel";
	  $silver_lock = 1 unless $medal and $medal eq CLASS_GOLD;
	  # exception: might already be called gold based on ClinVar match
	  push @reasons, generate_acmg_reason(ACMG_2015, "PM4");
	  # TO DO: any checks needed for repeat regions?:
	  # "By contrast, small in-frame deletions/insertions in
	  # repetitive regions, or regions that are not well conserved
	  # in evolution, are less likely to be pathogenic."
	}

	foreach my $ref (
			 [ $vm_tp53_gl, "IARC_TP53_germline", 1 ],
			 [ $vm_tp53_somatic, "IARC_TP53_somatic" ],
			) {
	  my ($db, $reason, $is_germline) = @{$ref};
	  #
	  #  French TP53 databases (germline and somatic):
	  #  - sometimes we have the AA annotation (preferred)
	  #  - use genomic coordinates if not available
	  #
	  my @extra;
	  if ($is_germline) {
	    @extra = (
		      "-hit-callback" => \&tp53_germline_population_callback,
		      "-hit-callback-perfect" => 1,
		     );
	  }

	  $medal = assign_aa_match(
				   "-vm" => $db,
				   "-row" => $row,
				   "-is-silent" => $is_silent_intron_utr,
				   "-silent-allowed" => 0,
				   "-protect" => 1,
				   "-protect-flag" => \$medal_call_protected_from_nhlbi,
				   "-transcript-handshake" => $FIELD_REFSEQ,
				   "-current-medal" => $medal,
				   "-reasons" => \@reasons,
				   "-label" => $reason,
				   "-novel-flag" => \$is_novel_these_db,
				   "-medal-perfect" => CLASS_GOLD,
				   "-medal-codon" => CLASS_SILVER,
				   "-medal-site" => CLASS_SILVER,
				   "-genomic-lookup" => 1,
				   @extra
				  );
	}

	foreach my $ref (
			 [$vm_lovd_apc, "LOVD_APC"],
			 [$vm_lovd_msh2, "LOVD_MSH2"],
			) {
	  my ($vm, $label) = @{$ref};
	  #
	  #  LOVD:
	  #
	  $medal = assign_aa_match(
				   "-vm" => $vm,
				   "-row" => $row,
				   "-is-silent" => $is_silent_intron_utr,
				   "-silent-allowed" => 0,
				   "-protect" => 0,
				   "-transcript-handshake" => $FIELD_REFSEQ,
				   "-current-medal" => $medal,
				   "-reasons" => \@reasons,
				   "-label" => $label,
				   "-novel-flag" => \$is_novel_these_db,
				   "-medal-perfect" => CLASS_SILVER,
				   "-medal-codon" => CLASS_SILVER,
				   "-medal-site" => CLASS_SILVER,
				   "-hit-callback" => \&lovd_germline_callback,
				   "-hit-callback-perfect" => 1,
				  );
	}



	foreach my $ref (
			 [ $vm_nhgri_brca1, "NHGRI_BRCA1" ],
			 [ $vm_nhgri_brca2, "NHGRI_BRCA2" ],
			) {
	  #
	  #  NHGRI databases:
	  #
	  my ($db, $label) = @{$ref};
	  #	  dump_die($row, "start NHGRI lookup", 1);
	  $medal = assign_aa_match(
				   "-vm" => $db,
				   "-row" => $row,
				   "-is-silent" => $is_silent_intron_utr,
				   "-silent-allowed" => 0,
				   "-protect" => 0,
				   "-transcript-handshake" => $FIELD_REFSEQ,
				   "-reasons" => \@reasons,
				   "-label" => $label,
				   "-novel-flag" => \$is_novel_these_db,
				   "-genomic-lookup" => 1,
				   "-medal-override" => \&nhgri_medal_callback,

				   "-current-medal" => $medal,
				   "-medal-perfect" => CLASS_SILVER,
				   "-medal-codon" => CLASS_SILVER,
				   "-medal-site" => CLASS_SILVER,
				   # medals will be overridden in callback

				   "-hit-callback" => \&nhgri_info_callback,
				   "-hit-callback-perfect" => 1,
				  );
	}

	#
	#  COSMIC oncogenic sites (hotspots):
	#
	$medal = assign_aa_match(
				 "-vm" => $vm_cosmic_hotspots,
				 "-row" => $row,
				 "-is-silent" => $is_silent_intron_utr,
				 "-silent-allowed" => 0,
				 "-protect" => 0,
				 "-transcript-handshake" => "Accession Number",
				 # these annotations are frequently
				 # incompatible (ENST), but are sometimes
				 # available in NM_ format
				 "-current-medal" => $medal,
				 "-reasons" => \@reasons,
				 "-label" => "COSMIC_hotspot",
				 "-novel-flag" => \$is_novel_these_db,
				 "-medal-perfect" => CLASS_GOLD,
				 "-medal-codon" => CLASS_SILVER,
				 "-medal-site" => CLASS_SILVER,
				 "-genomic-lookup" => 1,
				 "-hit-callback" => \&assign_acmg_pm1,
				 "-hit-callback-perfect" => 1,

				);

	#
	#  RB1 database:
	#
	$medal = assign_aa_match(
				 "-vm" => $vm_rb1,
				 "-row" => $row,
				 "-is-silent" => $is_silent_intron_utr,
				 "-silent-allowed" => 0,
				 "-protect" => 0,
				 "-transcript-handshake" => "transcript_id",
				 "-current-medal" => $medal,
				 "-reasons" => \@reasons,
				 "-label" => "RB1",
				 "-novel-flag" => \$is_novel_these_db,
				 "-medal-perfect" => CLASS_SILVER,
				 "-medal-codon" => CLASS_SILVER,
				 "-medal-site" => CLASS_SILVER,
				 "-genomic-lookup" => 1,
				 "-pubmed-list" => \@pmid,
				 "-pubmed-field" => "PubMed",
				);

	#
	# TERT database:
	# contains codon info only, genomic coordinates are not provided.
	#
	my $tert_callback = sub {
	  my ($hits, $medal) = @_;
	  my $medal_orig = $medal;
	  my $has_hq = 0;
	  my $has_lq = 0;
	  foreach my $h (@{$hits}) {
	    # in rare cases multiple rows might be returned
	    # e.g. for a codon-only match to chr5.1278894.C.T / A716A
	    # codon-only
	    my $p = $hits->[0]->{Presentation};
	    if (not($p) or $p eq "polymorphism") {
	      # override non-disease-specific variants
	      $has_lq = 1;
	    } else {
	      $has_hq = 1;
	    }
	  }
	  $medal = CLASS_BRONZE unless $has_hq;

	  printf STDERR "TERT hit: presentation=%s, requested=%s final=%s\n",
	    ($hits->[0]->{Presentation} || ""),
	      $medal_orig, $medal;

	  return $medal;
	};

	$medal = assign_aa_match(
				 "-vm" => $vm_tert,
				 "-row" => $row,
				 "-is-silent" => $is_silent_intron_utr,
				 "-silent-allowed" => 0,
				 "-protect" => 0,
				 "-transcript-handshake" => $FIELD_REFSEQ,
				 "-current-medal" => $medal,
				 "-reasons" => \@reasons,
				 "-label" => "ASU_TERT",
				 "-novel-flag" => \$is_novel_these_db,

				 "-medal-perfect" => CLASS_GOLD,
				 "-medal-codon" => CLASS_SILVER,

				 "-medal-override" => $tert_callback,

				 "-pubmed-list" => \@pmid,
				 "-pubmed-field" => "PMID",
				);

	#
	#  HGMD:
	#  - always has an AA annotation
	#
	$medal = assign_aa_match(
				 "-vm" => $vm_hgmd,
				 "-row" => $row,
				 "-is-silent" => $is_silent_intron_utr,
				 "-silent-allowed" => 0,
				 "-protect" => 0,

				 "-transcript-handshake" => "refcore",
				 "-current-medal" => $medal,
				 "-reasons" => \@reasons,
				 "-label" => "HGMD",
				 "-novel-flag" => \$is_novel_these_db,
				 "-medal-perfect" => CLASS_SILVER,
				 "-medal-codon" => CLASS_SILVER,
				 # MW 11/6/2013:
				 # "By our current definitions you can apply a silver medal
				 # for the codon matches too."
				 "-genomic-lookup" => 1,
				 "-medal-site" => CLASS_SILVER,

				 "-pubmed-list" => \@pmid,
				 "-pubmed-field" => "pmid",
				 #				 "-pubmed-debug" => 1,
				);

	#
	#  ARUP MEN2 (RET):
	#
	# MW 11/6/2013:
	# "For the MEN2 link the pathologic variants we should include as gold"

	my $ret_handshake = $FIELD_REFSEQ;
	# THIS ESSENTIALLY DISABLES DATABASE:
	# SJ USES A DIFFERENT TRANSCRIPT!
	if (0) {
	  printf STDERR "FIX ME: RET HANDSHAKING TEMPORARILY DISABLED!!\n";
	  $ret_handshake = "";
	}

	$medal = assign_aa_match(
				 "-vm" => $vm_arup,
				 "-row" => $row,
				 "-is-silent" => $is_silent_intron_utr,
				 "-silent-allowed" => 0,
				 "-protect" => 0,
				 "-transcript-handshake" => $ret_handshake,
				 "-current-medal" => $medal,
				 "-reasons" => \@reasons,
				 "-label" => "ARUP",
				 "-novel-flag" => \$is_novel_these_db,
				 "-medal-perfect" => CLASS_GOLD,
				 "-medal-codon" => CLASS_SILVER,
				);

	#
	#  umd.be:
	#
	# MW 11/6/2013:
	# "these are called causal variants (BRCA1/BRCA2/APC) for simple SNVs
	# we should call as silver and the fs/truncating we should include as
	# gold"
	# - fs/truncating will be generically called gold elsewhere
	# - other matches will be silver
	$medal = assign_aa_match(
				 "-vm" => $vm_umd,
				 "-row" => $row,
				 "-is-silent" => $is_silent_intron_utr,
				 "-silent-allowed" => 1,
				 # APC p.Leu410Leu
				 "-protect" => 0,
				 "-transcript-handshake" => $FIELD_REFSEQ,
				 "-current-medal" => $medal,
				 "-reasons" => \@reasons,
				 "-label" => "UMD",
				 "-novel-flag" => \$is_novel_these_db,
				 "-medal-perfect" => CLASS_SILVER,
				 # NOT gold for this db
				 "-medal-codon" => CLASS_SILVER,
				 "-hit-callback" => \&umd_info_callback,
				 "-hit-callback-perfect" => 1,
				);


	#
	#  COSMIC somatic variants:
	#
	if ($FLAGS{"tabix-cosmic"}) {
	  # tabix rewrite.  only works for perfect matches,
	  # however previous site-only matches are arguably a bug
	  # since this feature was meant to flag COSMIC truncating
	  # matches only.
	  if ($row->{$FIELD_COSMIC_TRUNCATING_FLAG}) {
	    $medal = CLASS_GOLD;
	    push @reasons, "COSMIC_somatic_truncating";
	  }
	} else {
	  $medal = assign_aa_match(
				   "-vm" => $vm_cosmic,
				   "-row" => $row,
				   "-is-silent" => $is_silent_intron_utr,
				   "-silent-allowed" => 0,
				   "-protect" => 0,
				   "-transcript-handshake" => "Accession Number",
				   "-current-medal" => $medal,
				   "-reasons" => \@reasons,
				   #				 "-label" => "COSMIC_somatic",
				   "-label" => "COSMIC_somatic_truncating",
				   # 3/2014: update to be more specific:
				   # COSMIC lookup done here is NOT universal
				   "-novel-flag" => \$is_novel_these_db,
				   "-medal-perfect" => CLASS_GOLD,
				   "-medal-codon" => CLASS_SILVER,
				   "-medal-site" => CLASS_SILVER,
				   "-genomic-lookup" => 1,
				  );
	}

	#
	#  ExAC (passive for now):
	#
	unless ($exac_vcf2tab) {
	  # old version
	  my $vm_exac = get_vm_exac(
				    "-row" => $row,
				    "-exac" => $tf_exac
				   );
	  # build a VariantMatcher including only variants in region

	  #	dump_die($row, "before exac lookup", 1);
	  my $hits = $vm_exac->find_snv(
					"-sj" => $row
				       );
	  $hits = indel_search($vm_exac, $row) unless $hits;

	  if ($hits) {
	    if (@{$hits} == 1) {
	      # usable single match
	      #	    dump_die($row, "debug before exac", 1);
	      populate_exac(
			    "-hit" => $hits->[0],
			    "-reasons" => \@reasons,
			   );
	    } else {
	      push @reasons, "exac_multiple_indel_matches";
	    }
	  }
	}
	push @reasons, $row->{$FIELD_EXAC_COVERAGE} if $row->{$FIELD_EXAC_COVERAGE};

	#
	# ExAc end
	#

	#
	#  private/internal/unpublished variants:
	#
	$medal = assign_aa_match(
				 "-vm" => $vm_private,
				 "-row" => $row,
				 "-is-silent" => $is_silent_intron_utr,

				 "-silent-allowed" => 1,
				 "-protect" => 1,
				 "-protect-flag" => \$medal_call_protected_from_nhlbi,
				 # chr12   12037444        C       A
				 # this is both common in NHLBI AND silent!!

				 "-transcript-handshake" => $FIELD_REFSEQ,
				 "-current-medal" => $medal,
				 "-reasons" => \@reasons,
				 "-label" => "SJ_private",
				 "-novel-flag" => \$is_novel_these_db,
				 "-medal-perfect" => CLASS_SILVER,
				 "-medal-codon" => CLASS_BRONZE,
				 "-medal-site" => CLASS_BRONZE,
				 "-genomic-lookup" => 1,
				);

	#
	#  ALSoD:
	#
	$medal = assign_aa_match(
				 "-vm" => $vm_alsod,
				 "-row" => $row,
				 "-is-silent" => $is_silent_intron_utr,
				 "-silent-allowed" => 0,
				 "-protect" => 0,
#				 "-transcript-handshake" => $FIELD_REFSEQ,
				 "-current-medal" => $medal,
				 "-reasons" => \@reasons,
				 "-label" => "ALSoD",
				 "-novel-flag" => \$is_novel_these_db,
				 "-medal-perfect" => CLASS_SILVER,
				 "-medal-codon" => CLASS_SILVER,
				 "-medal-site" => CLASS_BRONZE,
				 "-genomic-lookup" => 1,
				);

	#
	#  IKZF1:
	#
	$medal = assign_aa_match(
				 "-vm" => $vm_ikzf1,
				 "-row" => $row,
				 "-is-silent" => $is_silent_intron_utr,
				 "-silent-allowed" => 0,
				 "-protect" => 0,
#				 "-transcript-handshake" => $FIELD_REFSEQ,
				 "-current-medal" => $medal,
				 "-reasons" => \@reasons,
				 "-label" => "IKZF1_deleterious",
				 "-novel-flag" => \$is_novel_these_db,
				 "-medal-perfect" => CLASS_SILVER,
				 "-medal-codon" => CLASS_SILVER,
				 "-medal-site" => CLASS_BRONZE,
				 "-genomic-lookup" => 1,
				 "-pubmed-list" => \@pmid,
				 "-pubmed-field" => F_VDB_PMID()
				);

	#
	#  TP53 functional data:
	#
	if (my $hits = $vm_tp53_functional->find_snv("-sj" => $row)) {
	  push @reasons, sprintf "TP53_functional_data=%s", join ",", @{$hits->[0]}{qw(
											WAF1_Act
											MDM2_Act
											BAX_Act
											_14_3_3_s_Act
											AIP_Act
											GADD45_Act
											NOXA_Act
											p53R2_Act

											Giac_A549_WT_Nut_norm
											Giac_A549_Null_Nut_norm
											Giac_A549_Null_Eto_norm

											Kolt_RFS_H1299_norm
										     )
										  };

	  push @reasons, sprintf "TP53_functional_call=%s", $hits->[0]->{functional_call};
	}

	$medal = assign_aa_match(
				 "-vm" => $vm_user_variants,
				 "-row" => $row,
				 "-is-silent" => $is_silent_intron_utr,
				 "-silent-allowed" => 0,
				 "-protect" => 0,
#				 "-transcript-handshake" => $FIELD_REFSEQ,
				 "-current-medal" => $medal,
				 "-reasons" => \@reasons,
				 "-label" => "user_variant",
				 "-novel-flag" => \$is_novel_these_db,
				 "-medal-perfect" => CLASS_GOLD,
				 "-medal-codon" => CLASS_SILVER,
				 "-medal-site" => CLASS_SILVER,
				 "-genomic-lookup" => 1,
#				 "-pubmed-list" => \@pmid,
#				 "-pubmed-field" => F_VDB_PMID()
				);

	if ($vm_user_blacklist->find_snv("-sj" => $row) or
	    $vm_user_blacklist->find_indel(
					   "-sj" => $row,
					   "-match-basic-type" => 1,
					   "-match-size" => 1,
					   "-fuzz-bases" => 0
					  )) {
	  $user_blacklisted = 1;
	}


	if ($is_rare and not($is_silent_intron_utr)) {
	  #
	  #  deleterious predictions: call for rare non-silent variants only
	  #
	  if ($is_novel) {
	    if ($is_damaging) {
	      # Silver:
	      # "A novel variant not seen in > 3 NHLBI and predicted by one or more platforms to be damaging or deleterious"
	      $medal = request_medal($medal, CLASS_SILVER);
	      push @reasons, "rare_novel_predicted_deleterious";
	    } elsif ($is_tolerated) {
	      $medal = request_medal($medal, CLASS_BRONZE);
	      push @reasons, "rare_novel_predicted_tolerated";
	      # MW 11/22/2013:
	      # Bronze if > platforms benign . Bobby's structural platform also has to be incorporated
	      # => how many platforms??
	    }
	  } else {
	    if ($is_damaging) {
	      $medal = request_medal($medal, CLASS_SILVER);
	      push @reasons, "rare_known_predicted_deleterious";
	      # MW 11/22/2013
	    } elsif ($is_tolerated) {
	      # Bronze:
	      # "Rare variant tolerated by prediction algorithms"
	      # QUESTION: novel check too??
	      $medal = request_medal($medal, CLASS_BRONZE);
	      push @reasons, "rare_known_predicted_tolerated";
	    }
	  }
	}
      }

      promoter_check(
		     "-row" => $row,
		     "-medal-ref" => \$medal,
		     "-reasons-ref" => \@reasons,
		    );
      # independent from reviewable gene logic

      $medal = CLASS_SILVER if $silver_lock;
      # for in-frame indels, prevent gold assignment based on
      # unfortunate cases like off-by-1 matches to TP53_somatic, e.g.
      # TP53    In35gene        SJHGG048A       chr17   7578457 proteinDel      R156_R158>R     81      109     13      152     CGGACG  -       GGCCATGGCG(CGGACG>------)CGGGTGCCGG                                     718                            Gold     event=in-frame indel;gene_medal=gold;gene_class=TS,recurrent/oncogene  COSMIC=21901162;1;1,19739123;2;0,11044641;1;0,7852189;1;0,10690522;1;0,0,18772397;1;1,7882357;1;0,21901162;1;1,10690522;1;0,1324794;1;0,14688025;1;0,9816045;1;01235     Gold    in-frame indel;IARC_TP53_germline_site_only;IARC_TP53_somatic;COSMIC_hotspot_site_only;HGMD_site_only;needs_committee_review    PUBMED=10486318,20127978

      unless ($is_rare) {
#	unshift @reasons, "NHLBI_common";
	unshift @reasons, "ExAC_common";
	$medal = CLASS_UNKNOWN unless $medal_call_protected_from_nhlbi;
	# variants that are common in NHLBI should generally not
	# receive medals, unless they have received a medal based on
	# a match to a protected/gold database (e.g. TP53)
      }

      if ($user_blacklisted) {
	$medal = CLASS_UNKNOWN;
	push @reasons, "user_blacklisted";
      }


      if ($GERMLINE_GSB_LOCKED_TO_COMMITTEE and $committee_call_medal) {
	# if we have a GSB call based on a committee call, this should
	# always override auto classification.
	$medal = $committee_call_medal;
      }

      if ($vm_bad->find_snv("-sj" => $row) or
	  $vm_bad->find_aa_substitution("-sj" => $row)) {
	# manually curated list of known bad variants, e.g.
	# TP53 P72R.  These always override medal logic.
	$medal = CLASS_UNKNOWN;
	unshift @reasons, "variant_blacklist";
      }

      $medal = CLASS_UNKNOWN unless $medal;

      push @reasons, "needs_committee_review"
	if $medal ne CLASS_UNKNOWN and not($panel_decision);

      my ($freq, $pop_db_name) = get_is_rare($row, 1);
      if ($freq <= $GERMLINE_ACMG_PM2_MAX_POPULATION_FREQ) {
	# PM2 Absent from controls (or at extremely low frequency if
	# recessive) (Table 6) in Exome Sequencing Project,
	# 1000 Genomes Project, or Exome Aggregation Consortium
	push @reasons, generate_acmg_reason(ACMG_2015, "PM2", $pop_db_name, $GERMLINE_ACMG_PM2_MAX_POPULATION_FREQ);
      } elsif ($freq > $GERMLINE_ACMG_BA1_MAX_POPULATION_FREQ) {
	push @reasons, generate_acmg_reason(ACMG_2015, "BA1", $pop_db_name);
	# BA1 Allele frequency is >5% in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium
      }

      $row->{GSBClass} = $medal;
      $row->{sort_score} = $MEDAL_SORT_WEIGHT{$medal} || die "no weight for $medal";
      $row->{Reason} = join_reason_list(\@reasons);
      $row->{$FIELD_PANEL_DECISION} = $panel_decision;

      if (my $pm = $row->{$FIELD_PUBMED_IDS}) {
	# copy PubMed IDs from somatic
	push @pmid, split /,/, $pm;
      }

      my @evidence;

      @pmid = grep {/^\d+$/} @pmid;
      # exclude bogus values e.g. LSDB, NO ID, etc.
      if (@pmid) {
	#	push @evidence, sprintf 'PUBMED=%s', join ",", @{unique_ordered_list(\@pmid)};
	push @evidence, sprintf 'PUBMED=%s', join ",", sort {$a <=> $b} @{unique_list(\@pmid)};
      }
      #      $row->{Evidence} = join ":", @evidence;
      $row->{Evidence} = join_evidence_field(\@evidence);
      # use ":" as delimiter for better compatibility with SITH
      push @out_rows, $row;
    }
    ram_debug("germline_snv/after row loop");

    #
    #  write output:
    #
    my %have_headers = map {$_, 1} @{$df->{headers_raw}};
    die unless %have_headers;

    my @extra;
    push @extra, (
		  $FIELD_PANEL_DECISION,
		  $FIELD_PCGP_COUNT_GERMLINE,
		  $FIELD_PCGP_DENOMINATOR_GERMLINE,
		  "GSBClass",
		  "Reason",
		  "Evidence",
		 );

    my $rpt = $df->get_reporter(
				"-file" => $outfile,
				"-rename" => \%field_rename,
				"-extra" => \@extra
			       );
    $rpt->auto_qc(1);
    $rpt->write_headers();
    # 2/2014: force writing header line even if file contains
    # no data rows (per Matt Parker)

    my @sorted;
    if ($FLAGS{"no-sort"}) {
      @sorted = @out_rows;
    } else {
      @sorted = sort {$b->{sort_score} <=> $a->{sort_score}} @out_rows;
    }

    foreach my $row (@sorted) {
      $rpt->end_row($row);
    }
    $rpt->finish();
  }


  my $classification_end_time = time;
  log_msg(sprintf "germline reports done, startup overhead:%d classification time:%d",
	  ($classification_start_time - $start_time),
	  ($classification_end_time - $classification_start_time)
	 );
  ram_debug("germline_snv/done");

}

sub load_asu_tert {
  #
  #  parse ASU data for TERT:
  #  http://telomerase.asu.edu/diseases.html
  #
  my $asu_tert = $FLAGS{"asu-tert"} || die "-asu-tert";
  my $lines = read_simple_file($asu_tert);
  my %nm;
  foreach my $line (@{$lines}) {
    if ($line =~ /(NM_\d+)/) {
      $nm{$1} = 1;
    }
  }
  die unless scalar keys %nm == 1;
  my ($nm) = keys %nm;

  my $tables = parse_html_tables("-file" => $asu_tert, "-keep-html" => 1);
  die unless @{$tables} == 1;
  my $rows = $tables->[0];

  my $GENE = "TERT";
  # hack

  my $last_domain;
  foreach my $row (@{$rows}) {
    $row->{Gene} = $GENE;
    $row->{$FIELD_REFSEQ} = $nm;
    if (my $d = $row->{Domains}) {
      $last_domain = $d;
    } else {
      $row->{Domains} = $last_domain;
    }

    my $aa = $row->{"AA substitution"};
    if ($aa and $aa =~ /</) {
      # hack: strip HTML from AA (still needed for References to parse info)
      $aa =~ s/<[^>]+>//g;
      $row->{"AA substitution"} = $aa;
    }
    my @pmid = $row->{References} =~ /pubmed\/(\d+)/g;
    $row->{PMID} = join ",", @pmid;
  }
  return $rows;
}


sub parse_iarc {
  my (%options) = @_;
  my $gene = "TP53";
  my $df = new DelimitedFile(
			     "-file" => ($options{"-file"} || die "-file"),
			     "-headers" => 1,
			    );
#  my $vm = get_new_vm(1);
  my $vm = get_new_vm();
  my $blank_genomic = 0;

  my @f_trimmed = (
		   $FIELD_REFSEQ,
		   "Generation",
		   "Individual_ID",
		  );

  my @f_hgvs = (
		 "g_description_hg19",
		 # R16/R17

		 "g_description"
		 # R18
		 );

  my $f_aa = "ProtDescription";
  my $f_hgvs;

  while (my $row = $df->get_hash()) {
    $row->{$FIELD_REFSEQ} = $IARC_TP53_NM;

    my $save_row;
    if ($TRIM_TP53_INFO) {
      my %save;
      @save{@f_trimmed} = @{$row}{@f_trimmed};
      $save_row = \%save;
    } else {
      $save_row = $row;
    }

    #
    #  amino acid details:
    #
    dump_die($row, "where is $f_aa") unless exists $row->{$f_aa};
    my $aa_raw = $row->{$f_aa};
    #    printf STDERR "IARC AA %s\n", $aa_raw;
    my @aa;
    if ($aa_raw =~ /^\[p\.(.*)[\;\/](.*)\]$/) {
      # [p.Q167Y;H168L] => Q167Y H168L
      # [p.H168L/M169I] => H168L M169I
      @aa = ($1, $2);
#      printf STDERR "TP53: double AA for %s: %s %s\n", $genomic, @aa;
    } else {
      @aa = $aa_raw;
    }

    foreach my $aa (@aa) {
      # known formatting, either substitution or codon-based
      $vm->add_aa(
		  "-aa" => $aa,
		  "-gene" => $gene,
		  "-row" => $save_row
		 );
    }

    #      printf STDERR "can't parse $aa\n";
    #      my $genomic = $row->{g_description_hg19} || dump_die($row, "WTF: no g_description_19");

    if (exists $row->{sj_pos}) {
      # refactor/liftover for hg38+
      if (my $pos = $row->{sj_pos}) {
	my $ra = $row->{sj_ref_allele} || dump_die($row, "no ref allele");
	my $va = $row->{sj_alt_allele} || dump_die($row, "no variant allele");
	my $v = new Variant();
	$v->import_generic(
		     "-reference-name" => ($row->{Chr} || die),
		     "-base-number" => $pos,
		     "-reference-allele" => $ra,
		     "-variant-allele" => $va,
		    );
	$vm->add_variant($v, "-row" => $row);
      }

      next;
      # skip old parser
    }

    #
    #  genomic variant details:
    #
    unless ($f_hgvs) {
      # detect field
      foreach my $f (@f_hgvs) {
	if (exists $row->{$f}) {
	  $f_hgvs = $f;
	  last;
	}
      }
      die "can't find hgvs field" unless $f_hgvs;
    }
    my $genomic = $row->{$f_hgvs};
    unless ($genomic) {
      printf STDERR "WARNING: no IARC hg19 entry, skipping\n";
      next;
    }

    if (0) {
      printf STDERR "IARC DEBUG\n";
      $genomic = "NC_000017.10:g.7579472C>T";
    }

    my ($reference, $stuff);
    if ($genomic =~ /:/) {
      # R16: NC_000017.10:g.7577545T>C
      ($reference, $stuff) = split /:/, $genomic;
    } else {
      # R18: g.7578502A
      $stuff = $genomic;
      $reference = 17;
      # entire database is for single gene
    }
    $stuff =~ s/\([^\(]+\)$//;
    # strip interfering comment:
    # NC_000017.10:g.7577157_7577498del342(del intron7)

    if ($stuff =~ /^g\.(\d+)([A-Z])\>([A-Z])$/) {
      # SNV: e.g. NC_000017.10:g.7579601G>C
      #      printf STDERR "import %s.%s.%s.%s\n", $reference, $1, $2, $3;
      $vm->add_snv(
		   "-row" => $save_row,
		   "-reference" => $reference,
		   "-base-number" => $1,
		   "-reference-base" => $2,
		   "-variant-base" => $3
		  );
      # TO DO: sanity-check reference base
    } elsif ($stuff =~ /^g\.(\d+)_(\d+)([A-Z]+)\>([A-Z]+)$/) {
      # dinucleotide, e.g. NC_000017.10:g.7578533_7578534TC>AG
      my ($start, $end, $ref, $var) = ($1, $2, $3, $4);
      my $len = length $ref;
      die unless length($ref) == length($var);
      if ($end == $start + $len - 1) {
	$vm->add_snv(
		     "-row" => $save_row,
		     "-reference" => $reference,
		     "-base-number" => $start,
		     "-reference-base" => $ref,
		     "-variant-base" => $var
		    );
      } else {
	cluck "WTF: dinucleotide consistency epic fail for $stuff";
      }
    } elsif ($stuff =~ /^g\.(\d+)del1$/) {
      # single-base deletion, e.g.
      # NC_000017.10:g.7576852del1
      my $base = $1;
      $vm->add_deletion(
			"-row" => $save_row,
			"-reference" => $reference,
			"-start" => $base,
			"-end" => $base,
		       );
    } elsif ($stuff =~ /^g\.(\d+)_(\d+)del\d+$/) {
      # deleted range, e.g.
      # NC_000017.10:g.7578294_7578304del11
      my ($start, $end) = ($1, $2);
      $vm->add_deletion(
			"-row" => $save_row,
			"-reference" => $reference,
			"-start" => $start,
			"-end" => $end,
		       );
      #      } elsif ($stuff =~ /^g\.(\d+)_(\d+)ins(\d+)$/) {
    } elsif ($stuff =~ /^g\.(\d+)_(\d+)ins(\w+)$/) {
      # insertion, e.g.
      # NC_000017.10:g.7578476_7578477ins1
      # NC_000017.10:g.7577089_7577090insCGCC
      # - db seems to mention both reference bases (before + after insertion)
      # - for now, just track both a-la deletion
      my ($start, $end, $thing) = ($1, $2, $3);
      my $insert_count;
      if ($thing =~ /^\d+$/) {
	$insert_count = $thing;
      } else {
	$insert_count = length($thing);
      }
      $vm->add_insertion(
			 "-row" => $save_row,
			 "-reference" => $reference,
			 "-start" => $start,
			 "-end" => $end,
			 "-count" => $insert_count
			);
    } elsif ($stuff =~ /^g\.(\d+)_(\d+)delins/) {
      # complex: deleted region, plus insertion
      # NC_000017.10:g.7577106_7577113delins5
      # - HACK: just track deletion portion for now
      my ($start, $end) = ($1, $2);
      $vm->add_deletion(
			"-row" => $save_row,
			"-reference" => $reference,
			"-start" => $start,
			"-end" => $end,
		       );
    } elsif ($stuff =~ /^g\.(\d+)delins/) {
      # complex: deleted base, plus insertion
      # NC_000017.10:g.7578230delins2
      # - HACK: just track deletion portion for now
      my ($base) = $1;
      $vm->add_deletion(
			"-row" => $save_row,
			"-reference" => $reference,
			"-start" => $base,
			"-end" => $base,
		       );
    } elsif ($genomic =~ /:g\.\?$/) {
      $blank_genomic++;
    } else {
      printf STDERR "can't parse $genomic, SKIPPING\n" unless $genomic =~ /\?/;
    }
  }				# $row

  printf STDERR "blank genomic entries: %d\n", $blank_genomic;

  return $vm;
}

sub parse_hgmd {
  # OBSOLETE

  # shout-out to
  # /nfs_exports/genomes/1/PCGP/BucketIntermediate/germlineSNVs/clinical_sequencing/get_matchingSNPs.pl
  my (%options) = @_;
  my $file = $options{"-file"} || die;

  my $df = new DelimitedFile(
			     "-file" => $file,
			     "-headers" => 1,
			    );

  my $vm = get_new_vm();
  #  my $vm = get_new_vm(1);

  my $nsp = new NucleotideSubstitutionParser();
  while (my $row = $df->get_hash()) {
    if ($vm->add_aa(
		    "-aa" => $row->{descr},
		    "-gene" => $row->{gene},
		    "-row" => $row
		   )) {
      #      printf STDERR "OK parsing %s in %s\n", @row{qw(aachange gene)};
    } else {
      printf STDERR "ERROR: can't parse HGMD %s in %s\n", @{$row}{qw(descr gene)};
    }

    if (0) {
      # DISABLED for now: need update from Erin w/strand info,
      # these are in transcript space and sometimes need to be RC'd

      #
      #  nucleotide coordinates:
      #
      my ($chrom, $start, $end, $change) = @{$row}{qw(chromosome startCoord endCoord hgvs)};
      die "not SNV" unless $start == $end;
      $nsp->parse($change) || die "can't parse $change";
      if ($nsp->is_substitution()) {
	my $ref_base = $nsp->reference_sequence();
	my $var_base = $nsp->variant_sequence();
	my $ok = $vm->add_snv(
			      "-row" => $row,
			      "-reference" => $chrom,
			      "-base-number" => $start,
			      "-reference-base" => $ref_base,
			      "-variant-base" => $var_base
			     );
	#	printf STDERR "HGMD import: %s\n", join " ", $chrom, $start, $ref_base, $var_base, $ok;
      } else {
	die "unhandled match type";
      }
    }

  }

  return $vm;
}

sub get_new_vm {
  my ($force_rsc) = @_;
  my $rsc = $force_rsc || $REFERENCE_SANITY_CHECK;
  my @options;
  push @options, ("-fasta_dir" => $FLAGS{"fasta-dir"} || die "-fasta-dir") if $rsc;
  return new VariantMatcher(@options);
}

sub dump_hgmd {
  #
  # 12/2013: use new view provided by Erin Hedlund
  #
  my $dbi = get_dbi_gedi("-type" => "research");

  #  my $cmd = 'select * from hgmd.classification_dm_mp';
  #  my $cmd = 'select * from hgmd.classification_dm_mp order by chromosome';
  my $cmd = sprintf 'select * from hgmd.%s order by chromosome', $FLAGS{"hgmd-table"} || "classification_dm_mp";

  # ordering by chrom optimizes reference sequence verification check
  #  my $outfile = "export_hgmd_classification_dm_mp.tab";
  #  my $outfile = "export_hgmd_classification_dm_mp_with_indels.tab";
  #  my $outfile = "export_hgmd_classification_dm_mp_with_indels_and_splices.tab";
  my $outfile = "export_hgmd_classification_dm_mp_everything.tab";
  # 3/4/14: additional variant codes added, including some possibly
  # unusable codes

  printf STDERR "SQL: %s\n", $cmd;
  printf STDERR "exporting %s to %s...\n", $cmd, $outfile;

  export_query_to_flatfile(
			   "-dbi" => $dbi,
			   "-sql" => $cmd,
			   "-outfile" => $outfile,
			   "-trim-whitespace" => 1,
			  );
  # -trim-whitespace: fixes issues with tab chars in a few rows, e.g.
  # Fc fragment of IgG, high affinity Ia, receptor (CD64) \^I^I135911^I
  # (quoted tab at end of description field)

  $dbi->disconnect();

  my %counts;
  open(QC, $outfile) || die;
  while (<QC>) {
    chomp;
    my @f = split /\t/, $_;
    $counts{scalar @f}++;
  }
  close QC;
  die "QC fail: column count mismatch!" unless scalar keys %counts == 1;
}

sub dump_hgmd_old {
  my $dbi = get_dbi_hgmd() || die;
  #  my $cmd = "select chromosome, startCoord, endCoord, gene,chrom, genename, gdbid,omimid, amino, deletion, insertion, codon, descr, hgvs, hgvsAll, dbsnp, chromosome from allmut where tag='DM' and base in ('M','P') and StartCoord IS NOT NULL and E.ndCoord IS NOT NULL";
  my $cmd = "select chromosome, startCoord, endCoord, gene,chrom, genename, gdbid,omimid, amino, deletion, insertion, codon, descr, hgvs, hgvsAll, dbsnp, pmid from allmut where tag='DM' and base in ('M','P') and StartCoord IS NOT NULL and EndCoord IS NOT NULL";
  # add: pmid
  # remove: duplicate chromosome field
  printf STDERR "exporting %s...\n", $cmd;
  my $outfile = "export_hgmd_dm.tab";
  export_query_to_flatfile(
			   "-dbi" => $dbi,
			   "-sql" => $cmd,
			   "-outfile" => $outfile,
			   "-trim-whitespace" => 1,
			  );
  # -trim-whitespace: fixes issues with tab chars in a few rows, e.g.
  # Fc fragment of IgG, high affinity Ia, receptor (CD64) \^I^I135911^I
  # (quoted tab at end of description field)

  $dbi->disconnect();

  my %counts;
  open(QC, $outfile) || die;
  while (<QC>) {
    chomp;
    my @f = split /\t/, $_;
    $counts{scalar @f}++;
  }
  close QC;
  die "QC fail: column count mismatch!" unless scalar keys %counts == 1;
}

sub indel_search {
  my ($vm, $sj_row) = @_;
  my @options = (
		 "-sj" => $sj_row,
		 "-match-basic-type" => $INDEL_MATCH_BASIC_TYPE,
		 "-match-size" => $INDEL_MATCH_SIZE,
		 "-fuzz-bases" => INDEL_NUCLEOTIDE_NT_FUZZY_MATCH,
		);
  return $vm->find_indel(@options);
}

sub join_reason_list {
  my ($reasons) = @_;
  return join ";", @{$reasons};
  # values in some key/value pairs may themselves be lists,
  # so use semicolons for outer list.
}

sub join_evidence_field {
  my ($evidence) = @_;
  return join ":", @{$evidence};
}

sub init_germline_tests {
  my ($row) = @_;
  if (my $test = $FLAGS{"gl-test"}) {
    if ($test == 1) {
      # ASU: in a disease
      $row->{GeneName} = "TERT";
      $row->{AAChange} = "P33S";
    } elsif ($test == 2) {
      # ASU: not in a disease
      $row->{GeneName} = "TERT";
      $row->{AAChange} = "Ala279Thr";
    } elsif ($test == 3) {
      # ASU: indel
      $row->{GeneName} = "TERT";
      $row->{AAChange} = "P112del";
    } elsif ($test == 4) {
      # French TP53 GL: SNV
      $row->{Chr} = "chr17";
      $row->{WU_HG19_Pos} = "7579601";
      $row->{ReferenceAllele} = "G";
      $row->{MutantAllele} = "C";
    } elsif ($test == 5) {
      # French TP53 GL: substitution
      $row->{GeneName} = "TP53";
      #	  $row->{AAChange} = "G105C";
      #      $row->{AAChange} = "Gly105Cys";

      #      $row->{AAChange} = "P151T";
      # still in curated set, perfect match

      $row->{AAChange} = "P151A";
      # fuzzy match: still in curated set, different variant codon

      #      $row->{AAChange} = "P72R";
      # REMOVED from curated set

    } elsif ($test == 6) {
      # French TP53 GL: deletion
      $row->{Chr} = "chr17";
      $row->{WU_HG19_Pos} = "7576852";
      # this is the actual site

      #	  $row->{WU_HG19_Pos} = "7576850";
      # fuzzy location tests

      $row->{ReferenceAllele} = "G"; # bogus
      $row->{MutantAllele} = "-";
    } elsif ($test == 7) {
      # French TP53 GL: insertion
      $row->{Chr} = "chr17";
      $row->{WU_HG19_Pos} = 7578476; # actual
      #	  $row->{WU_HG19_Pos} = 7578470;
      $row->{ReferenceAllele} = "-";
      $row->{MutantAllele} = "C"; # bogus
    } elsif ($test == 8) {
      # French TP53 somatic
      $row->{GeneName} = "TP53";
      $row->{AAChange} = "V143A";
    } elsif ($test == 9) {
      # HGMD
      $row->{GeneName} = "AAAS";
      $row->{AAChange} = "Ile482Ser";
    } elsif ($test == 10) {
      # ClinVar: deletion test (this site has both!)
      # CV: 6 31323501 CT C,CCT
      $row->{Chr} = "chr6";
      $row->{WU_HG19_Pos} = 31323502;
      $row->{ReferenceAllele} = "T";
      $row->{MutantAllele} = "-";
      $row->{ClinicalVarNum} = 4;
    } elsif ($test == 11) {
      # ClinVar: insertion test (this site has both!)
      # CV: 6 31323501 CT C,CCT
      $row->{Chr} = "chr6";
      $row->{WU_HG19_Pos} = 31323502;
      $row->{ReferenceAllele} = "-";
      $row->{MutantAllele} = "C";
    } elsif ($test == 12) {
      # ClinVar: deletion
      # CV: 1 94517224 AAG A
      $row->{Chr} = "chr1";
      $row->{WU_HG19_Pos} = 94517225;
      $row->{ReferenceAllele} = "AG";
      $row->{MutantAllele} = "-";
    } elsif ($test == 13) {
      # ClinVar: deletion
      # CV: 1 94517224 AAG A
      # this is an insertion in the same vicininty,
      # test with $INDEL_MATCH_BASIC_TYPE set to 0
      $row->{Chr} = "chr1";
      $row->{WU_HG19_Pos} = 94517225;
      $row->{ReferenceAllele} = "-";
      $row->{MutantAllele} = "CC";
      # indel near same site as ClinVar deletion
    } elsif ($test == 14) {
      # ClinVar: deletion
      # CV: 1 94517224 AAG A
      # outside of fuzzy match range
      $row->{Chr} = "chr1";
      $row->{WU_HG19_Pos} = 94517225 + INDEL_NUCLEOTIDE_NT_FUZZY_MATCH + 5;
      $row->{ReferenceAllele} = "AG";
      $row->{MutantAllele} = "-";
    } elsif ($test == 15) {
      # HGMD stop codon format vs SJ stop codon format
      $row->{GeneName} = "AAAS";
      $row->{AAChange} = "Gln145*";
      # Gln145Term
      #	  $row->{GeneName} = "ABCA4";
      #	  $row->{AAChange} = "Gln1199*";
      # Gln1199Term
    } elsif ($test == 16) {
      # TP53 P72R chr17.7579472.G.C
      # in ClinVar, but w/CLNSIG=2 (insignificant)
      $row->{GeneName} = "TP53";
      $row->{Chr} = "chr17";
      $row->{WU_HG19_Pos} = 7579472;
      $row->{ReferenceAllele} = "G";
      $row->{MutantAllele} = "C";
    } elsif ($test == 17) {
      # known SJ gold variant (perfect match)
      # (Dr. Downing: "gold is gold" / "simple table lookup")
      $row->{GeneName} = "ALK";
      $row->{AAChange} = "R1275Q";
    } elsif ($test == 18) {
      # truncating event in SJ gene
      $row->{GeneName} = "TP53";
      $row->{Class} = "frameshift";
    } elsif ($test == 19) {
      # in-frame insertion in SJ gene
      $row->{GeneName} = "TP53";
      $row->{Class} = "proteinins";
    } elsif ($test == 20) {
      # known SJ gold variant (match to codon number only)
      # per Dr. Downing 10/24/2013
      $row->{GeneName} = "ALK";
      $row->{AAChange} = "R1275A";
    } elsif ($test == 21) {
      # PCGP
      $row->{GeneName} = "TP53";
      $row->{AAChange} = "R248Q";
    } else {
      die "unknown gl test";
    }
  }
}

sub request_medal {
  my ($old_medal, $new_medal) = @_;
  confess "new medal not defined!" unless defined $new_medal;
  my $ok;
  my $old_weight = $old_medal ? $MEDAL_SORT_WEIGHT{$old_medal} : 0;
  my $new_weight = $MEDAL_SORT_WEIGHT{$new_medal} || confess "no weight for $new_medal";
  $ok = $new_weight >= $old_weight ? 1 : 0;
  #  printf STDERR "medal request: old:%s new:%s accept:%s\n", $old_medal, $new_medal, $ok if $old_medal;
  return $ok ? $new_medal : $old_medal;
}

sub parse_iarc_curated {
  my (%options) = @_;
  my $gene = "TP53";
  my $vm = get_new_vm();
  open(IN, $options{"-file"} || die "-file") || die;
  while (my $aa = <IN>) {
    chomp $aa;
    $vm->add_aa(
		"-aa" => $aa,
		"-gene" => $gene,
		"-row" => $aa
	       );

  }
  return $vm;
}

sub assign_aa_match {
  my (%options) = @_;
  my $db = $options{"-vm"} || confess "-vm";
  my $row = $options{"-row"} || confess "-row";

  my $gene = $row->{GeneName} || die "no GeneName field";
#  my $aa = $row->{AAChange} || die "no AAChange field";
  my $aa = $row->{AAChange} || "";
  # may not be present, e.g. promoter
  #  my $gene = $options{"-gene"} || die "-gene";
  #  my $aa = $options{"-aa"} || die "-aa";

  my $no_medal = $options{"-no-medal"};
  my $reasons = $options{"-reasons"} || die;
  my $label_perfect = $options{"-label"} || die;
  my $debug_hits = $options{"-debug-hits"};

  my $medal = $options{"-current-medal"};
  confess "-current-medal" unless $no_medal or exists $options{"-current-medal"};

  my $medal_perfect = $options{"-medal-perfect"};
  die "-medal-perfect" unless $no_medal or $medal_perfect;

  my $medal_codon = $options{"-medal-codon"};
  die "-medal-codon" unless $no_medal or $medal_codon;

  my $medal_site;
  my $genomic_lookup = $options{"-genomic-lookup"};
  if ($genomic_lookup) {
    $medal_site = $options{"-medal-site"};
    confess "need -medal-site" unless $no_medal or $medal_site;
  }

  my $hit_callback = $options{"-hit-callback"};
  my $hit_callback_perfect = $options{"-hit-callback-perfect"};
  my $medal_override = $options{"-medal-override"};

  my $protect = $options{"-protect"};
  die "specify -protect" unless defined $protect;

  my $is_silent = $options{"-is-silent"};
  confess "specify -is-silent" unless defined $is_silent;

  my $silent_allowed = $options{"-silent-allowed"};
  confess "specify -silent-allowed" unless defined $silent_allowed;

  my $silent_and_suppressed = ($is_silent and not($silent_allowed));
  # user does not want to report matches to silent variants,
  # even if they perfectly match the database
  printf STDERR "silent and suppressed for %s\n", join " ", @{$row}{qw(GeneName Class AAChange)} if $silent_and_suppressed and $VERBOSE_SILENT_AND_SUPPRESSED;

  my $protect_flag;
  if ($protect) {
    $protect_flag = $options{"-protect-flag"} || confess "-protect-flag";
  }


  my $VERBOSE = $FLAGS{"aa-verbose"};
  my $found;
  my $hits;
  my $pubmed_list = $options{"-pubmed-list"};
  my $pubmed_field = $options{"-pubmed-field"};
  my $pubmed_debug = $options{"-pubmed-debug"};

  my @passed;
  my $perfect;

  if (not($found) and
      $genomic_lookup and
      $db->enable_literal_variant_lookup() and
      $hits = $db->find_literal_variant("-sj" => $row)
     ) {
    # if literal variant lookup mode is enabled, for indels the
    # reference and variant bases must match exactly.
    # Primarily used for committee matches.
    $found = 1;
    $perfect = 1;
    my $request = $medal_override ? &$medal_override($hits, $medal_perfect, $perfect) : $medal_perfect;
    $medal = request_medal($medal, $request) unless $no_medal or $silent_and_suppressed;
    push @{$reasons}, $label_perfect;
    push @passed, @{$hits};
  }

  if (not($found) and
      $aa and
      $db->enable_literal_aa_lookup() and
      $hits = $db->find_aa_literal(
				   "-gene" => $gene,
				   "-aa" => $aa,
				  )
     ) {
    #
    #  literal/trusted AA matches, e.g. SJ-to-SJ,
    #  used particuarly for committee decisions.
    #  This lookup takes priority over genomic search in case
    #  the AA is identical but the genomic lookup is slightly different.
    #
    my $handshake_ok = check_transcript_handshake(%options, "-hits" => $hits);
    if ($handshake_ok) {
      $found = 1;
      $perfect = 1;
      my $request = $medal_override ? &$medal_override($hits, $medal_perfect, $perfect) : $medal_perfect;
      $medal = request_medal($medal, $request) unless $no_medal or $silent_and_suppressed;
      push @{$reasons}, $label_perfect;
      push @passed, @{$hits};
    }
  }

  if (not($found) and $genomic_lookup) {
    #
    # genomic information available:
    #
    if ($hits = $db->find_snv(
			      "-sj" => $row
			     )) {
      # perfect lookup of SNV by genomic location
      $found = 1;
      $perfect = 1;
      my $request = $medal_override ? &$medal_override($hits, $medal_perfect, $perfect) : $medal_perfect;
      $medal = request_medal($medal, $request) unless $no_medal or $silent_and_suppressed;
      push @{$reasons}, $label_perfect;
      push @passed, @{$hits};
    } elsif ($hits = indel_search($db, $row)) {
      # perfect (within tolerances) match to indel by genomic location
      $found = 1;
      $perfect = 1 unless $db->enable_literal_variant_lookup();
      # in literal variant mode match will already have been made above.
      # for committee variants, a perfect match triggers assignment
      # of a draft committee call (intended e.g. if a different
      # variant triggers the same AA call).  However the indel
      # matching rules are inherently fuzzier, e.g. +/- 3 nt.
      # We still want to detect these kinds of matches for committee
      # variants, however they should NOT be considered perfect because
      # we don't want to speak for the committee if the indel
      # is different (i.e. even if it overlaps, unless the AA matches
      # too we shouldn't consider it an identical match).
      my $request = $medal_override ? &$medal_override($hits, $medal_perfect, $perfect) : $medal_perfect;
      $medal = request_medal($medal, $request) unless
	$no_medal or
	  ($is_silent and not($perfect)) or $silent_and_suppressed;
      # if variant is silent, only request a medal if match is perfect
      push @{$reasons}, $label_perfect;
      push @passed, @{$hits};
    }
  }

  if (not($found) and
      $aa and
      $hits = $db->find_aa_substitution(
					"-gene" => $gene,
					"-aa" => $aa,
				       )) {
    #
    # perfect AA substitution:
    #
    my $handshake_ok = check_transcript_handshake(%options, "-hits" => $hits);
    #    die "found substitution $aa handshake=$handshake_ok";
    if ($handshake_ok) {
      $found = 1;
      $perfect = 1;
      # some databases only provide AA annotations (e.g. UMD),
      # so if we don't consider this category of match perfect the
      # annotation callback below will never be invoked
      my $request = $medal_override ? &$medal_override($hits, $medal_perfect, $perfect) : $medal_perfect;
      $medal = request_medal($medal, $request) unless $no_medal or $silent_and_suppressed;
      push @{$reasons}, $label_perfect;
      push @passed, @{$hits};
    }
  }

  #
  #  begin IMPERFECT matches:
  #
  if (not($found) and
      $genomic_lookup and
      $hits = $db->find_snv_site("-sj" => $row)) {
    # see if SJ site matches database variant position alone.
    my $request = $medal_override ? &$medal_override($hits, $medal_site, $perfect) : $medal_site;
    $medal = request_medal($medal, $request) unless $no_medal or $is_silent;
    $found = 1;
    push @{$reasons}, $label_perfect . "_site_only";
    push @passed, @{$hits};
  }

  if (not($found) and
      $aa and
      $hits = $db->find_aa_specific_codon(
					  "-gene" => $gene,
					  "-aa" => $aa,
					 )) {
    #
    # match to specific codon
    #
    my $handshake_ok = check_transcript_handshake(%options, "-hits" => $hits);
    if ($handshake_ok) {
      $found = 1;
      #      dump_die($hits->[0], "codon-only hit for $aa", 1);
      my $request = $medal_override ? &$medal_override($hits, $medal_codon, $perfect) : $medal_codon;
      $medal = request_medal($medal, $request) unless $no_medal or $is_silent;
      # MW "Definitions and TP53a.pptx", 10/28/2013:
      # "Silver:
      # A variant at the same genomic location as a gold germline variant"
      # HOWEVER: not sure if this applies to all databases
      push @{$reasons}, $label_perfect . "_codon_only";
      push @passed, @{$hits};
    }
  }

  if (not($found) and
      $aa and
      $hits = $db->find_aa_codons(
				  "-gene" => $gene,
				  "-aa" => $aa,
				 )) {
    #
    #  broad codon overlap
    #
    my $handshake_ok = check_transcript_handshake(%options, "-hits" => $hits);
    if ($handshake_ok) {
      $found = 1;
      my $request = $medal_override ? &$medal_override($hits, $medal_codon, $perfect) : $medal_codon;
      $medal = request_medal($medal, $request) unless $no_medal or $is_silent;
      # FIX ME: does this deserve a medal?  separate/lesser class?
      push @{$reasons}, $label_perfect . "_broad_codon_overlap";
      push @passed, @{$hits};
    }
  }

  if (@passed and $pubmed_field) {
    die unless $pubmed_list;
    cluck sprintf "***** ERROR: missing pubmed field %s for %s!", $pubmed_field, $label_perfect unless exists $passed[0]->{$pubmed_field};

    my @pm = grep {$_} map {$_->{$pubmed_field}} @passed;
    #    printf STDERR "found PubMed for %s: %s\n", $label_perfect, join ",", @pm if @pm;
    push @{$pubmed_list}, @pm;
    if ($pubmed_debug) {
      foreach my $ref (@passed) {
	dump_die($ref, "PubMed lookup debug $pubmed_field", 1);
      }
    }
  }

  if (@passed and $hit_callback) {
    &$hit_callback(
		   "-reasons" => $reasons,
		   "-hits" => \@passed,
		   "-label" => $label_perfect
		  ) if $hit_callback_perfect ? $perfect : 1;
  }

  if ($found) {
    $$protect_flag = 1 if $protect;
  }

  if (@passed and ($debug_hits or $FLAGS{"debug-db-hit"})) {
    printf STDERR "DEBUG: %d hits to %s, perfect=%d\n", scalar(@passed), $label_perfect, $perfect || 0;
    foreach my $h (@passed) {
      dump_die($h,"dump hit", 1);
    }
  }

  if ($options{"-return-hits-and-perfect"}) {
    # ugh
    return (\@passed, $perfect);
  } else {
    return $medal;
  }
}

sub parse_damaging {
  # see http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFP2.1.readme.txt
  my ($row) = @_;
  my ($call_available, $is_damaging, $is_tolerated);

  my $is_broken;

  my @all_fields;

  if ($GERMLINE_USE_MUTATIONASSESSOR_OVER_SIFT) {
    @all_fields = (
		   "MutationAssessor_pred",
		   # "SIFT_pred",
		   "Polyphen2_HDIV_pred",
		   #		    "Polyphen2_HVAR_pred",
		  );
  } else {
    @all_fields = (
		   "SIFT_pred",
		   "Polyphen2_HDIV_pred",
		   #		    "Polyphen2_HVAR_pred",
		  );
  }
  # Bobby Huether 1/21/2014:
  # From my experience HumDiv is reported in literature.  I cannot
  # recall a time when I have seen both reported.  I suspect this
  # might be due to the fact that this is the default Polyphen run.
  # In addition, I am more comfortable with the HumDiv model as I can
  # "see" how the data that generate this model is applied to predict
  # the variant effect.  I feel HumVar is ok but b/c the data used to
  # train the model is built from variants loosely associated with
  # disease and you risk missing obvious deleterious mutants.  We
  # generate both predictions through the sfp pipeline but only report
  # the HumDiv.

  $call_available = 0;
  foreach my $f (@all_fields) {
    $call_available = undef unless exists $row->{$f};
  }

  my $platforms_with_call = 0;
  if (defined($call_available)) {
    # required columns are present
    foreach my $f (@all_fields) {
      my $v_raw = $row->{$f};
      if ($v_raw) {
	$call_available = 1;
	$platforms_with_call++;
	foreach my $v (split /;/, $v_raw) {
	  # for PolyPhen2, multiple entries separated by ";"
	  # (isoform-specific)
	  if ($v eq "D") {
	    # "D(amaging)"
	    $is_damaging = 1;
	  } elsif ($v eq "T" or $v eq "B") {
	    # Tolerated
	    # Benign
	    $is_tolerated = 1;
	  } elsif ($v eq "P") {
	    # possibly damaging
	    # takes up a large chunk of score space;
	    # since not definitively damaging, consider tolerated?
	    # if we don't set one of these flags that creates holes in logic
	    $is_tolerated = 1;
	  } elsif ($v eq "H" or $v eq "M") {
	    # MutationAssessor: predicted high or medium functional impact
	    # recommended by Bobby Huether 1/31/2014
	    $is_damaging = 1;
	  } elsif ($v eq "L" or $v eq "N") {
	    # MutationAssessor:
	    # - predicted non-functional impact: (L)ow or (N)eutral
	    $is_tolerated = 1;
	  } elsif ($v eq "ambiguous_dbnsfp_lookup") {
	    $is_broken = 1;
	  } elsif ($v eq ".") {
	    # null value, may still sometimes appear in lists,
	    # e.g. "B;."
	  } else {
	    die "unhandled tolerance value \"$v\" in $f";
	  }
	}
      }
    }
  }

#  printf STDERR "WARNING: conflicting damage prediction results!\n" if $is_damaging and $is_tolerated;
  # this can definitely happen!

  if ($is_broken) {
    printf STDERR "ERROR: ambiguous dbNSFP hits!! FIX ME\n";
    $call_available = $is_damaging = $is_tolerated = 0;
  }

  return ($call_available, $is_damaging, $is_tolerated);
}

sub parse_nhlbi_blacklist_BROKEN {
  #
  #  database of sites we DON'T want to report
  #  (sites found in > 3 entries)
  #
  my $blacklist = $FLAGS{"nhlbi-blacklist"} || die "-nhlbi-blacklist";
  die "hopelessly broken format";

  my $vm = get_new_vm();

  my $np = new NHLBIParser("-file" => $blacklist);
  my $count = 0;
  while (my $row = $np->next()) {
    print STDERR "." if ++$count % 25000 == 0;
    if (my $aal = $row->{aa_list}) {
      # AA tracking
      foreach my $gene (@{$row->{genes}}) {
	foreach my $aa (@{$aal}) {
	  $vm->add_aa(
		      "-row" => $row,
		      "-gene" => $gene,
		      "-aa" => $aa
		     );
	}
      }
    }

    # nucleotide tracking:
    my $ref_seq = $row->{reference_sequence} || die;
    my $ref_position = $row->{reference_position} || die;
    my $ref_base = $row->{reference_base} || die;
    foreach my $var_base (@{$row->{variant_bases}}) {
      $vm->add_snv(
		   "-row" => $row,
		   "-reference" => $ref_seq,
		   "-base-number" => $ref_position,
		   "-reference-base" => $ref_base,
		   "-variant-base" => $var_base
		  );
    }

  }

  return $vm;
}

sub parse_nhlbi_blacklist {
  #
  #  database of sites we DON'T want to report
  #  (sites found at freqency of > .001+)
  #  (...taking it again from the top)
  #
  print STDERR "parsing NHLBI blacklist (mkII)...";
  my $blacklist = $FLAGS{"nhlbi-blacklist-mk2"} || die "-nhlbi-blacklist-mk2";
  my $vm = get_new_vm();

  my $df = new DelimitedFile("-file" => $blacklist,
			     "-headers" => 1,
			    );
  my $count = 0;
  while (my $row = $df->get_hash()) {
    print STDERR "." if ++$count % 25000 == 0;
    if (0) {
      foreach (sort keys %{$row}) {
	printf "%s: %s\n", $_, $row->{$_};
      }
    }

    if (my $aa = $row->{aa}) {
      # AA tracking:
      my @genes = split /,/, $row->{genes};
      foreach my $gene (@genes) {
	$vm->add_aa(
		    "-row" => $row,
		    "-gene" => $gene,
		    "-aa" => $aa
		   );
      }
    }

    # nucleotide tracking:
    my $ref_seq = $row->{chrom} || die;
    my $ref_position = $row->{pos} || die;
    my $ref_base = $row->{reference_base} || die;
    my $var_base = $row->{variant_base} || die;

    $vm->add_snv(
		 "-row" => $row,
		 "-reference" => $ref_seq,
		 "-base-number" => $ref_position,
		 "-reference-base" => $ref_base,
		 "-variant-base" => $var_base
		);
  }

  return $vm;
}

sub parse_arup {
  my (%options) = @_;
  my $fn = $FLAGS{"gl-arup"} || die "-gl-arup";

  my $tables = parse_html_tables(
				 "-file" => $fn,
				 #				 "-require-header" => "Protein nomenclature",
				 #				 "-dump" => 1,
				 "-keep-html" => 1,
				 # required to extract PubMed
				);

  my $GENE = "RET";
  # MEN2 database -> NM_020630.4 -> RET

  my $vm = get_new_vm();
  my @raw_rows;
  foreach my $row (@{$tables->[0]}) {
    # MW 11/6/2013:
    # "For the MEN2 link the pathologic variants we should include as gold"

    $row->{PubMed} = parse_pmid($row->{"First Reference"} || die "no reference?");

    $row->{$FIELD_REFSEQ} = $ARUP_MEN2_RET_NM;
    # HACK

    if (1) {
      $row->{"Genotype (cDNA)"} = html_strip($row->{"Genotype (cDNA)"});
    } else {
      # nested tags in comment field, FAIL
      foreach my $k (keys %{$row}) {
	$row->{$k} = html_strip($row->{$k});
      }
    }
    $row->{Gene} = $GENE;

    push @raw_rows, $row;
    next unless lc($row->{Classification} || die) eq "pathogenic";
    my $aa = $row->{"Protein Change"} || next;
    # FIX ME: other entries without this annotation??
    $vm->add_aa(
		"-row" => $row,
		"-gene" => $GENE,
		"-aa" => $aa
	       );
  }
  if ($options{"-raw"}) {
    return \@raw_rows;
  } else {
    return $vm;
  }
}

sub parse_cosmic {
  my %options = @_;
  my $fn = $options{"-file"} || die "-file";
  my $all_variants = $options{"-all-variants"};
  my $all_genes = $options{"-all-genes"};
  my $wanted_genes;
  unless ($all_variants) {
    $wanted_genes = $options{"-genes"} || die "-genes";
  }
  my $VERBOSE = $FLAGS{"verbose-cosmic"};
  printf STDERR "loading COSMIC (%s)...\n", $fn;

  my %saw_genes;

  #  my $vm = get_new_vm(1);
  # test COSMIC parsing for correct reference sequence
  # - SNVs seem ok (+ and -)
  # - TO DO: multi-base substitutions on - strand!!!
  my $vm = get_new_vm();

  my $STATUS_NO = 0;
  my $STATUS_MAYBE = 1;
  my $STATUS_YES = 2;

  my %wanted_desc = (
		     "Substitution - Missense" => $STATUS_NO,
		     "Unknown" => $STATUS_MAYBE,
		     "Complex" => $STATUS_MAYBE,
		     # for unknown/complex, AA annotation may
		     # reveal frameshift status

		     "Substitution - Nonsense" => $STATUS_YES,
		     "Deletion - Frameshift" => $STATUS_YES,
		     "Insertion - Frameshift" => $STATUS_YES,
		     "Complex - frameshift" => $STATUS_YES,

		     "Nonstop extension" => $STATUS_YES,
		     # e.g. p.*2844E (run-on)

		     "Whole gene deletion" => $STATUS_NO,

		     "Deletion - In frame" => $STATUS_NO,
		     "Insertion - In frame" => $STATUS_NO,
		     "Complex - deletion inframe" => $STATUS_NO,
		     "Complex - insertion inframe" => $STATUS_NO,
		     # QUESTION: what about inframe indels???
		     "Complex - compound substitution" => $STATUS_NO,
		     # are these always inframe?

		     "No detectable mRNA/protein" => $STATUS_NO

		    );
  # variant categories we definitely want to include

  my $df = new DelimitedFile(
			     "-file" => $fn,
			     "-delimiter" => "\t",
			     "-headers" => 1,
			    );
  while (my $row = $df->get_hash()) {
    my $gene = $row->{"Gene name"};
    next unless $all_variants or $wanted_genes->{$gene};
    $saw_genes{$gene}++;

    my $aa = $row->{"Mutation AA"};
    my $desc = $row->{"Mutation Description"};

    my $wanted_desc = $wanted_desc{$desc};
    die "unhandled COSMIC Mutation Description $desc gene=$gene aa=$aa" unless defined $wanted_desc;

    my $usable;
    if ($wanted_desc == $STATUS_YES) {
      $usable = 1;
    } elsif ($wanted_desc == $STATUS_MAYBE) {
      if ($aa =~ /fs/) {
	# desc=Unknown aa=p.L35fs*10
	# desc=Complex aa=p.R158fs
	# desc=Complex aa=p.R158fs
	# desc=Unknown aa=p.L35fs*10
	# desc=Complex aa=p.R158fs
	# desc=Complex aa=p.R158fs
	# desc=Unknown aa=p.L35fs*10
	# desc=Complex aa=p.R26fs
	# desc=Complex aa=p.R26fs
	# desc=Complex aa=p.R65fs
	# desc=Complex aa=p.R65fs
	$usable = 1;
      }
    }

    $usable = 1 if $all_variants;

    if ($usable) {
      my $cds = $row->{"Mutation CDS"};
      my $pos = $row->{"Mutation GRCh37 genome position"} || die;
      # 19:58862934-58862934
      my $strand = $row->{"Mutation GRCh37 strand"} || die;
      printf STDERR "adding COSMIC %s AA:%s NT:%s:%s %s\n", $gene, $aa, $pos, $cds, $strand if $VERBOSE;

      $vm->add_aa(
		  "-gene" => $gene,
		  "-aa" => $aa,
		  "-row" => $row
		 ) if $aa;

      #
      # also add nucleotide-level tracking:
      #
      $pos =~ /^(\w+):(\d+)\-(\d+)$/ || die $pos;
      my ($chrom, $start, $end) = ($1, $2, $3);
      my $np = new NucleotideSubstitutionParser();
      $np->auto_strand_fix($strand);

      if ($np->parse($cds)) {
	my $category;
	if ($np->is_substitution()) {
	  # simple substitution
	  $category = "simple_substitution";
	  my $ref_base = $np->reference_sequence();
	  my $var_base = $np->variant_sequence();

	  $vm->add_snv(
		       "-reference" => $chrom,
		       "-base-number" => $start,
		       "-reference-base" => $ref_base,
		       "-variant-base" => $var_base,
		       "-row" => $row,
		      );
	} elsif ($np->is_deletion()) {
	  # simple deletion
	  $category = "simple_deletion";
	  $vm->add_deletion(
			    "-reference" => $chrom,
			    "-start" => $start,
			    "-end" => $end,
			    "-row" => $row
			   );
	} elsif ($np->is_insertion()) {
	  # simple insertion
	  $category = "simple_insertion";
	  $vm->add_insertion(
			     "-reference" => $chrom,
			     "-start" => $start,
			     "-end" => $end,
			     "-row" => $row,
			     "-count" => $np->event_length()
			    );
	} elsif ($np->is_complex_indel()) {
	  # complex indel, record as both an insertion
	  # and a deletion at the same site (hack)
	  $category = "complex_indel";
	  $vm->add_deletion(
			    "-reference" => $chrom,
			    "-start" => $start,
			    "-end" => $end,
			    "-row" => $row
			   );
	  $vm->add_insertion(
			     "-reference" => $chrom,
			     "-start" => $start,
			     "-end" => $end,
			     "-row" => $row,
			     "-count" => $np->complex_insertion_length()
			    );
	  # for insertions, not sure what position to use.
	  # if a multi-base indel and the event is on the - strand,
	  # should the "last" base be used instead?
	  # To cast the widest net, just use the entire deleted range.
	} else {
	  die "ERROR: unhandled event type for $cds!";
	}
	printf STDERR "adding COSMIC NT %s %s %s.%d-%d (%s)\n",
	  $gene, $cds, $chrom, $start, $end, $category if $category and $VERBOSE;
      } else {
	printf STDERR "WARNING: skipping unparsable COSMIC event %s, strand %s\n", $cds, $strand;
      }
    }
  }

  unless ($all_variants) {
    foreach my $gene (keys %{$wanted_genes}) {
      printf STDERR "raw COSMIC count for %s = %d\n", $gene, $saw_genes{$gene} || 0;
    }
  }

  return $vm;
}

sub parse_hgmd_dm_mp {
  #
  #  HGMD
  #
  #  NOTE: appears to contain some reference base inconsistencies
  #
  my (%options) = @_;
  my $file = $options{"-file"} || die "-file";

  my $df = new DelimitedFile(
			     "-file" => $file,
			     "-headers" => 1,
			    );

  my $vm;
  if (0) {
    print STDERR "DEBUG: sanity-checking HGMD reference\n";
    $vm = get_new_vm(1);
  } else {
    $vm = get_new_vm();
  }

  my $nsp = new NucleotideSubstitutionParser();

  #  my %new_types = map {$_, 1} qw(E G N R X);
  #  my %new_types = map {$_, 1} qw(E G N R);
  my %new_types = map {$_, 1} qw(E G N);

  my @f_trimmed = qw(
		      refcore
		      pmid
		   );

  while (my $row = $df->get_hash()) {
    my $type = $row->{base} || die;
    dump_die($row, "FIX ME: unhandled HGMD base $type") if $new_types{$type};

    my $save_row;
    if ($TRIM_HGMD_INFO) {
      my %save;
      @save{@f_trimmed} = @{$row}{@f_trimmed};
      $save_row = \%save;
    } else {
      $save_row = $row;
    }


    if ($row->{amino} and $row->{amino} ne "-") {
      # if this field is present, descr has info
      # "-": unparsable data, e.g. chr2.32379567 "618" (codon # only)???
      if ($vm->add_aa(
		      "-aa" => $row->{descr},
		      "-gene" => $row->{gene},
		      "-row" => $save_row
		     )) {
	#      printf STDERR "OK parsing %s in %s\n", @row{qw(aachange gene)};
      } else {
	#	printf STDERR "ERROR: can't parse HGMD amino desc %s in %s\n", @{$row}{qw(descr gene)};
	dump_die($row, "can't parse HGMD AA change in amino field");
      }
    }

    #
    #  nucleotide coordinates:
    #
    if (0) {
      print STDERR "DEBUG: HGMD base lookup disabled\n";
    } elsif (exists $row->{WU_HG19_Pos}) {
      # refactor: pre-parsed variants, possibly lifted for hg38+
      if (my $pos = $row->{WU_HG19_Pos}) {
	my $chrom = $row->{Chr} || die;
	my $ra = $row->{ReferenceAllele} || dump_die($row, "no ref allele");
	my $va = $row->{MutantAllele} || dump_die($row, "no variant allele");
	my $v = new Variant();
	$v->import_generic(
		     "-reference-name" => $chrom,
		     "-base-number" => $pos,
		     "-reference-allele" => $ra,
		     "-variant-allele" => $va,
		    );
	$vm->add_variant($v, "-row" => $save_row);
      }
    } else {
      my $is_indel;

      #
      # indels: handle separately as it's possible for a record
      # to contain both types (type X = complex indel)
      #

      if ($row->{insertion}) {
	$is_indel = 1;
	my $start = $row->{coordstart} || die;
	my $end = $row->{coordend} || die;
	my $count = length($row->{mutbase}) || die;
	$vm->add_insertion(
			   "-row" => $save_row,
			   "-reference" => $row->{chromosome},
			   "-start" => $start,
			   "-end" => $end,
			   "-count" => $count
			  );
      }

      if ($row->{deletion}) {
	$is_indel = 1;
	my $start = $row->{coordstart} || die;
	my $end = $row->{coordend} || die;
	my $count = $end - $start + 1;
	$vm->add_deletion(
			  "-row" => $save_row,
			  "-reference" => $row->{chromosome},
			  "-start" => $start,
			  "-end" => $end,
			  "-count" => $count
			 );

      }

      unless ($is_indel) {
	# substitution
	my ($chrom, $start, $end, $change, $strand, $base_ref, $base_var) = @{$row}{qw(chromosome coordstart coordend hgvs strand wildbase mutbase)};
	die "not SNV" unless $start == $end;

	if (1) {
	  # use pre-parsed wildbase/mutbase
	  if ($strand eq "-") {
	    foreach ($base_ref, $base_var) {
	      $_ = reverse_complement($_);
	    }
	  }
	  $vm->add_snv(
		       "-row" => $save_row,
		       "-reference" => $chrom,
		       "-base-number" => $start,
		       "-reference-base" => $base_ref,
		       "-variant-base" => $base_var
		      );
	} else {
	  # original version:
	  # parsing problems for e.g. 2055+18G>A
	  $nsp->auto_strand_fix($strand);
	  $nsp->parse($change) || die "can't parse $change";
	  if ($nsp->is_substitution()) {
	    my $ref_base = $nsp->reference_sequence();
	    my $var_base = $nsp->variant_sequence();
	    my $ok = $vm->add_snv(
				  "-row" => $save_row,
				  "-reference" => $chrom,
				  "-base-number" => $start,
				  "-reference-base" => $ref_base,
				  "-variant-base" => $var_base
				 );
	    #      printf STDERR "HGMD import: %s\n", join " ", $chrom, $start, $ref_base, $var_base, $strand, $ok;
	  } else {
	    dump_die($row, "HGMD: substitution parse failed for $change");
	  }
	}
      }
    }
  }
  return $vm;
}

sub check_transcript_handshake {
  my (%options) = @_;
  my $row = $options{"-row"} || die "-row";
  my $transcript_handshake = $options{"-transcript-handshake"};
  # specify the field name rather than the value, so if the
  # value is missing we can complain
  my $hits = $options{"-hits"};

  my $handshake_ok;
  if ($transcript_handshake) {
    # check required
    if (my $sj_mrna = $row->{$FIELD_REFSEQ}) {
      $sj_mrna =~ s/\.\d+$//;
      # strip version # (not present now, but just in case)
      my %saw;
      foreach my $hit (@{$hits}) {
	my $db_mrna = $hit->{$transcript_handshake};
	dump_die($hit, "database doesn't contain data in transcript field $transcript_handshake") unless $db_mrna;
	$db_mrna =~ s/\.\d+$//;
	# strip version number if present
	$saw{$db_mrna} = 1;
	if ($sj_mrna eq $db_mrna) {
	  $handshake_ok = 1;
	} else {
	  printf STDERR "handshake failed: %s vs %s\n", $sj_mrna, $db_mrna;
	}
      }
      unless ($handshake_ok) {
	my $msg = sprintf "transcript handshake failed: label=%s SJ=%s db=%s",
	  $options{"-label"}, $sj_mrna, join(",", sort keys %saw);
	cluck $msg;
	#	dump_die(\%options, $msg, 1);
      }
    } else {
      printf STDERR "ERROR: can't perform transcript handshake, no mRNA_acc field!\n";
    }
  } else {
    # not required for this database
    $handshake_ok = 1;
  }

  return $handshake_ok;
}

sub dump_raw_lists {
  #
  # export fairly raw dumps of source variants, with some extended
  # annotation (SIFT/PPH2 from dbNSFP.pm, NHLBI blacklist status).
  #
  my ($gene_list) = @_;

  my $ENABLE_UMD;
  my $ENABLE_ARUP_MEN2;
  my $ENABLE_ASU_TERT;
  my $ENABLE_CLINVAR;
  my $ENABLE_TP53;
  my $ENABLE_HGMD;
  my $ENABLE_RB1;
  my $ENABLE_NHGRI;
  my $ENABLE_LOVD;

  if ($FLAGS{"debug-one"}) {
    printf STDERR "DEBUG: limited databases!\n";
    $ENABLE_TP53 = 1;
  } else {
    # everything
    $ENABLE_UMD = 1;
    $ENABLE_ARUP_MEN2 = 1;
    $ENABLE_ASU_TERT = 1;
    $ENABLE_CLINVAR = 1;
    $ENABLE_TP53 = 1;
    $ENABLE_HGMD = 0;
    print STDERR "*** HGMD disabled: need parsing update!\n";
    $ENABLE_RB1 = 1;
    $ENABLE_NHGRI = 1;
    $ENABLE_LOVD = 1;
  }

  my $VERIFY_REF_SEQ_CLINVAR = 1;
  my $VERIFY_REF_SEQ_TP53 = 1;
  my $VERIFY_REF_SEQ_HGMD = 1;
  my $VERIFY_REF_SEQ_RB1 = 1;
  my $VERIFY_REF_SEQ_NHGRI = 1;

  if ($FLAGS{"debug-no-rsc"}) {
    foreach ($VERIFY_REF_SEQ_CLINVAR,
	     $VERIFY_REF_SEQ_TP53,
	     $VERIFY_REF_SEQ_HGMD,
	     $VERIFY_REF_SEQ_RB1,
	     $VERIFY_REF_SEQ_NHGRI) {
      $_ = 0;
    }
  }

  my $DEBUG_HGMD = 0;

  my $gold_genes = read_simple_file($gene_list);
  my %wanted_genes = map {$_, 1} @{$gold_genes};

  my $dbnsfp = get_nsfp();
  $dbnsfp->unparsable_aa_ok(1);
  # source AA data is not always available/parsable,
  # don't die if encountered

  my $vm_nhlbi;
  if ($FLAGS{"no-nhlbi"}) {
    printf STDERR "DEBUG: no NHLBI\n";
    $vm_nhlbi = new VariantMatcher();
  } else {
    printf STDERR "FIX ME: update NHLBI code to remove blacklist file!!!\n";
#    $FLAGS{"nhlbi-blacklist-mk2"} = "/nfs_exports/genomes/1/projects/ClinicalSeq/germline/nhlbi_blacklist_max_freq_0.001.tab";
    $FLAGS{"nhlbi-blacklist-mk2"} = basename("/nfs_exports/genomes/1/projects/ClinicalSeq/germline/nhlbi_blacklist_max_freq_0.001.tab");
    # TEMPORARY HACK:
    # this code should be removed and replaced with same lookup code
    # used elsewhere (binary file search)
    $vm_nhlbi = parse_nhlbi_blacklist();
  }

  my @labels = (
		"database",
		"Gene",
		$FIELD_REFSEQ,
		"genomic_acc",
		"site_info",

		"Class",
		# not sure about this one; distraction?

		"Chr",
		"start",
		"end",
		"AAchange",
		"AAchange_cleaned",
		"codon_start",
		"codon_end",
		"raw_ref_base",
		"raw_var_base",
		"genome_ref_base",
		"genome_var_base",
		"NHLBI_common",
		"dbNSFP_found",
		"PubMed",
		"variant_link",
		"SJ_note",
	       );
  push @labels, @NSFP_FIELDS;

  my $rpt = new Reporter(
			 "-file" => "meta_variants_list.tab",
			 "-delimiter" => "\t",
			 "-labels" => \@labels,
			 "-auto_qc" => 1,
			);

  my $nsp = new NucleotideSubstitutionParser();
  my $rsc = get_rsc();

  if ($ENABLE_LOVD) {
    #
    #  LOVD (APC/MSH2)
    #
    foreach my $ref (
		     [ "APC", $FLAGS{"gl-apc-flatfile"} || die ],
		     [ "MSH2", $FLAGS{"gl-msh2-flatfile"} || die ]
		    ) {
      my ($gene, $fn) = @{$ref};
      printf STDERR "processing LOVD %s (%s)...\n", $gene, $fn;

      my $df = new DelimitedFile(
				 "-file" => $fn,
				 "-headers" => 1,
				);
      while (my $row = $df->get_hash()) {
	my %r;
	my $gene = $row->{$FIELD_GENE} || die;
	$r{Gene} = $gene;
	$r{$FIELD_REFSEQ} = $row->{$FIELD_REFSEQ} || die;
	$r{database} = sprintf 'LOVD_%s', $gene;
	$r{genomic_acc} = $row->{"Genomic refseq ID"};
	$r{site_info} = $row->{DNA_change};

	foreach my $f (@NSFP_FIELDS,
		       qw(
			   Class
			   SJ_note
			   Chr
			   start
			   end
			   raw_ref_base
			   raw_var_base
			   genome_ref_base
			   genome_var_base
			   NHLBI_common
			   dbNSFP_found
			)) {
	  $r{$f} = "";
	}

	#
	#  import variant by AA:
	#
	$r{AAchange} = $row->{$FIELD_AACHANGE};
	add_cleaned_aa(\%r);
	$r{PubMed} = "";
	$r{variant_link} = "";
	$rpt->end_row(\%r);
      }
    }
  }

  if ($ENABLE_NHGRI) {
    #
    #  NHGRI BRCA1/2
    #
    foreach my $ref (
		     [ "BRCA1", $FLAGS{"nhgri-brca1"} || die ],
		     [ "BRCA2", $FLAGS{"nhgri-brca2"} || die ]
		    ) {
      my ($gene, $fn) = @{$ref};
      printf STDERR "processing NHGRI %s (%s)...\n", $gene, $fn;

      my $chrom;
      my $nm;
      if ($gene eq "BRCA1") {
	$chrom = 17;
	$nm = $NHGRI_BRCA1_NM;
      } elsif ($gene eq "BRCA2") {
	$chrom = 13;
	$nm = $NHGRI_BRCA2_NM;
      } else {
	die;
      }

      my $df = new DelimitedFile(
				 "-file" => $fn,
				 "-headers" => 1,
				);
      my $nsp = new NucleotideSubstitutionParser();
      while (my $row = $df->get_hash()) {
	my %r;
	$r{$FIELD_REFSEQ} = $nm;
	# record for transcript handshaking during later lookup
	$r{database} = sprintf 'NHGRI_%s', $gene;
	$r{Gene} = $gene;
	$r{genomic_acc} = "";
	$r{site_info} = "";
	$r{Class} = "";
	$r{SJ_note} = "";
	$r{Chr} = $chrom;

	my $hgvs = $row->{"HGVS Genomic (hg19)"} || die;

	#	printf STDERR "debug hgvs %s %s\n", $gene, $hgvs;
	populate_hgvs(\%r, $hgvs, $rsc, $VERIFY_REF_SEQ_NHGRI, $vm_nhlbi, $dbnsfp, 1);

	#
	#  import variant by AA:
	#
	my $aa = $row->{"HGVS Protein"} || die;
	$aa = trim_flanking_whitespace($aa);
	# ugh, e.g. "p.Gln148del "
	$aa = "" if $aa eq "-";

	$r{AAchange} = $aa;
	add_cleaned_aa(\%r);

	$r{PubMed} = "";
	$r{variant_link} = "";

	$rpt->end_row(\%r);
      }
    }
  }

  if ($ENABLE_UMD) {
    #
    #  umd.be variants for APC, BRCA1, BRCA2
    #
    my $umd_ff = $FLAGS{"gl-umd-flatfile"} || die "-gl-umd-flatfile";
    my $df = new DelimitedFile(
			       "-file" => $umd_ff,
			       "-headers" => 1,
			      );

    while (my $row = $df->get_hash()) {
      my %r;
      $r{database} = "UMD";
      $r{Gene} = $row->{$FIELD_GENE} || die;
      $r{variant_link} = "";
      $r{AAchange} = $row->{$FIELD_AACHANGE};
      add_cleaned_aa(\%r);
      $r{site_info} = $row->{"cDNA Nomenclature"};
      $r{$FIELD_REFSEQ} = $row->{$FIELD_REFSEQ} || die;

      $r{genomic_acc} = "";
      $r{Class} = "";
      $r{PubMed} = "";
      $r{SJ_note} = "";

      foreach my $f (@NSFP_FIELDS,
		     qw(
			 Chr
			 start
			 end
			 raw_ref_base
			 raw_var_base
			 genome_ref_base
			 genome_var_base
			 NHLBI_common
			 dbNSFP_found
		      )) {
	$r{$f} = "";
      }
      $rpt->end_row(\%r);
    }
  }

  if ($ENABLE_ARUP_MEN2) {
    #
    #  ARUP MEN2 database (RET)
    #
    my $CHR_ACC = "NC_000010.9";

    printf STDERR "processing ARUP MEN2...\n";
    my $rows = parse_arup("-raw" => 1);
    foreach my $row (@{$rows}) {
      #      dump_die($row, "debug", 1);

      my %r;
      $r{database} = "ARUP_MEN2";
      $r{Gene} = $row->{Gene} || dump_die($row, "no Gene");
      $r{Chr} = cook_chromosome_name($CHR_ACC) || die;
      # NOTE: "genomic position" field is hidden from display!
      # coordinates probably need some work, reference sequence verification,
      # etc.
      $r{$FIELD_REFSEQ} = $ARUP_MEN2_RET_NM;
      $r{genomic_acc} = "";
      # not sure if this is NM_ again??
      $r{site_info} = $row->{"Genotype (cDNA)"};
      $r{Class} = "";
      $r{SJ_note} = sprintf "Classification=%s", $row->{Classification};

      $r{AAchange} = $row->{"Protein Change"};
      $r{PubMed} = $row->{PubMed};
      $r{variant_link} = "";

      add_cleaned_aa(\%r);

      foreach my $f (@NSFP_FIELDS,
		     qw(
			 start
			 end
			 raw_ref_base
			 raw_var_base
			 genome_ref_base
			 genome_var_base
			 NHLBI_common
			 dbNSFP_found
		      )) {
	$r{$f} = "";
      }
      # requires genomic mapping for SNVs, not done yet

      $rpt->end_row(\%r);
    }
  }

  #
  #  import ASU TERT:
  #
  if ($ENABLE_ASU_TERT) {
    # STILL NEEDED:
    # - genomic mapping
    # - genomic accession

    die "ASU tert needs update";
    my $ff_tert = $FLAGS{"asu-tert-flatfile"} || die "-asu-tert-flatfile";
    printf STDERR "processing ASU TERT (%s)...\n", $ff_tert;

    my $df = new DelimitedFile(
			       "-file" => $ff_tert,
			       "-headers" => 1,
			      );
    my $TERT_GENE = "TERT";
    while (my $row = $df->get_hash()) {
      #      dump_die($row, "debug", 1);
      my %r;
      $r{database} = "ASU_TERT";
      $r{Gene} = $TERT_GENE;
      $r{$FIELD_REFSEQ} = $row->{$FIELD_REFSEQ} || die;
      $r{AAchange} = $row->{"AA substitution"} || die;
      add_cleaned_aa(\%r);
      $r{PubMed} = $row->{PMID};
      $r{variant_link} = 'http://telomerase.asu.edu/diseases.html#tert';

      my $mutation = $row->{Mutation} || "";
      $mutation =~ s/<[^>]+>//g;
      # HTML
      $r{site_info} = $mutation;
      # may be blank

      $r{genomic_acc} = "";
      # FIX ME: add to parser!!
      # not sure if translation was done with NM??:
      #
      # "The Genbank accession numbers used are NT_006576 for the
      # 41881 bp gene and NM_198253 for the cDNA and amino acid
      # sequences."

      $r{Class} = "";
      $r{SJ_note} = "SNVs need genomic mapping";

      foreach my $f (@NSFP_FIELDS,
		     qw(
			 Chr
			 start
			 end
			 raw_ref_base
			 raw_var_base
			 genome_ref_base
			 genome_var_base
			 NHLBI_common
			 dbNSFP_found
		      )) {
	$r{$f} = "";
      }
      # requires genomic mapping for SNVs, not done yet

      $rpt->end_row(\%r);
    }
  }

  #
  #  import ClinVar:
  #
  if ($ENABLE_CLINVAR) {
#    my $ff_cv = $FLAGS{"clinvar-gedi-flatfile"} || die "-clinvar-gedi-flatfile";
    die "clinvar flatfile obsolete";
    my $ff_cv = $FLAGS{"clinvar-flatfile"} || die "-clinvar-flatfile";
    printf STDERR "processing ClinVar (%s)...\n", $ff_cv;

    my %clnsig = (
		  "0" => "unknown",
		  "1" => "untested",
		  "2" => "non-pathogenic",
		  "3" => "probable-non-pathogenic",
		  "4" => "probable-pathogenic",
		  "5" => "pathogenic",
		  "6" => "drug-response",
		  "7" => "histocompatibility",
		  "255" => "other"
		 );

    my $df = new DelimitedFile(
			       "-file" => $ff_cv,
			       "-headers" => 1,
			      );
    while (my $row = $df->get_hash()) {
      my $gene = $row->{gene_symbol} || die;
      next unless $wanted_genes{$gene};
      #      dump_die($row, "debug", 1);

      my %r;
      $r{database} = "ClinVar";
      $r{Gene} = $gene;
      $r{$FIELD_REFSEQ} = "";
      $r{genomic_acc} = "";
      $r{AAchange} = "";
      $r{AAchange_cleaned} = "";
      add_cleaned_aa(\%r);
      $r{PubMed} = "";
      $r{variant_link} = "";
      my $cs = $row->{clinsig_raw};
      die unless defined $cs;
      $r{SJ_note} = sprintf "clinical_significance_code=%d (%s)", $cs, $clnsig{$cs};
      # source is a VCF file; not provided

      $r{Class} = "";
      # I think there's a bit somewhere for stops

      my $hgvs = $row->{hgvs_string} || die;
      $r{site_info} = $hgvs;
      populate_hgvs(\%r, $hgvs, $rsc, $VERIFY_REF_SEQ_CLINVAR, $vm_nhlbi, $dbnsfp, 1);

      $rpt->end_row(\%r);
    }
  }


  #
  #  import TP53:
  #
  if ($ENABLE_TP53) {
    my $ff_gl = $FLAGS{"iarc-tp53-germline"} || die "-iarc-tp53-germline";
    my $ff_somatic = $FLAGS{"iarc-tp53-somatic"} || die "-iarc-tp53-somatic";
    my $TP53_SYM = "TP53";

    foreach my $dr (
		    [ "IARC_TP53_germline", $ff_gl ],
		    [ "IARC_TP53_somatic", $ff_somatic ],
		   ) {
      # this will ultimately be replaced by Erin's GeDI version once done
      my ($db_label, $flatfile) = @{$dr};
      printf STDERR "processing %s (%s)...\n", $db_label, $flatfile;
      my $df = new DelimitedFile(
				 "-file" => $flatfile,
				 "-headers" => 1,
				);
      while (my $row = $df->get_hash()) {
	#	dump_die($row, "debug", 1);
	my %r;

	$r{database} = $db_label;
	$r{Gene} = $TP53_SYM;
	$r{$FIELD_REFSEQ} = $IARC_TP53_NM;
	$r{site_info} = $row->{g_description_hg19};

	$r{AAchange} = $row->{ProtDescription};
	add_cleaned_aa(\%r);
	$r{Class} = $row->{Effect};

	$r{SJ_note} = "";

	$r{genomic_acc} = "";
	$r{variant_link} = "";
	$r{PubMed} = "";
	# FIX ME if possible

	if (my $g_raw = $row->{g_description_hg19}) {
	  populate_hgvs(\%r, $g_raw, $rsc, $VERIFY_REF_SEQ_TP53, $vm_nhlbi, $dbnsfp);
	} else {
	  # hg19 coordinates not provided
	  $r{SJ_note} = "no hg19 field data";
	  $r{Chr} = "17";
	  # no field data, can't parse accession
	  foreach my $f (@NSFP_FIELDS,
			 qw(
			     start
			     end
			     raw_ref_base
			     raw_var_base
			     genome_ref_base
			     genome_var_base
			     NHLBI_common
			     dbNSFP_found
			  )) {
	    $r{$f} = "";
	  }
	}

	$rpt->end_row(\%r);
      }
    }
  }


  #
  #  import RB1:
  #
  if ($ENABLE_RB1) {
    my $rb1 = $FLAGS{"rb1-flatfile"} || die "-rb1-flatfile";
    printf STDERR "processing RB1 (%s)...\n", $rb1;
    my $RB1_SYM = "RB1";
    my $df = new DelimitedFile(
			       "-file" => $rb1,
			       "-headers" => 1,
			      );
    while (my $row = $df->get_hash()) {
      my %r;
      $r{database} = "RB1";
      $r{Gene} = $RB1_SYM;
      $r{$FIELD_REFSEQ} = $row->{transcript_id} || die;
      $r{genomic_acc} = $row->{genbank_id} || die;
      $r{site_info} = join " ", @{$row}{qw(g-position cDNA_change)};
      $r{Class} = $row->{expected_consequence_type};
      #      dump_die($row, "before", 1);
      $r{AAchange} = $row->{Protein};
      add_cleaned_aa(\%r);
      #      dump_die(\%r, "after");
      $r{Chr} = $row->{hg19_chrom} || die;
      $r{variant_link} = $row->{RB1_db_link};
      $r{PubMed} = $row->{PubMed};
      $r{SJ_note} = "";

      if ($row->{hg19_base_number}) {
	#
	# SNV that was successfully genomically mapped
	#
	$r{start} = $r{end} = $row->{hg19_base_number} || die;
	$r{raw_ref_base} = $row->{db_reference_base} || die;
	$r{raw_var_base} = $row->{db_variant_base} || die;
	$r{genome_ref_base} = $row->{hg19_reference_base} || die;
	$r{genome_var_base} = $row->{hg19_variant_base} || die;

	if (($row->{hg19_possible_base_flip} || 0) == -1) {
	  # -1: base at hg19 site doesn't match either the reference
	  # or the variant base specified by the database entry (!)
	  add_sj_note($row, "genomic site is neither reference or variant base");
	}

	#	dump_die($row, "debug", 1);

	my $rsc_ok = verify_ref_seq($VERIFY_REF_SEQ_RB1, $rsc, \%r);
	if ($rsc_ok) {
	  add_nhlbi(\%r, $vm_nhlbi);
	  add_nsfp(
		   "-row" => \%r,
		   "-dbnsfp" => $dbnsfp,
		   "-chr" => $r{Chr},
		   "-pos" => $r{start},
		   "-reference-base" => $r{genome_ref_base},
		   "-variant-base" => $r{genome_var_base},
		   "-aa" => $r{AAchange},
		   "-nm" => ($r{$FIELD_REFSEQ} || die),
		  );
	  $rpt->end_row(\%r);
	} else {
	  dump_die($row, "RB1 genome sanity check fail!", 1) unless (($row->{hg19_possible_base_flip} || "") eq "-1");
	  # this should be dealt with in mapping, however there are some
	  # oddities indicated by hg19_possible_base_flip: -1
	}
      } else {
	# - SNVs without genomic mappings
	# - indels
	# - complex/other
	foreach my $f (qw(
			   start
			   end
			   raw_ref_base
			   raw_var_base
			   genome_ref_base
			   genome_var_base
			   NHLBI_common
			   dbNSFP_found
			),
		       @NSFP_FIELDS
		      ) {
	  $r{$f} = "";
	}
	$rpt->end_row(\%r);
      }
    }
  }

  #
  #  import HGMD:
  #
  if ($ENABLE_HGMD) {
    die "FIX ME: HGMD PARSING UPDATE";
#    my $hgmd = $FLAGS{"hgmd-dm-mp-indels-and-splices"} || die "-hgmd-dm-mp-indels-and-splices";
    my $hgmd = $FLAGS{"hgmd-dm-mp-all"} || die "-hgmd-dm-mp-all";

    printf STDERR "processing HGMD (%s)...\n", $hgmd;
    my $df = new DelimitedFile(
			       "-file" => $hgmd,
			       "-headers" => 1,
			      );
    while (my $row = $df->get_hash()) {
      dump_die($row, "debug", 1) if $DEBUG_HGMD;
      my $gene = $row->{gene};
      next unless $wanted_genes{$gene};

      my %r;
      $r{database} = "HGMD";
      $r{Gene} = $gene;
      $r{Chr} = cook_chromosome_name($row->{chromosome});
      $r{start} = $row->{coordstart};
      $r{end} = $row->{coordend};

      $r{AAchange} = $row->{amino} ? $row->{descr} : "";
      add_cleaned_aa(\%r);

      if ($row->{refcore}) {
	$r{$FIELD_REFSEQ} = join ".", @{$row}{qw(refcore refver)};
      } else {
	die "no refGene annotation";
      }
      $r{PubMed} = $row->{pmid};
      $r{SJ_note} = "";
      $r{variant_link} = "";

      $r{genomic_acc} = "";
      $r{Class} = "";
      $r{site_info} = "";
      # not provided, but not needed

      my $ref_base = $r{raw_ref_base} = $row->{wildbase};
      my $var_base = $r{raw_var_base} = $row->{mutbase};
      die unless $ref_base or $var_base;
      # for indels one will be blank, but both shouldn't be
      foreach ($ref_base, $var_base) {
	$_ = "-" unless $_;
      }

      my $is_indel = ($row->{insertion} or $row->{deletion}) ? 1 : 0;

      my $strand = $row->{strand};
      if ($strand eq "+") {
	# ok
      } elsif ($strand eq "-") {
	$ref_base = reverse_complement($ref_base);
	$var_base = reverse_complement($var_base);
      } else {
	die;
      }

      $r{genome_ref_base} = $ref_base;
      $r{genome_var_base} = $var_base;

      if ($is_indel) {
	foreach my $f (qw(
			   NHLBI_common
			   dbNSFP_found
			),
		       @NSFP_FIELDS
		      ) {
	  $r{$f} = "";
	}
      } else {
	# these steps will only work for SNVs:
	add_nhlbi(\%r, $vm_nhlbi);
	add_nsfp(
		 "-row" => \%r,
		 "-dbnsfp" => $dbnsfp,
		 "-chr" => ($r{Chr} || die "Chr"),
		 "-pos" => ($r{start} || die "start"),
		 "-reference-base" => ($r{genome_ref_base} || die "grb"),
		 "-variant-base" => ($r{genome_var_base} || die "gvb"),
		 #		 "-aa" => ($r{AAchange} || die "no AA"),
		 "-aa" => $r{AAchange} || "",
		 # may not be present for splices
		 "-nm" => ($r{$FIELD_REFSEQ} || die "no NM"),
		);
	my $rsc_ok = verify_ref_seq($VERIFY_REF_SEQ_HGMD, $rsc, \%r);
	if ($rsc_ok) {
	  $rpt->end_row(\%r);
	} else {
	  die "refseq check failed!";
	}
      }

    }
  }

  $rpt->finish();
}

sub verify_ref_seq {
  my ($check, $rsc, $row) = @_;
  my $ok;
  if ($check) {
    $ok = $rsc->check(
		      "-ref-name" => ($row->{chrom} || $row->{Chr} || confess "chrom"),
		      "-base-number" => ($row->{start} || $row->{pos} || die "start"),
		      #			"-ref-base" => $ref_base)) {
		      "-ref-base" => substr(($row->{genome_ref_base} || die "grb"), 0 ,1),
		      # just check 1st base
		      # (might be di/tri etc.)
		     );
  } else {
    $ok = 1;
  }
  return $ok;
}

sub add_nhlbi {
  my ($row, $vm_nhlbi) = @_;
  my $blacklisted = $vm_nhlbi->find_snv(
					"-reference" => ($row->{chrom} || $row->{Chr} || die),
					"-base-number" => $row->{start},
					# hack
					"-reference-base" => $row->{genome_ref_base},
					"-variant-base" => $row->{genome_var_base}
				       );
  # just a simple SNV check:
  # more sophisticated AA-based lookup needed in medals version

  $row->{NHLBI_common} = $blacklisted ? 1 : 0;
  #    add_sj_note($row, "NHLBI_common");

  return $blacklisted;
}

sub add_nsfp {
  my (%options) = @_;
  my $row = $options{"-row"} || die "-row";
  my $dbnsfp = $options{"-dbnsfp"} || die "-dbnsfp";
  my $pre_hits = $options{"-hits"};

  #  my $fields = $options{"-fields"} || die;
  # passthrough
  my $fields = \@NSFP_FIELDS;

  my $aa = $options{"-aa"};
  #  printf STDERR "WARNING: no AA in NSFP lookup!\n" unless $aa;
  my $preferred = $options{"-preferred"};
  if ($preferred) {
    my $nm = $options{"-nm"};
    if (not($nm) or not($nm =~ /^NM_/)) {
      # if data not available or malformatted,
      # use SJ preferred isoform for this gene
      $nm = $preferred->get_preferred_isoform($row->{GeneName} || die "no GeneName");
      printf STDERR "WARNING: invalid mRNA_acc field, using SJ preferred %s\n", $nm if $nm;
      $options{"-nm"} = $nm if $nm;
    }
  }

  if (my $sj = $options{"-sj"}) {
    $options{"-chr"} = $row->{Chr} || die;
    $options{"-pos"} = get_sj_pos($row) || die;
    $options{"-reference-base"} = $row->{ReferenceAllele} || die;
    $options{"-variant-base"} = $row->{MutantAllele} || die;
  }

  foreach my $f (@{$fields}) {
    $row->{$f} = "";
  }

  my $hits = $dbnsfp->find(%options);

  if ($pre_hits) {
    # provided, i.e. tabix->VariantMatcher
    $row->{dbNSFP_raw_row_count} = scalar @{$pre_hits};
  } else {
    # old method: binary search
    $row->{dbNSFP_raw_row_count} = scalar @{$dbnsfp->raw_hits()};
  }

  my $result;

  if (@{$hits} == 1) {
    # unique match
    @{$row}{@{$fields}} = map {$hits->[0]->{$_}} @{$fields};
    $row->{dbNSFP_found} = 1;

    if ($options{"-add-aachange"}) {
      $row->{AAchange} = $dbnsfp->get_uniprot_aachange("-row" => $hits->[0]);
    }

  } elsif (@{$hits} == 0) {
    $row->{dbNSFP_found} = 0;
  } elsif (@{$hits} > 1) {
    # ambiguous results, for now just complain and leave blank
    # (dbNSFP.pm should already have reported debug info)
    printf STDERR "ERROR: ambiguous dbnsfp lookup! aa=%s nm=%s!\n",
      ($aa || ""),
	($options{"-nm"} || "");
    $row->{dbNSFP_found} = "ambiguous matches";
  }

  return $result;
}

sub add_sj_note {
  my ($row, $note) = @_;
  my @things;
  if ($row->{SJ_note}) {
    @things = split /,/, $row->{SJ_note};
  }
  push @things, $note;
  $row->{SJ_note} = join ",", @things;
}

sub get_rsc {
  return new ReferenceSanityCheck("-fasta_dir" =>
				  $FLAGS{"fasta-dir"} || die "-fasta-dir"
				 );
}

sub germline_gene_survey {
  # 12/18/2013:
  # Re-compile reference genomic location list for potentially valuable germline SNVs
  # scan various databases for entries for genes of interest
  my ($gene_list) = @_;
  # QUESTIONS:
  # - we are assuming that the HUGO gene symbols are standardized/comparable
  #   THEY ARE NOT! how to to deal with ambiguous synonyms, e.g. HNPCC
  #   maps to MLH1 and MSH2
  # - include dinucleotides? (TP53/ClinVar)
  # - new ClinVar version
  # - COSMIC meeds rebuild + revist
  # - RB1 db
  # - how to report multiple AA events for different isoforms? [hgmd etc.]
  # - standardization of AA annotations?
  # - TO DO: PCGP somatic?

  my $VERIFY_REF_SEQ_TP53 = 0;	  # checked
  my $VERIFY_REF_SEQ_HGMD = 0;	  # checked
  my $VERIFY_REF_SEQ_CLINVAR = 0; # checked
  my $VERIFY_REF_SEQ_COSMIC = 0;
  # still needs complete check, looks mostly OK

  my $STANDARDIZE_GENES = 0;
  # disabled in favor of Gang's cleaned list

  my $gss = new GeneSymbolStandardizer("-filename" => $FLAGS{"gene-info"} || die "-gene-info");

  #
  #  search genes setup:
  #
  my $lines = read_simple_file($gene_list);
  my $genes = [];
  my %wanted_loci;
  foreach my $line (@{$lines}) {
    my @f = split /\t/, $line;
    die unless @f == 2;
    my ($gene, $locus) = @f;
    push @{$genes}, $gene;
    my ($chrom, $range) = split /:/, $locus;
    $chrom = cook_chromosome_name($chrom) || die;
    my ($start, $end) = split /\-/, $range;
    push @{$wanted_loci{$chrom}}, {
				   gene => $gene,
				   chrom => $chrom,
				   start => $start,
				   end => $end,
				  };
  }
  my %wanted_genes;

  if ($STANDARDIZE_GENES) {
    #
    #  preprocess gene symbols (both in user list and parsed databases)
    #  to standardized versions
    #
    foreach my $g_raw (@{$genes}) {
      my $g_cleaned = $gss->find($g_raw);
      if (not($g_cleaned)) {
	printf STDERR "ERROR: can't find gene symbol for %s!\n", $g_raw;
	$g_cleaned = $g_raw;
	# FIX ME
      } elsif ($g_raw ne $g_cleaned) {
	printf STDERR "WARNING: user gene %s maps to %s\n", $g_raw, $g_cleaned;
      }
      printf STDERR "map user gene %s => %s\n", $g_raw, $g_cleaned;
      printf STDERR "WARNING: apparent synonym in source list, duplicate hit to %s\n", $g_cleaned if $wanted_genes{$g_cleaned};
      $wanted_genes{$g_cleaned} = 1;
    }
  } else {
    %wanted_genes = map {$_, 1} @{$genes};
  }

  #
  #  report init:
  #
  my @labels = qw(
		   database
		   gene
		   chrom
		   pos
		   AA
		   genome_ref_base
		   genome_var_base
		);

  my $rsc = get_rsc();

  my $outfile = "germline_variants_list.tab";
  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => \@labels
			);

  my %counts;

  #
  #  COSMIC (cleaned version):
  #
  my $cc = $FLAGS{"cosmic-cleaned"} || die "-cosmic-cleaned";
  my $df = new DelimitedFile(
			     "-file" => $cc,
			     "-delimiter" => "\t",
			     "-headers" => 1,
			    );

  my $nsp = new NucleotideSubstitutionParser();
  while (my $row = $df->get_hash()) {
    my $gene = $row->{"Gene name"} || die;
    $gene = $gss->find($gene) || $gene if $STANDARDIZE_GENES;

    my %r;
    $r{database} = "COSMIC";
    $r{gene} = $gene;
    $r{AA} = $row->{"Mutation AA"};
    my $pos_s = $row->{"Mutation GRCh37 genome position"} || die;
    # 23:153008476-153008476

    my $desc = $row->{"Mutation Description"};
    next unless $desc =~ /substitution/i;
    next if $desc =~ /complex/i;
    # only simple substitutions
    die if $desc =~ /silent/i;
    # should have already been pruned from this file

    $nsp->parse($row->{"Mutation CDS"} || die "no cds") || die sprintf "parse fail for %s, desc=%s", $row->{"Mutation CDS"}, $desc;

    my $ref_base = $nsp->reference_sequence;
    my $var_base = $nsp->variant_sequence;

    my ($c_raw, $stuff) = split /:/, $pos_s;
    my ($start, $end) = split /\-/, $stuff;

    if ($start != $end) {
      printf STDERR "ignoring multi-base substitution\n";
    } else {
      # SNV
      $r{chrom} = cook_chromosome_name($c_raw);
      $r{pos} = $start;

      my $strand = $row->{"Mutation GRCh37 strand"} || die;
      if ($strand eq "+") {
	# ok
      } elsif ($strand eq "-") {
	if (1) {
	  $ref_base = reverse_complement($ref_base);
	  $var_base = reverse_complement($var_base);
	}
      } else {
	die;
      }

      $r{genome_ref_base} = $ref_base;
      $r{genome_var_base} = $var_base;

      if (my $g = is_site_wanted(\%r, \%wanted_genes, \%wanted_loci)) {
	$r{gene} = $g;
	# returns DESIRED gene symbol, e.g.
	# COSMIC entry for AB019437.1 is actually IGH
	my $rsc_ok = verify_ref_seq($VERIFY_REF_SEQ_COSMIC, $rsc, \%r);
	if ($rsc_ok) {
	  #	    printf STDERR "import OK, strand=%s\n", $strand;
	  $counts{$gene}++;
	  $rpt->end_row(\%r);
	} else {
	  printf STDERR "ERROR: import failed, strand=%s\n", $strand;
	}
      }
    }
  }

  #
  #  ClinVar (TEMPORARY):
  #
  my $vcf = $FLAGS{"clinvar-manual"} || die "-clinvar-manual";

  $df = new DelimitedFile(
			  "-file" => $vcf,
			  "-delimiter" => "\t",
			  "-headers" => 1,
			  "-skip_until" => '#CHROM',
			 );

  my %wanted_clnsig = map {$_, 1} GERMLINE_CLINVAR_CLNSIG_TYPES;

  while (my $row = $df->get_hash()) {
    my @things = split /;/, $row->{INFO} || die "no INFO";
    my %info;
    foreach my $thing (@things) {
      my @f = split /=/, $thing;
      if (@f == 1) {
	# flag
	$info{$f[0]} = 1;
      } elsif (@f == 2) {
	$info{$f[0]} = $f[1];
      } else {
	die "WTF $thing";
      }
    }

    my $cs = $info{CLNSIG} || 0;
    $cs =~ s/\|/,/g;
    # appears to sometimes be lists(s) delimited by , and/or |
    my @cs = split /,/, $cs;
    my $cs_ok = 0;
    foreach (@cs) {
      # this is super hacky: use Erin's parsed version instead!
      $cs_ok = 1 if $wanted_clnsig{$_};
    }
    if ($cs_ok) {
      if (my $str = $info{GENEINFO}) {
	my ($gene, $gi) = split /:/, $str;
	$gene = $gss->find($gene) || $gene if $STANDARDIZE_GENES;

	my $ref_base = $row->{REF} || die;
	foreach my $var_base (split /,/, ($row->{ALT} || die)) {
	  if (length($ref_base) == length($var_base)) {
	    my %r;
	    $r{database} = "ClinVar";
	    $r{gene} = $gene;
	    $r{chrom} = cook_chromosome_name($row->{"#CHROM"} || die);
	    $r{pos} = $row->{POS} || die;
	    $r{genome_ref_base} = $ref_base;
	    $r{genome_var_base} = $var_base;
	    $r{AA} = "n/a";

	    if (my $g = is_site_wanted(\%r, \%wanted_genes, \%wanted_loci)) {
	      $r{gene} = $g;
	      my $rsc_ok = verify_ref_seq($VERIFY_REF_SEQ_CLINVAR, $rsc, \%r);
	      if ($rsc_ok) {
		$counts{$gene}++;
		$rpt->end_row(\%r);
	      } else {
		die "refseq check failed!";
	      }
	    }
	  } else {
	    printf STDERR "ignoring ClinVar indel: ref=%s var=%s\n", $ref_base, $var_base;
	  }
	}

      }
    }
  }

  #
  #  HGMD DM:
  #
  my $vm_hgmd = get_vm_hgmd();
  my $hgmd_db = $vm_hgmd->db_snv();
  foreach my $key (sort keys %{$hgmd_db}) {
    my $rows = $hgmd_db->{$key};
    printf STDERR "WARNING: %d multiple HGMD rows for %s\n", scalar(@{$rows}), $key if @{$rows} > 1;
    my $row = $rows->[0];
    my $gene = $row->{gene} || die;
    $gene = $gss->find($gene) || $gene if $STANDARDIZE_GENES;

    my %r;
    $r{gene} = $gene;
    $r{database} = "HGMD_DM";

    my ($chr, $pos, $ref_base, $var_base) = split /\./, $key;
    $r{chrom} = $chr;
    $r{pos} = $pos;
    $r{genome_ref_base} = $ref_base;
    $r{genome_var_base} = $var_base;
    # FIX ME: centralize this chunk
    #    $r{AA} = $rows->[0]->{ProtDescription} || "";
    $r{AA} = $row->{descr} || "";

    if (my $g = is_site_wanted(\%r, \%wanted_genes, \%wanted_loci)) {
      $r{gene} = $g;
      my $rsc_ok = verify_ref_seq($VERIFY_REF_SEQ_HGMD, $rsc, \%r);
      if ($rsc_ok) {
	$counts{$gene}++;
	$rpt->end_row(\%r);
      } else {
	die "refseq check failed!";
      }
    }

  }

  #
  #  IARC TP53 (germline + somatic):
  #
  my $vm_tp53_gl = get_vm_tp53_gl();
  my $vm_tp53_somatic = get_vm_tp53_somatic();

  foreach my $r (
		 [ $vm_tp53_gl, "TP53_germline" ],
		 [ $vm_tp53_somatic, "TP53_somatic" ],
		) {
    my $gene = "TP53";
    die unless $wanted_genes{$gene};
    my ($vm, $db_label) = @{$r};
    my $db = $vm->db_snv();
    foreach my $key (keys %{$db}) {
      my %r;
      my $rows = $db->{$key};

      printf STDERR "WARNING: %d entries for %s\n", scalar @{$rows}, $key unless @{$rows} == 1;
      #      dump_die($rows->[0], "debug");
      my ($chr, $pos, $ref_base, $var_base) = split /\./, $key;
      if (length($ref_base) > 1) {
	printf STDERR "ignoring IARC multi-base %s\n", $key;
	next;
      }

      $r{database} = $db_label;
      $r{gene} = $gene;
      $r{chrom} = $chr;
      $r{pos} = $pos;
      $r{genome_ref_base} = $ref_base;
      $r{genome_var_base} = $var_base;
      $r{AA} = $rows->[0]->{ProtDescription} || "";

      if (my $g = is_site_wanted(\%r, \%wanted_genes, \%wanted_loci)) {
	$r{gene} = $g;
	my $rsc_ok = verify_ref_seq($VERIFY_REF_SEQ_TP53, $rsc, \%r);
	if ($rsc_ok) {
	  $counts{$gene}++;
	  $rpt->end_row(\%r);
	} else {
	  die "refseq check failed!";
	}
      }
    }
  }

  $rpt->finish();

  printf "parsed counts:\n";
  my $missing = 0;
  foreach my $g (sort keys %wanted_genes) {
    my $count = $counts{$g} || 0;
    printf "  %s: %d\n", $g, $count;
    $missing++ if $count == 0;
  }
  printf "total_genes:%d  missing_genes: %d\n",
    scalar(keys %wanted_genes),
      $missing;

  condense_output($outfile);
}

sub get_vm_tp53_gl {
  printf STDERR "loading IARC (germline)...\n";
  return parse_iarc("-file" => $FLAGS{"iarc-tp53-germline"} || die "-iarc-tp53-germline");
}

sub get_vm_tp53_somatic {
  printf STDERR "loading IARC (somatic)...\n";
  return parse_iarc("-file" => $FLAGS{"iarc-tp53-somatic"} || die "-iarc-tp53-somatic");
}

sub get_vm_hgmd {
  # my $vm_hgmd = parse_hgmd("-file" => $FLAGS{"hgmd-dm"} || die "-hgmd-dm");
  # 1st iteration
  #  my $fn = $FLAGS{"hgmd-dm-mp"} || die "-hgmd-dm-mp";
  # 2nd iteration: SNVs only
  #  my $fn = $FLAGS{"hgmd-dm-mp-indels"} || die "-hgmd-dm-mp-indels";
  # 3rd iteration: full dump of both SNVs and indels
  #  my $fn = $FLAGS{"hgmd-dm-mp-indels-and-splices"} || die "-hgmd-dm-mp-indels-and-splices";
  # 4th iteration: SNVs, indels, and splice variants
  my $vm_hgmd;
  if ($FLAGS{"no-hgmd"}) {
    printf STDERR "HGMD: disabled\n";
    $vm_hgmd = get_new_vm();
  } else {
    my $fn = $FLAGS{"hgmd-dm-mp-all"} || die "-hgmd-dm-mp-all";
    printf STDERR "loading HGMD (%s)...\n", $fn;
    $vm_hgmd = parse_hgmd_dm_mp("-file" => $fn);
    # 2nd iteration: new view with expanded info
  }
  return $vm_hgmd;
}

sub is_site_wanted {
  my ($row, $wanted_genes, $wanted_loci) = @_;
  my $chr = cook_chromosome_name($row->{chrom});
  my $pos = $row->{pos} || die;
  my $wanted = 0;

  my $gene = $row->{gene} || die;
  if ($wanted_genes->{$gene}) {
    $wanted = $gene;
  } else {
    foreach my $wanted_range (@{$wanted_loci->{$chr}}) {
      if ($pos >= $wanted_range->{start} and
	  $pos <= $wanted_range->{end}) {
	$wanted = $wanted_range->{gene} || die;
	#	dump_die($row, "test $chr $pos $gene " . $wanted_range->{gene});
	last;
      }
    }
  }
  return $wanted;
}

sub condense_output {
  my ($infile) = @_;
  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			    );

  my %rows;
  my @headers;
  while (my $row = $df->get_hash()) {
    my $key = join "_", @{$row}{qw(chrom pos genome_ref_base genome_var_base)};
    push @{$rows{$key}}, $row;
    @headers = keys %{$row} unless @headers;
  }

  my $rpt = $df->get_reporter(
			      "-file" => $infile . ".condensed.tab"
			     );
  foreach my $key (sort keys %rows) {
    my %v;
    foreach my $row (@{$rows{$key}}) {
      foreach my $k (@headers) {
	$v{$k}{$row->{$k}} = 1;
      }
    }

    my %r;
    foreach my $k (@headers) {
      $r{$k} = join ",", sort keys %{$v{$k}};
    }
    $rpt->end_row(\%r);
  }

  $rpt->finish();
}

sub parse_rb1 {
  my ($dir) = @_;

  my @data_html;
  my $overview_html;

  my $type = $FLAGS{"rb1-type"} || die "-rb1-type [unique|all]";
  die "type must be unique or all" unless $type eq "unique" or $type eq "all";
  my $all_mode = $type eq "all";

  die "TO DO: ADD GENE NAME COLUMN!";
  die "TO DO: ENSURE REF AND VAR BASES ARE DIFFERENT";

  my @EXTRA_HEADERS = (
		       "transcript_id",
		       "genbank_id",
		       "Protein_raw",
		       "Protein_modified",
		       "RB1_db_link",
		      );
  push @EXTRA_HEADERS, "PubMed" if $all_mode;

  foreach my $f (glob($dir . "/*.htm")) {
    my $bn = basename($f);
    if ($bn eq "rb1_homepage.htm") {
      $overview_html = $f;
    } elsif ($bn =~ /rb1_(\w+)_page\d+\.htm/) {
      my $ht = $1;
      if ($ht eq "unique") {
	push @data_html, $f if $type eq "unique";
      } elsif ($ht eq "all") {
	push @data_html, $f if $type eq "all";
      } else {
	die;
      }
    }
  }
  die "can't identify homepage and data pages" unless $overview_html and @data_html;

  #
  #  parse general info from homepage:
  #
  my $tables = parse_html_tables(
				 "-file" => $overview_html,
				 "-headers" => 0,
				 #				 "-dump" => 1,
				);

  my @wanted = grep {$_ and @{$_} and $_->[0]->[0] =~ /general information/i} @{$tables};
  die unless @wanted == 1;
  # general info table

  my %info;
  foreach my $row (@{$wanted[0]}) {
    $info{$row->[0]} = $row->[1];
  }

  my $nm = $info{"Transcript refseq ID"} || die;

  #
  #  parse data rows:
  #
  my @final_rows;
  my $total_entries;
  my $final_headers;
  my $genbank_id;

  foreach my $f (@data_html) {

    open(TMP, $f) || die;
    while (<TMP>) {
      if (/Variation at DNA level according to GenBank accession number (\w+)/) {
	$genbank_id = $1;
      }
    }

    my @po = (
	      "-file" => $f,
	      "-headers" => 0,
	      # won't work here
	      #	"-require-header" => "g-position"
	     );

    my $tables = parse_html_tables(@po);
    my $tables_html = parse_html_tables(@po, "-keep-html" => 1);

    #
    #  extract headers
    #
    my @headers;
    my $enabled;
    foreach my $rows (@{$tables}) {
      next unless @{$rows};
      #      printf STDERR "new table\n";
      #      printf STDERR "  %d rows\n", scalar @{$rows};
      #      printf STDERR "  %s\n", join ",", @{$rows->[0]} if @{$rows};

      if (@{$rows} == 1) {
	my $thing = $rows->[0]->[0];
	if ($all_mode) {
	  $enabled = 1 if $thing =~ /Path\./;
	  if ($thing =~ /^(\d+) public entries/) {
	    # "all contents" view: appears first
	    $total_entries = $1;
	    last;
	  }
	} else {
	  $enabled = 1 if $thing =~ /g-position/;
	  if ($thing =~ /^(\d+) entries/) {
	    # "unique variants" view: appears at end
	    $total_entries = $1;
	    last;
	  }
	}
	$thing =~ s/[^\w\-\#\.]+/_/g;
	push @headers, $thing if $enabled;
      }
    }
    printf STDERR "headers: %s\n", join ",", @headers;
    $final_headers = [ @headers ] unless $final_headers;

    my $data_table;
    my $which = $all_mode ? $tables_html : $tables;
    foreach my $rows (@{$which}) {
      next unless @{$rows};
      if (@{$rows} > 3) {
	die "argh" if $data_table;
	$data_table = $rows;
      }
    }

    foreach my $r (@{$data_table}) {
      if (@headers == @{$r}) {
	# column counts match
	my %r;
	@r{@headers} = @{$r};

	my $usable = 1;

	printf STDERR "RAW:%s\n", join "\t", @{$r};

	if ($all_mode) {
	  $r{PubMed} = $r{Reference} =~ /nih\.gov\/pubmed\/(\d+)/ ? $1 : "";
	  my $blank = 1;
	  printf STDERR "PATH=%s\n", $r{"Path."};

	  foreach my $h (@headers) {
	    $r{$h} =~ s/<[^>]+>//g;
	    # crude HTML strip
	    $r{$h} =~ s/^\s+//;
	    $r{$h} =~ s/\s+$//;
	    $blank = 0 if $r{$h};
	  }
	  $usable = 0 if $blank;
	}


	if ($usable) {
	  $r{RB1_db_link} = 'http://rb1-lovd.d-lohmann.de/variants.php?select_db=RB1&action=search_all&search_Variant%2FDNA=' . $r{cDNA_change};
	  # generate callback link, e.g.
	  # http://rb1-lovd.d-lohmann.de/variants.php?select_db=RB1&action=search_all&search_Variant%2FDNA=c.-1865_1864insCTGATA

	  push @final_rows, \%r;
	}
      } else {
	die "header/data sync problem";
      }
    }
  }

  printf STDERR "headers: %s\n", join ",", @{$final_headers};

  die unless $genbank_id;

  die sprintf "expected %d entries, parsed %d", $total_entries, scalar @final_rows unless $total_entries == @final_rows;

  my $outfile = sprintf "export_RB1_%s.tab", $type;

  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => [
				       @{$final_headers},
				       @EXTRA_HEADERS
				      ]
			);
  foreach my $row (@final_rows) {

    if (my $filter = $FLAGS{"rb1-hack"}) {
      next unless $row->{"g-position"} eq $filter;
    }

    foreach my $h (@{$final_headers}) {
      $row->{$h} = "" unless defined $row->{$h} and $row->{$h} =~ /\w/;
      $row->{$h} =~ tr/\040-\176/ /c;
      # replace non-printable characters with spaces
    }
    $row->{transcript_id} = $nm;
    $row->{genbank_id} = $genbank_id;

    my $consequence = $row->{expected_consequence_type} || "";
    my $rna_change = $row->{RNA_change} || "";
    $row->{Protein_raw} = $row->{Protein};

    if ($consequence eq "nonsense" or
	lc($rna_change) eq "stop gained") {
      # protein change should reflect presence of stop codon.
      # for the latter, not marked nonsense (???) however:
      # NM_000321.2     L11910  g.2150G>T       c.91G>T Substitution    01_ex   Stop gained     p.Glu31X                unknown RB1_01492       somatic G       T      13       48878139        G       T       0       2150
      # - raw codon is GAG
      # - variant codon is TAG (stop)
      my $aap = new AAParser();
      my $aa = $row->{Protein};
      if (not($aa) and $rna_change =~ /^p\./) {
	printf STDERR "ERROR: Protein annotation blank, found in RNA_Change field for %s\n", $rna_change;
	$aa = $rna_change;
      }
      my $broken;
      if ($aa) {
	if ($aap->parse_substitution($aa)) {
	  my $codon_reference = $aap->codon_reference;
	  my $codon_number = $aap->codon_number;
	  my $codon_variant = $aap->codon_variant;
	  # some variants already a stop, e.g. p.Leu200*
	  # however others need reformatting, e.g. p.(Leu134*)
	  die "variant is $codon_variant in $aa" unless $codon_variant eq "X" or $codon_variant eq "*";
	  # they report as X instead of * for some reason
	  my $aa_new = join "", $codon_reference, $codon_number, "*";
	  $row->{Protein} = $aa_new;
	} else {
	  $broken = 1;
	}
      } else {
	$broken = 1;
      }

      if ($broken) {
	printf STDERR "ERROR: can't check AA annotation for nonsense, consequence=%s rna_change=%s AA=%s CDS=%s gpos=%s ",
	  $consequence, $rna_change, $aa, $row->{"cDNA_change"}, $row->{"g-position"};
	if ($aa) {
	  printf STDERR "(unparsable)\n";
	} else {
	  print STDERR "(blank)\n";
	}
      }
    }

    $row->{Protein_modified} = $row->{Protein} eq $row->{Protein_raw} ? 0 : 1;

    $rpt->end_row($row);
  }
  $rpt->finish();
}

sub map_rb1 {
  my ($infile) = @_;

  my $fa = Bio::SeqIO->new(
			   "-file" => ($FLAGS{"rb1-fasta"} || die "-rb1-fasta"),
			   "-format" => "fasta",
			  );
  my $RB1_CHROM = 13;
  my $RB1_START = 48875985;
  my $RB1_END = 49057518;
  # BLAT start and end to gate results

  my $rsc = get_rsc();

  my @sequences;
  while (1) {
    my $seq = $fa->next_seq() || last;
    push @sequences, $seq;
  }
  die unless @sequences == 1;
  my $src_sequence = $sequences[0]->seq;

  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			    );


  my @extra = qw(
		  db_reference_base
		  db_variant_base
		  hg19_chrom
		  hg19_base_number
		  hg19_reference_base
		  hg19_variant_base
		  hg19_possible_base_flip
		  genbank_base_number
	       );

  my $outfile = $infile . ".map.tab";

  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => \@extra
			     );

  my %wanted_consequence = map {$_, 1} qw(
					   missense
					   nonsense
					);
  my $nsp = new NucleotideSubstitutionParser();
  while (my $r = $df->get_hash()) {
    #    if (lc($r->{Type}) eq "substitution" and
    #	$wanted_consequence{$r->{expected_consequence_type}}
    #       ) {
    foreach my $f (@extra) {
      $r->{$f} = "";
    }
    $r->{hg19_chrom} = $RB1_CHROM;

    my $gpos = $r->{"g-position"};
    if (my $pattern = $FLAGS{"rb1-hack"}) {
      next unless $gpos eq $pattern;
    }

    if ($gpos =~ /^g\.\s+\d+/) {
      my $before = $gpos;
      printf STDERR "WARNING: patching busted gpos %s => ", $gpos;
      $gpos =~ s/\s+//g;
      printf STDERR "%s\n", $gpos;
    }


    my $usable;
    my ($seq_ref, $seq_var, $pos);
    if ($gpos =~ /^g\.(\d+)$/) {
      # genomic location has base only, doesn't contain bases
      # e.g. g.41946
      $pos = $1;
      if (my $change = $r->{cDNA_change}) {
	$change =~ s/\s+.*reported \d times.*//i;
	# e.g.:
	# c.2073G>T   (Reported 2 times)
	if ($nsp->parse($change) and $nsp->event_length() == 1) {
	  $seq_ref = uc($nsp->reference_sequence);
	  $seq_var = $nsp->variant_sequence;
	  $usable = 1;
	} else {
	  printf STDERR "ERROR: gpos=%s change=%s, can't determine bases\n", $gpos, $change;
	}
      }
    } elsif ($nsp->parse($gpos)) {
      # hopefully contains genome base number as well as ref/variant bases
      # e.g. g.2124C>A
      if ($nsp->is_substitution()) {
	die unless $nsp->event_length == 1;
	$seq_ref = uc($nsp->reference_sequence);
	$seq_var = $nsp->variant_sequence;
	$pos = $nsp->start;
	$usable = 1;
      } else {
	printf STDERR "ERROR: %s is not a substitution\n", $gpos;
      }
    } else {
      printf STDERR "ERROR: can't parse g-position $gpos\n";
    }

    if ($usable) {
      if (uc(substr($src_sequence, $pos - 1, 1)) eq $seq_ref) {
	printf STDERR "YAY: %s\n", join " ", $gpos, $pos, $seq_ref, $seq_var, substr($src_sequence, $pos - 1, 1);
	my $map = new MapSNVToGenome();
	if (my $mapped_base = $map->find(
					 "-sequence" => $src_sequence,
					 "-base-number" => $pos,
					 "-base" => $seq_ref,
					 "-chrom" => $RB1_CHROM,
					 "-min-start" => $RB1_START,
					 "-max-end" => $RB1_END,
					)
	   ) {
	  $r->{hg19_base_number} = $mapped_base;
	  $r->{genbank_base_number} = $pos;

	  $r->{db_reference_base} = $seq_ref;
	  $r->{db_variant_base} = $seq_var;
	  # raw reference and variant bases from the database entry

	  if (
	      $rsc->check("-ref-name" => $RB1_CHROM,
			  "-base-number" => $mapped_base,
			  "-ref-base" => $seq_ref
			 )) {
	    print STDERR "YAY2 $gpos $pos => $mapped_base\n";
	    $r->{hg19_possible_base_flip} = 0;
	    $r->{hg19_reference_base} = $seq_ref;
	    $r->{hg19_variant_base} = $seq_var;
	  } else {
	    # user reference base doesn't match
	    my $actual = uc($rsc->last_reference_nt);
	    if ($actual eq uc($seq_var)) {
	      # actual reference base at mapped position is
	      # user variant base rather than user reference base
	      $r->{hg19_possible_base_flip} = 1;
	      printf STDERR "UGH $gpos $pos => $mapped_base but sanity failed (variant), actual=" . $rsc->last_reference_nt . "\n";
	      $r->{hg19_reference_base} = $seq_var;
	      $r->{hg19_variant_base} = $seq_ref;
	      # swap bases for mapping
	    } else {
	      # actual reference is not user reference or variant!
	      $r->{hg19_possible_base_flip} = -1;
	      $r->{hg19_reference_base} = $seq_ref;
	      $r->{hg19_variant_base} = $seq_var;
	      printf STDERR "UGH $gpos $pos => $mapped_base but sanity failed (3rd base!), actual=" . $rsc->last_reference_nt . "\n";
	    }
	  }
	} else {
	  printf STDERR "ERROR: map failed\n";
	}
      } else {
	printf STDERR "ERROR: SANITY CHECK FAILED: %s\n", join " ", $gpos, $pos, $seq_ref, $seq_var, substr($src_sequence, $pos - 1, 1);
      }
    } else {
      printf STDERR "can't parse $gpos\n";
    }
    $rpt->end_row($r);
  }
  $rpt->finish();

  sanity_check_ordering($outfile);
}

  sub sanity_check_ordering {
    my ($fn) = @_;
    my $df = new DelimitedFile(
			       "-file" => $fn,
			       "-headers" => 1,
			      );
    my @rows;
    while (my $row = $df->get_hash()) {
      my $bn = $row->{hg19_base_number} || next;
      push @rows, $row;
    }

    my $last_bn = -1;
    foreach my $row (sort {$a->{genbank_base_number} <=> $b->{genbank_base_number}} @rows) {
      my $bn = $row->{hg19_base_number} || die;
      printf "genbank:%d genome:%d\n", $row->{genbank_base_number}, $bn;
      if ($last_bn != -1) {
	die "order fail last=$last_bn this=$bn" unless $bn >= $last_bn;
      }
      $last_bn = $bn;
    }
  }

sub log_msg {
  my ($msg) = @_;
  printf STDERR "LOG %s: %s\n", scalar(localtime()), $msg;
}

sub get_infile_and_outfile {
  my ($thing) = @_;
  my ($infile, $target) = split /\t/, $thing;
  my $outfile;
  my $suffix = ".medals.tab";

  if ($target) {
    if (-d $target) {
      # write output to this directory
      $outfile = sprintf '%s/%s%s', $target, basename($infile), $suffix;
    } else {
      # full output filename
      $outfile = $target;
    }
  } elsif ($FLAGS{"outfile-fq"}) {
    $outfile = $infile . $suffix;
  } else {
    $outfile = basename($infile) . $suffix;
  }
  die unless $outfile;
  return ($infile, $outfile);
}

sub outfile_safety_check {
  my ($list) = @_;

  my %all;
  foreach my $raw (@{$list}) {
    my ($infile, $outfile) = get_infile_and_outfile($raw);

    if ($all{$outfile}) {
      my $msg = sprintf "ERROR: output filenames not unique: %s uses the same output file as %s, so output for one file will overwrite the other.  You might need -outfile-fq or to check tabbed output spec in input list. Drive safely!",
	$infile,
	  $all{$outfile};
      if ($FLAGS{"seatbelts-are-for-sissies"}) {
	print "$msg\n";
      } else {
	die $msg;
      }
    }
    $all{$outfile} = $infile;
  }

}


sub get_nsfp {
  my @opts;
  push @opts, ("-uniprot_idmapping" => ($FLAGS{"uniprot-idmapping"} || die "-uniprot-idmapping"));
  push @opts, ("-f_pos" => $FLAGS{"tabix-dbnsfp-pos"}) if $FLAGS{"tabix-dbnsfp-pos"};

  my $dbnsfp = new dbNSFP(
			  @opts
			 );
  return $dbnsfp;
}

sub add_cleaned_aa {
  my ($row) = @_;
  confess "no AAchange field" unless exists $row->{AAchange};
  my $aa_raw = $row->{AAchange};
  my $aa_cleaned = $aa_raw;
  my $aap = new AAParser();
  my $codon_start = "";
  my $codon_end = "";
  if (my $cooked = $aap->parse_substitution($aa_raw)) {
    $aa_cleaned = sprintf 'p.%s', $cooked;
    $codon_start = $codon_end = $aap->codon_number;
  } elsif ($aap->parse($aa_raw)) {
    $codon_start = $aap->codon_start();
    $codon_end = $aap->codon_end();
  }
  $row->{codon_start} = $codon_start;
  $row->{codon_end} = $codon_end;
  $row->{AAchange_cleaned} = $aa_cleaned;
}

sub dump_clinvar {
  my $dbi = get_dbi_gedi(
			 "-type" => "research",
			 "-dev" => $FLAGS{"gedi-dev"},
			);

  #   my $query = <<EOS;
  # select distinct g.gene_symbol, a.hgvs_string, v.chromosome, v.position,
  # v.rv is_orientation_reversed, cs.clinsig_raw
  # from clinvar.cv_variant v, clinvar.cv_variant_gene vg,
  # clinvar.cv_entrez_gene g, clinvar.cv_variant_allele va,
  # clinvar.cv_allele a, clinvar.cv_clinsig cs
  # where v.cv_variant_id = vg.cv_variant_id
  # and vg.entrez_gene_id = g.entrez_gene_Id
  # and  v.cv_variant_id = va.cv_variant_id
  # and va.cv_allele_id = a.cv_allele_id
  # and cs.cv_allele_id = a.cv_allele_id
  # and VC='SNV'
  # EOS

  my $query = <<EOS;
select distinct VC, g.gene_symbol, a.hgvs_string, v.chromosome, v.position,
v.rv is_orientation_reversed, cs.clinsig_raw
from clinvar.cv_variant v, clinvar.cv_variant_gene vg,
clinvar.cv_entrez_gene g, clinvar.cv_variant_allele va,
clinvar.cv_allele a, clinvar.cv_clinsig cs
where v.cv_variant_id = vg.cv_variant_id
and vg.entrez_gene_id = g.entrez_gene_Id
and  v.cv_variant_id = va.cv_variant_id
and va.cv_allele_id = a.cv_allele_id
and cs.cv_allele_id = a.cv_allele_id
order by v.chromosome, v.position
EOS

  # order by: optimization for reference sanity check

  export_query_to_flatfile(
			   "-dbi" => $dbi,
			   "-sql" => $query,
			   "-outfile" => "export_clinvar.tab",
			  );
  $dbi->disconnect();
}

sub populate_hgvs {
  # populate row from a HGVS string
  my ($row, $g_raw, $rsc, $verify_ref, $vm_nhlbi, $dbnsfp, $no_mrna_ok) = @_;

  my @f = split /:/, $g_raw;
  my $gpos;
  if (@f == 2) {
    my $chrom_acc;
    ($chrom_acc, $gpos) = @f;
    $row->{Chr} = cook_chromosome_name($chrom_acc) || die;
  } elsif (@f == 1) {
    # NHGRI: no chrom
    die unless $row->{Chr};
    # should be already
    ($gpos) = @f;
  } else {
    die;
  }

  $gpos =~ s/\(.*\)$//;
  # strip interfering parens, e.g.
  # TP53: g.7577157_7577498del342(del intron7)

  my $nsp = new NucleotideSubstitutionParser();
  if ($nsp->parse($gpos)) {
    $row->{start} = $nsp->start() || die;
    $row->{end} = $nsp->end() || die;

    if ($nsp->is_substitution()) {
      my $elen = $nsp->event_length();
      my $seq_ref = $nsp->reference_sequence || die;
      my $seq_var = $nsp->variant_sequence || die;
      dump_die($row, sprintf("elen %d ref=%s var=%s", $elen, $nsp->reference_sequence, $nsp->variant_sequence), 1) unless $elen == 1;

      $row->{raw_ref_base} = $seq_ref;
      $row->{raw_var_base} = $seq_var;

      $row->{genome_ref_base} = $seq_ref;
      $row->{genome_var_base} = $seq_var;

      unless (verify_ref_seq($verify_ref, $rsc, $row)) {
	dump_die($row, sprintf("TP53 reference validation failed for %s.%s.%s.%s, actual=%s", @{$row}{qw(Chr start genome_ref_base genome_var_base)}, $rsc->last_reference_nt), 1);
	$row->{SJ_note} = sprintf "reference mismatch, actual=%s", $rsc->last_reference_nt;
	# 17.7579472.C.T, actual=G
	# failure to complement transcript minus-oriented C to G?
      }

      # checking only possible if bases parsable

      if ($elen == 1) {
	# these steps will only work for SNVs:
	add_nhlbi($row, $vm_nhlbi);
	my $nm = $row->{$FIELD_REFSEQ};
	die "no mRNA" if not($nm) and not($no_mrna_ok);

	add_nsfp(
		 "-row" => $row,
		 "-dbnsfp" => $dbnsfp,
		 "-chr" => $row->{Chr},
		 "-pos" => $row->{start},
		 "-reference-base" => $row->{genome_ref_base},
		 "-variant-base" => $row->{genome_var_base},
		 "-aa" => $row->{AAchange},
		 "-nm" => $nm
		);
	# FIX ME: STANDARDIZE NSFP INVOCATION
      } else {
	#
	# multi-nucleotide substitution
	#
	foreach my $f (@NSFP_FIELDS,
		       qw(
			   NHLBI_common
			   dbNSFP_found
			)) {
	  $row->{$f} = "";
	}
      }
    } elsif ($nsp->is_deletion() and $nsp->reference_sequence()) {
      # deletion
      $row->{raw_ref_base} = $row->{genome_ref_base} = $nsp->reference_sequence();
      $row->{raw_var_base} = $row->{genome_var_base} = "";

      foreach my $f (@NSFP_FIELDS,
		     "NHLBI_common",
		     "dbNSFP_found"
		    ) {
	$row->{$f} = "";
      }
      dump_die($row, "deletion reference sanity fail") unless (verify_ref_seq($verify_ref, $rsc, $row));
    } elsif ($nsp->is_insertion() and $nsp->variant_sequence()) {
      # insertion with sequence
      $row->{raw_ref_base} = $row->{genome_ref_base} = "";
      $row->{raw_var_base} = $row->{genome_var_base} = $nsp->variant_sequence();
      foreach my $f (@NSFP_FIELDS,
		     "NHLBI_common",
		     "dbNSFP_found"
		    ) {
	$row->{$f} = "";
      }
    } elsif (
	     $nsp->is_complex_indel() and
	     $nsp->reference_sequence() and
	     $nsp->variant_sequence()
	    ) {
      # complex indel with sequence specified
      printf STDERR "hey now complex variant %s\n", $gpos;
      $row->{raw_ref_base} = $row->{genome_ref_base} = $nsp->reference_sequence();
      $row->{raw_var_base} = $row->{genome_var_base} = $nsp->variant_sequence();
      dump_die($row, "deletion complex indel sanity sanity fail") unless (verify_ref_seq($verify_ref, $rsc, $row));
      foreach my $f (@NSFP_FIELDS,
		     "NHLBI_common",
		     "dbNSFP_found"
		    ) {
	$row->{$f} = "";
      }
    } else {
      #
      # not enough information to populate SNV/indel base fields
      #
      foreach my $f (@NSFP_FIELDS,
		     qw(
			 raw_ref_base
			 raw_var_base
			 genome_ref_base
			 genome_var_base
			 NHLBI_common
			 dbNSFP_found
		      )) {
	$row->{$f} = "";
      }
    }
  } else {
    printf STDERR "WARNING: can't parse gpos %s\n", $gpos;
    foreach my $f (@NSFP_FIELDS,
		   qw(
		       start
		       end
		       raw_ref_base
		       raw_var_base
		       genome_ref_base
		       genome_var_base
		       NHLBI_common
		       dbNSFP_found
		    )) {
      $row->{$f} = "";
    }
  }
}

sub parse_pmid {
  my ($string) = @_;
  return $string =~ /nih\.gov\/pubmed\/(\d+)/ ? $1 : "";
}

sub html_strip {
  # crude HTML stripped: fails miserably if nested!
  my ($string) = @_;
  $string =~ s/<[^>]+>//g;
  return $string;
}

sub collapse_meta_list {
  #
  #  uniquify entries for variants with complete genomic annotations.
  #  all other variants go into "other" list.
  #  - 2 flavors?:
  #    - unique "cleaned" AA annotation
  #    - unique complete genomic annotation
  #
  my ($infile) = @_;
  my $TYPE_AA_CLEAN = 1;
  my $TYPE_AA_DIRTY = 2;
  my $TYPE_POS = 3;

  foreach my $bucket_type ($TYPE_AA_CLEAN, $TYPE_AA_DIRTY, $TYPE_POS) {
    my $tag;
    if ($bucket_type == $TYPE_AA_CLEAN) {
      $tag = "aa_clean";
    } elsif ($bucket_type == $TYPE_AA_DIRTY) {
      $tag = "aa_all";
    } elsif ($bucket_type == $TYPE_POS) {
      $tag = "position";
    } else {
      die;
    }
    printf STDERR "processing type %s...\n", $tag;
    my $outfile_main = sprintf "%s.%s.unique.tab", $infile, $tag;
    my $outfile_other = sprintf "%s.%s.other.tab", $infile, $tag;

    my $df = new DelimitedFile(
			       "-file" => $infile,
			       "-headers" => 1,
			      );

    my @extra = qw(
		    AA_parsable
		    database_count
		    mRNA_acc_count
		 );

    my $rpt_main = $df->get_reporter("-file" => $outfile_main,
				     "-extra" => \@extra
				    );
    my $rpt_other = $df->get_reporter("-file" => $outfile_other);

    my %bucket;

    my %aa_junk = map {$_, 1} (
			       "p.?",
			       "variant",
			      );


    while (my $row = $df->get_hash()) {
      my $aa_parsable = 0;

      if (my $aa = $row->{AAchange_cleaned}) {
	my $aap = new AAParser();
	if ($aap->parse_substitution($aa)) {
	  $aa_parsable = 1;
	} elsif ($aap->parse($aa)) {
	  $aa_parsable = 1;
	} else {
	  printf STDERR "unbucketable AA %s\n", $aa unless $aa_junk{$aa};
	}
      }
      $row->{AA_parsable} = $aa_parsable;

      my $key;
      if ($bucket_type == $TYPE_AA_CLEAN or
	  $bucket_type == $TYPE_AA_DIRTY) {
	if (my $aa = $row->{AAchange_cleaned}) {
	  my $usable = 0;
	  if ($bucket_type == $TYPE_AA_DIRTY) {
	    $usable = $aa_junk{$aa} ? 0 : 1;
	    # bucket by just about any AA value.  Required because
	    # there are still some weird cases we want
	    # to keep that are not dealt with by current parsers,
	    # e.g. lists such as [p.H168L/M169I]
	    # - upside: catches unhandled cases
	    # - downside: buckets junk annotations
	  } elsif ($bucket_type == $TYPE_AA_CLEAN) {
	    #
	    # require some parsability of AA to include
	    #
	    # - upside: clean data
	    # - downside: misses some odd formatting cases (lists, etc.)
	    #
	    $usable = $row->{AA_parsable};
	  } else {
	    die;
	  }
	  $key = join ".", $row->{Gene}, $aa if $usable;
	  # only bucket if the AA code has some parsability;
	  # some sources put junk values here, e.g. RB1 "variant"
	}
      } elsif ($bucket_type == $TYPE_POS) {
	if ($row->{Chr} and
	    $row->{start} and
	    $row->{end} and
	    (
	     $row->{genome_ref_base} or $row->{genome_var_base}
	     # indels may have only one or the other
	    )) {
	  $key = join ".", map {$_ || "-"} @{$row}{
	    qw(
		Chr
		start
		end
		genome_ref_base
		genome_var_base
	     )};
	  printf STDERR "key=%s\n", $key;
	}
      } else {
	die;
      }

      if ($key) {
	push @{$bucket{$key}}, $row;
      } else {
	$rpt_other->end_row($row);
      }
    }
    $rpt_other->finish();

    foreach my $key (sort keys %bucket) {
      my $rows = $bucket{$key};
      my @columns = keys %{$rows->[0]};
      push @columns, "AA_parsable";
      my %r;
      foreach my $col (@columns) {
	my $unique = uniquify($rows, $col);
	$r{$col} = join ", ", @{$unique};
      }

      foreach my $field (qw(database mRNA_acc)) {
	my $list = uniquify($rows, $field);
	my $column = $field . "_count";
	$r{$column} = scalar @{$list};
      }

      #      dump_die(\%r, "debug", 1);

      $rpt_main->end_row(\%r);
    }
    $rpt_main->finish();

  }
}

sub uniquify {
  my ($rows, $column) = @_;
  # report values only once across a set of rows,
  # skipping empty values and preserving order.
  # may need work

  my %saw;
  my @values;
  foreach my $row (@{$rows}) {
    my $value = $row->{$column};
    if (defined($value) and $value =~ /\w/ and not($saw{$value})) {
      push @values, $value;
      $saw{$value } = 1;
    }
  }
  return \@values;
}

sub meta_reannotate {
  #
  # - reannotate unique variants list:
  #   - use SJ preferred mRNA accession
  #   - add predicted AA change if possible
  #
  my $infile = $FLAGS{"meta-reannotate"} || die;

  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			    );

  my $outfile = basename($infile) . ".reannotated.tab";

  my $rpt = new Reporter(
			 "-delimiter" => "\t",
			 "-file" => $outfile,
			 "-labels" => [
				       qw(
					   database
					   Gene
					   mRNA_acc
					   Chr
					   start
					   end
					   AAchange_orig
					   AAchange_orig_cleaned
					   AAchange
					   AAchange_cleaned
					   AAchange_discrepancy
					   genome_ref_base
					   genome_var_base
					   NHLBI_common
					   dbNSFP_found
					   PubMed
					   SJ_note
					),
				       @NSFP_FIELDS
				      ]
			);

  my $sjpi = get_sjpi();
  my $nsfp = get_nsfp();

  while (my $row = $df->get_hash()) {
    #    dump_die($row, "debug", 1);
    my $chr = $row->{Chr} || die;
    my $pos = $row->{start} || die;
    my $base_ref = $row->{genome_ref_base};
    my $base_var = $row->{genome_var_base};

    my %r = %{$row};
    delete @r{@NSFP_FIELDS};
    delete $r{dbNSFP_found};

    add_cleaned_aa(\%r);
    $r{AAchange_orig} = $r{AAchange} || "";
    $r{AAchange_orig_cleaned} = $r{AAchange_cleaned} || "";

    delete $r{AAchange};
    delete $r{AAchange_cleaned};
    # reset as we are re-annotating from dbNSFP

    my $sj_rna = $sjpi->get_preferred_isoform($row->{Gene} || die) || die;
    $r{$FIELD_REFSEQ} = $sj_rna;

    my $medal = "";

    if ($base_ref and
	$base_var and
	length($base_ref) == 1 and
	length($base_var) == 1) {
      # will work for SNVs only
      add_nsfp(
	       "-row" => \%r,
	       "-dbnsfp" => $nsfp,
	       "-chr" => ($row->{Chr} || die),
	       "-pos" => ($row->{start} || die),
	       "-reference-base" => $base_ref,
	       "-variant-base" => $base_var,
	       "-nm" => $sj_rna,
	       "-add-aachange" => 1,
	      );
    } else {
      # indel
      foreach my $f (@NSFP_FIELDS, "dbNSFP_found") {
	$r{$f} = "";
      }
    }
    $r{AAchange} = "" unless $r{AAchange};
    $r{AAchange} = $r{AAchange_orig} unless $r{AAchange};
    add_cleaned_aa(\%r);

    my $discrepancy;
    if ($r{AAchange_cleaned} and $r{AAchange_orig_cleaned}) {
      $discrepancy = $r{AAchange_cleaned} eq $r{AAchange_orig_cleaned} ? 0 : 1;
    } else {
      $discrepancy = 0;
    }
    $r{AAchange_discrepancy} = $discrepancy;

    $rpt->end_row(\%r);
  }
  $rpt->finish();

}

sub run_mini_germline {
  my $variants_db = $FLAGS{"mini-gl-run"};
  # variants to match against
  # matches to SNVs only for now!
  die unless @INPUT_GL_FILES;

  my @copy_fields = qw(
			database
			PubMed
			SJ_note
		     );

  my $dbnsfp = get_nsfp();

  my $df = new DelimitedFile(
			     "-file" => $variants_db,
			     "-headers" => 1,
			    );
  #
  #  prep SNV database:
  #
  my $vm = new VariantMatcher();
  while (my $row = $df->get_hash()) {
    my $chrom = $row->{Chr};
    my $pos = $row->{start};
    my $base_ref = $row->{genome_ref_base};
    my $base_var = $row->{genome_var_base};
    die unless $base_ref or $base_var;
    if (length($base_ref) == 1 and
	length($base_var) == 1) {
      foreach ($base_ref, $base_var) {
	die unless /\w/;
      }
      $vm->add_snv(
		   "-reference" => $chrom,
		   "-base-number" => $pos,
		   "-reference-base" => $base_ref,
		   "-variant-base" => $base_var,
		   "-row" => $row
		  );
    }
  }

  my $vm_nhlbi = parse_nhlbi_blacklist();

  foreach my $gl (@INPUT_GL_FILES) {
    printf STDERR "processing %s...\n", $gl;
    my $df = new DelimitedFile(
			       "-file" => $gl,
			       "-headers" => 1,
			      );

    my $outfile = sprintf '%s.medals.tab', basename($gl);
    my $rpt = $df->get_reporter(
				"-extra" =>
				[
				 @NSFP_FIELDS,
				 @copy_fields,
				 "dbNSFP_found",
				 #				 "dbNSFP_fuzzy",
				 "dbNSFP_raw_row_count",
				 "NHLBI_common",
				 "Medal",
				 "Reason",
				],
				"-file" => $outfile
			       );
    while (my $row = $df->get_hash()) {
      #
      #  find matches to curated variant list, report if so:
      #
      #      dump_die($row, "debug_raw", 1);
      if (my $hits = $vm->find_snv("-sj" => $row)) {
	die unless @{$hits} == 1;
	my $hit = $hits->[0];
	foreach my $f (@copy_fields) {
	  $row->{$f} = $hit->{$f};
	}
      } else {
	foreach my $f (@copy_fields) {
	  $row->{$f} = "";
	}
      }

      my $chrom = $row->{Chr} || die "Chr";
      my $pos = get_sj_pos($row);
      my $base_ref = $row->{ReferenceAllele} || die "ra";
      my $base_var = $row->{MutantAllele} || die "ma";

      #
      #  add dbNSFP annotations:
      #
      my $nm_cleaned = $row->{$FIELD_REFSEQ};
      $nm_cleaned = "" unless $nm_cleaned =~ /^NM_/;
      # will break NSFP unless valid
      # some examples: "Unknown", "MUTYH"

      #      dump_die($row, "debug", 1);

      add_nsfp(
	       "-row" => $row,
	       "-dbnsfp" => $dbnsfp,
	       "-chr" => $chrom,
	       "-pos" => $pos,
	       "-reference-base" => $base_ref,
	       "-variant-base" => $base_var,
	       "-aa" => ($row->{AAChange} || die "aac"),
	       "-nm" => $nm_cleaned,
	      );

      #
      # add NHLBI:
      #
      my $blacklisted = $vm_nhlbi->find_snv(
					    "-reference" => $chrom,
					    "-base-number" => $pos,
					    "-reference-base" => $base_ref,
					    "-variant-base" => $base_var
					   );
      $row->{NHLBI_common} = $blacklisted ? 1 : 0;

      #
      #  "medal" assignment
      #
      my $medal;
      my @reasons;
      my $nsfp_status = $row->{dbNSFP_found};
      if ($nsfp_status eq "") {
	push @reasons, "dbNSFP_no_lookup";
      } elsif ($nsfp_status eq "0") {
	my $reason = "dbNSFP_not_found";
	my $aap = new AAParser();
	if ($aap->parse_substitution($row->{AAChange})) {
	  $reason .= "_SILENT" if $aap->is_silent;
	}

	push @reasons, $reason;
      } elsif ($nsfp_status eq "1") {
	# single record found
	my $ma_damaging = ($row->{MutationAssessor_pred} || "") eq "H";
	# dbNSFP:
	# MutationAssessor's functional impact of a variant : predicted
	# functional, i.e. high ("H") or medium ("M"), or predicted
	# non-functional, i.e. low ("L") or neutral ("N").

	my @pp = grep {$_} map {$row->{$_}} qw(Polyphen2_HDIV_pred Polyphen2_HVAR_pred);
	my $pp_damaging;
	# dbNSFP:
	# Polyphen2_HDIV_pred: Polyphen2 prediction based on HumDiv, "D"
	# ("porobably damaging"), "P" ("possibly damaging") and "B"
	# ("benign"). Multiple entries separated by ";".
	# Polyphen2_HVAR_pred: Polyphen2 prediction based on HumVar, "D"
	# ("porobably damaging"), "P" ("possibly damaging") and "B"
	# ("benign"). Multiple entries separated by ";".
	foreach my $calls (@pp) {
	  foreach my $call (split /;/, $calls) {
	    $pp_damaging = 1 if $call eq "D";
	  }
	}

	if ($ma_damaging and $pp_damaging) {
	  $medal = "Gold";
	  push @reasons, "2_damaging_platforms";
	}
      } else {
	die $nsfp_status;
      }

      $medal = "unknown" unless $medal;

      $row->{Medal} = $medal;
      $row->{Reason} = join ",", @reasons;

      $rpt->end_row($row);
    }
    $rpt->finish();
  }
}

sub get_vm_bad {
  my $bad_snvs = $FLAGS{"gl-bad-snv"} || die "-gl-bad-snv";
  my $df = new DelimitedFile(
			     "-file" => $bad_snvs,
			     "-headers" => 1,
			    );
  my $vm = new VariantMatcher();
  while (my $row = $df->get_hash()) {
    $vm->add_snv(
		 "-row" => $row,
		 "-sj" => 1,
		);

    $vm->add_aa_substitution(
			     "-gene" => ($row->{GeneName} || die),
			     "-aa" => ($row->{AAChange} || die),
			     "-row" => $row,
			    );
  }
  return $vm;
}

sub get_sjpi {
#  my $sjpi = new SJPreferredIsoform("-file" => $FLAGS{"preferred-isoforms"} || die "-preferred-isoforms");
  my @opts = ("-file" => $FLAGS{"preferred-isoforms"} || die "-preferred-isoforms");
#  my $sjpi = new SJPreferredIsoform("-file" => $FLAGS{"preferred-isoforms"} || die "-preferred-isoforms");
  push @opts, "-gsm" => new_gsm_lite() if $ENABLE_SJPI_GSM;
  my $sjpi = new SJPreferredIsoform(@opts);
  return $sjpi;
}

sub get_vm_rb1 {
  my $infile = $FLAGS{"rb1-flatfile"} || die "-rb1-flatfile";
  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			    );
  # die $infile;

  my $vm = new VariantMatcher();
  my $GENE = "RB1";
  # FIX ME: this should be added to export report

  while (my $row = $df->get_hash()) {
    #    dump_die($row, "debug", 1);

    my $aa = $row->{Protein};
    if ($aa) {
      $vm->add_aa(
		  "-gene" => $GENE,
		  "-aa" => $aa,
		  "-row" => $row
		 ) unless $aa eq "variant";
    }

    if ($row->{hg19_base_number}) {
      #
      # SNV that was successfully genomically mapped
      #
      my $flip = $row->{hg19_possible_base_flip};
      die unless defined $flip;

      if ($flip == 0 or $flip == 1) {
	# 0 = genome bases are consistent with bases in RB1 genomic accession
	# 1 = genome reference/variant bases swapped in genome
	#     (probably innocuous, e.g. hg18/hg19 reference base change,
	#     or similar swap in RB1 genomic accession).
	#     Bases have been swapped to genomic in this report.
	# -1 = neither base matches, don't use for now
	my $ref_base = $row->{hg19_reference_base} || die;
	my $var_base = $row->{hg19_variant_base} || die;
	my $chrom = $row->{hg19_chrom} || die;
	my $pos = $row->{hg19_base_number} || die;

	if ($ref_base eq $var_base) {
	  printf STDERR "WARNING: RB1 ref/var bases identical!\n";
	} elsif ($var_base eq "-") {
	  # deletion
	  $vm->add_deletion(
			    "-row" => $row,
			    "-reference" => $chrom,
			    "-start" => $pos,
			    "-end" => $pos,
			   );
	} else {
	  # SNV
	  $vm->add_snv(
		       "-row" => $row,
		       "-reference" => $chrom,
		       "-base-number" => $pos,
		       "-reference-base" => $ref_base,
		       "-variant-base" => $var_base
		      );
	}
      }
    }
  }

  return $vm;
}

sub generate_cosmic_hotspot_list {
  die "needs update for TSOncoDB";
  my $gene2class = load_gene_class_files();

  my $df = new DelimitedFile(
			     "-file" => ($FLAGS{"cosmic-cleaned"} || die),
			     "-headers" => 1,
			    );

  # questions:
  # 1. use "confirmed somatic variant" only?

  my %wanted_status = (
		       "Reported in another cancer sample as somatic" => 0,
		       "Variant of unknown origin" => 0,
		       "Confirmed somatic variant" => 1
		      );

  my %bucket_aa;
  my %bucket_gpos_cds;
  my %bucket_gpos_bare;

  my $headers = $df->headers_raw();

  my @raw_rows;

  my %warned;

  while (my $row = $df->get_hash()) {
    my $status = $row->{"Mutation somatic status"};
    my $wanted = $wanted_status{$status};
    die "no preference info for $status" unless defined $wanted;

    my $gene = $row->{"Gene name"} || die;
    $gene = clean_sharp_gene_symbol($gene);
    my $gene_class = $gene2class->{$gene};
    if ($gene_class) {
      # TO DO: need formal list of oncogenes
      my ($is_ts, $is_recurrent) = parse_gene_class($gene_class);
      $wanted = 0 unless $is_recurrent;
    } else {
      unless ($warned{$gene}) {
	printf STDERR "WARNING: no TS/oncogene info for %s!\n", $gene;
	$warned{$gene} = 1;
      }
      next;
    }

    if ($wanted) {
      push @raw_rows, $row;
      my $gpos = $row->{"Mutation GRCh37 genome position"} || $row->{"Mutation genome position"} || dump_die($row, "no genome position");
      my @f = split /:/, $gpos;
      die unless @f == 2;
      my $chrom = $f[0];
      $row->{hg19_chrom} = $chrom;

      my $cds = $row->{"Mutation CDS"} || dump_die($row, "no CDS");
      my $transcript = $row->{"Accession Number"} || die;
      my $aa = $row->{"Mutation AA"};
      unless (is_valid_cosmic_aa($aa)) {
	printf STDERR "rejecting AA code \"%s\"\n", $aa;
	$aa = "";
      }

      my $key_gpos_cds = join ".", $gpos, $cds;
      my $key_gpos_bare = $gpos;
      push @{$bucket_gpos_cds{$key_gpos_cds}}, $row;
      push @{$bucket_gpos_bare{$key_gpos_bare}}, $row;

      if ($aa) {
	die $aa unless $aa =~ /\w/;
	my $key_aa = join ".", $transcript, $aa;
	push @{$bucket_aa{$key_aa}}, $row;
	# FIX ME: do this by CODON NUMBER, see e.g. SNRPD3 G96
	# which has several variant amino acids
      }
    }
  }

  #
  #  find hotspots as defined by different categories and tag:
  #
  foreach my $r (
		 [ \%bucket_aa, "aa" ],
		 [ \%bucket_gpos_cds, "gpos_cds" ],
		 [ \%bucket_gpos_bare, "gpos_base" ]
		) {
    my ($bucket, $type) = @{$r};
    foreach my $key (keys %{$bucket}) {
      my $set = $bucket->{$key};
      if (@{$set} >= $COSMIC_HOTSPOT_MIN_CONFIRMED_PATIENTS) {
	foreach my $row (@{$set}) {
	  $row->{hotspot_type}{$type} = 1;
	}
      }
    }
  }

  #
  #  write hotspot rows:
  #
  my $rpt = new Reporter(
			 "-file" => sprintf("cosmic_hotspots_%d_patients_sorted.tab", $COSMIC_HOTSPOT_MIN_CONFIRMED_PATIENTS),
			 "-delimiter" => "\t",
			 "-labels" => [
				       @{$headers},
				       "hotspot_type"
				      ],
			);

  #  foreach my $row (@raw_rows) {
  foreach my $row (sort {$a->{hg19_chrom} cmp $b->{hg19_chrom}} @raw_rows) {
    # write sorted by chrom to aid later reference base sanity checking
    if ($row->{hotspot_type}) {
      my %r = %{$row};
      $r{hotspot_type} = join ",", sort keys %{$row->{hotspot_type}};
      $rpt->end_row(\%r);
    }
  }

  $rpt->finish();
}

sub get_vm_cosmic_hotspots {
  my $infile = $FLAGS{"cosmic-hotspots"} || die "-cosmic-hotspots";
  #  my $VERIFY_REF_SEQ = 1;
  # tested 1/27/2014, OK
  my $VERIFY_REF_SEQ = 0;

  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			    );

  my $vm = new VariantMatcher();
  if ($VERIFY_REF_SEQ) {
    $vm->fasta_dir($FLAGS{"fasta-dir"} || die);
  }

  printf STDERR "loading COSMIC hotspots (%s)...", $infile;

  my $count = 0;
  my $PING = 1000;
  my %import;

  while (my $row = $df->get_hash()) {
    #    printf STDERR "." if ++$count % $PING == 0;
    my $gene = $row->{"Gene name"} || die;
    my $gpos = $row->{"Mutation GRCh37 genome position"} || $row->{"Mutation genome position"} || die;
    my $strand = $row->{"Mutation GRCh37 strand"} || $row->{"Mutation strand"} || dump_die($row, "no strand");
    my $cds = $row->{"Mutation CDS"} || dump_die($row, "no CDS");
    my $transcript = $row->{"Accession Number"} || die;
    my $aa = $row->{"Mutation AA"};
    $aa = "" unless is_valid_cosmic_aa($aa);

    if ($aa) {
      $vm->add_aa(
		  "-gene" => $gene,
		  "-aa" => $aa,
		  "-row" => $row
		 );
      $import{aa}++;
    }

    my ($chrom, $range) = split /:/, $gpos;
    my @f = split /\-/, $range;
    die unless @f == 2;
    my ($start, $end) = @f;

    my $nsp = new NucleotideSubstitutionParser();

    if ($nsp->parse($cds)) {
      if ($nsp->is_substitution()) {
	my $reference_sequence = $nsp->reference_sequence();
	my $variant_sequence = $nsp->variant_sequence();
	if ($strand eq "+") {
	  # OK
	} elsif ($strand eq "-") {
	  foreach ($reference_sequence, $variant_sequence) {
	    $_ = reverse_complement($_);
	  }
	} else {
	  die;
	}

	$vm->add_snv(
		     "-row" => $row,
		     "-reference" => $chrom,
		     "-base-number" => $start,
		     "-reference-base" => $reference_sequence,
		     "-variant-base" => $variant_sequence,
		    );
	$import{substitition}++;
      } elsif ($nsp->is_deletion()) {
	$vm->add_deletion(
			  "-row" => $row,
			  "-reference" => $chrom,
			  "-start" => $start,
			  "-end" => $end
			 );
	$import{deletion}++;
      } elsif ($nsp->is_insertion()) {
	$vm->add_deletion(
			  "-row" => $row,
			  "-reference" => $chrom,
			  "-start" => $start,
			  "-end" => $end,
			  "-count" => ($nsp->event_length || die),
			 );
	$import{insertion}++;
      } elsif ($nsp->is_complex_indel) {
	printf STDERR "WARNING: not importing complex indel %s\n", join ",", $gpos, $cds, $start, $end, $nsp->complex_insertion_length(), $nsp->complex_deletion_length();
	$import{skipped_complex}++;
      } else {
	printf STDERR "WARNING: not importing event %s\n", $cds;
	$import{skipped_misc}++;
      }
    } else {
#      printf STDERR "ERROR: can't parse $cds\n";
    }
  }

  printf STDERR "COSMIC hotspot summary:\n";
  foreach (sort keys %import) {
    printf STDERR "   %s: %d\n", $_, $import{$_};
  }

  return $vm;
}

sub is_valid_cosmic_aa {
  my ($aa) = @_;
  my $valid = 1;
  $aa =~ s/^p\.//;
  $valid = 0 unless $aa =~ /\w/;
  return $valid;
}

sub base_match_test {
  my $vm = get_vm_cosmic_hotspots();

  my $hits = $vm->find_snv_site(
				"-reference" => 1,
				"-base-number" => 226252135,

				#			   "-reference-base" => "A",
				#			   "-variant-base" => "T",
			       );

  die scalar @{$hits};

}

sub tp53_hack {
  my $vm_tp53_gl = get_vm_tp53_gl();
  my $vm_tp53_somatic;
  #  my $vm_tp53_somatic = get_vm_tp53_somatic();
  die unless @INPUT_GL_FILES;

  foreach my $glf (@INPUT_GL_FILES) {
    printf STDERR "processing %s...\n", $glf;
    my $df = new DelimitedFile(
			       "-file" => $glf,
			       "-headers" => 1
			      );

    while (my $row = $df->get_hash()) {
      #      dump_die($row, "start variant lookup", 1);

      my @reasons;
      my $medal;

      my $is_novel_these_db;
      my $medal_call_protected_from_nhlbi;

      foreach my $ref (
		       [ $vm_tp53_gl, "IARC_TP53_germline" ],
		       #		       [ $vm_tp53_somatic, "IARC_TP53_somatic" ],
		      ) {
	my ($db, $reason) = @{$ref};
	#
	#  French TP53 databases (germline and somatic):
	#  - sometimes we have the AA annotation (preferred)
	#  - use genomic coordinates if not available
	#
	$medal = assign_aa_match(
				 "-vm" => $db,
				 "-row" => $row,
				 "-transcript-handshake" => $FIELD_REFSEQ,
				 "-current-medal" => $medal,
				 "-reasons" => \@reasons,
				 "-label" => $reason,
				 "-novel-flag" => \$is_novel_these_db,
				 "-medal-perfect" => CLASS_GOLD,
				 "-medal-codon" => CLASS_SILVER,
				 "-medal-site" => CLASS_SILVER,
				 "-protect" => \$medal_call_protected_from_nhlbi,
				 "-genomic-lookup" => 1,
				 #				 "-hit-callback" => \&tp53_germline_population_callback,
				 #				 "-hit-callback-perfect" => 1,
				);

      }

      next unless $medal;
      @reasons = "" unless @reasons;

      printf "%s\n",
	join ",", $medal, @reasons;

    }
  }
}

sub hack_cosmic_orig {
  die unless @INPUT_SNV_INDEL_FILES;
  my $vm_cosmic_somatic = parse_cosmic(
				     "-file" => ($FLAGS{"cosmic-cleaned"} || die "-cosmic-cleaned"),
				     "-all-genes" => 1,
				     "-all-variants" => 1
				    );

  foreach my $infile (@INPUT_SNV_INDEL_FILES) {
    printf STDERR "processing %s...\n", $infile;
    my $df = new DelimitedFile(
			       "-file" => $infile,
			       "-headers" => 1
			      );
    while (my $row = $df->get_hash()) {
      add_cosmic($row, $vm_cosmic_somatic);
      printf "%s\n", join " ", @{$row}{$FIELD_VPOS, $FIELD_COSMIC};
    }
  }
}

sub get_tabix_nsfp {
  my $fn = $FLAGS{"tabix-dbnsfp"} || die;
  printf STDERR "dbNSFP: %s\n", $fn;
  my $tf = new TabixFile(
			 "-file" => $fn,
			 "-indel_wiggle_bases" => 0,
			);
  # database only contains SNVs so wiggle of 0 is OK
  return $tf;
}


sub hack_cosmic_tabix {
  my %options = @_;
  die unless @INPUT_SNV_INDEL_FILES;

  foreach my $infile (@INPUT_SNV_INDEL_FILES) {
    printf STDERR "processing %s...\n", $infile;
    my $df = new DelimitedFile(
			       "-file" => $infile,
			       "-headers" => 1
			      );
    my @rows;
    while (my $row = $df->get_hash()) {
      push @rows, $row;
    }

    add_batch_cosmic(
		     "-rows" => \@rows,
		    );
    foreach my $row (@rows) {
      printf "%s\n", join " ", @{$row}{$FIELD_VPOS, $FIELD_COSMIC};
    }
  }

}

sub hack_cosmic_tabix_germline {
  my %options = @_;
  die unless @INPUT_GL_FILES;

  foreach my $infile (@INPUT_GL_FILES) {
    printf STDERR "processing %s...\n", $infile;
    my $df = new DelimitedFile(
			       "-file" => $infile,
			       "-headers" => 1
			      );
    my @rows;
    while (my $row = $df->get_hash()) {
      push @rows, $row;
    }

    add_batch_cosmic(
		     "-rows" => \@rows,
		     "-germline" => 1,
		    );
    foreach my $row (@rows) {
      printf "%s\n", join " ", @{$row}{$FIELD_VPOS, $FIELD_COSMIC_TRUNCATING_FLAG};
    }
  }

}

sub hack_nsfp_tabix {
  my %options = @_;
  die unless @INPUT_SNV_INDEL_FILES;

  my $tf_nsfp = get_tabix_nsfp();
  my $sjpi = get_sjpi();

#  my @pf = qw(SIFT_pred Polyphen2_HDIV_pred MutationAssessor_pred);
  my @pf = qw(SIFT_pred);
  foreach my $infile (@INPUT_SNV_INDEL_FILES) {
    printf STDERR "processing %s...\n", $infile;
    my $df = new DelimitedFile(
			       "-file" => $infile,
			       "-headers" => 1
			      );
    my @rows;
    while (my $row = $df->get_hash()) {
      push @rows, $row;
    }

    add_batch_nsfp(
		   "-rows" => \@rows,
		   "-nsfp" => $tf_nsfp,
		   "-sjpi" => $sjpi,
		  );
    foreach my $row (@rows) {
      printf "%s\n", join "\t", @{$row}{($FIELD_VPOS, @pf)};
    }
  }

}

sub hack_clinvar_tabix {
  my %options = @_;
  die unless @INPUT_GL_FILES;

  foreach my $infile (@INPUT_GL_FILES) {
    printf STDERR "processing %s...\n", $infile;
    my $df = new DelimitedFile(
			       "-file" => $infile,
			       "-headers" => 1
			      );
    my @rows;
    while (my $row = $df->get_hash()) {
      push @rows, $row;
    }

    add_batch_clinvar(
		      "-rows" => \@rows,
		     );

    dump_die($rows[0]);
  }

}

sub hack_taylor {
  my %options = @_;
  die unless @INPUT_GL_FILES;

  foreach my $infile (@INPUT_GL_FILES) {
    printf STDERR "processing %s...\n", $infile;
    my $df = new DelimitedFile(
			       "-file" => $infile,
			       "-headers" => 1
			      );
    my @rows;
    while (my $row = $df->get_hash()) {
      push @rows, $row;
    }

    add_batch_taylor_hotspots(
			      "-rows" => \@rows,
			     );

    my $outfile = $infile . ".taylor.tab";
    my $rpt = $df->get_reporter("-file" => $outfile,
				"-extra" => [ "Reason" ]
			       );
    foreach my $row (@rows) {
      my $reason = "";
      $reason = $row->{$INTERNAL_FIELD_QUEUE_REASONS}->[0] if $row->{$INTERNAL_FIELD_QUEUE_REASONS};
      $row->{Reason} = $reason;
      $rpt->end_row($row);
    }
    $rpt->finish();

  }

}

sub nhlbi_hack_tabix {
  my $tabix_nhlbi = get_nhlbi_tabix();

  die unless @INPUT_GL_FILES;
  foreach my $glf (@INPUT_GL_FILES) {
    printf STDERR "processing %s...\n", $glf;
    my $df = new DelimitedFile(
			       "-file" => $glf,
			       "-headers" => 1
			      );
    my @rows;
    while (my $row = $df->get_hash()) {
      push @rows, $row;
    }

    if (1) {
      add_batch_nhlbi_frequency(
				"-rows" => \@rows,
				"-nhlbi" => $tabix_nhlbi
			       );
    } else {
      # single
      foreach my $row (@rows) {
	add_nhlbi_frequency_tabix($row, $tabix_nhlbi);
      }
    }

    # report
    foreach my $row (@rows) {
      printf STDERR "%s\n", $row->{NHLBI_frequency};
    }

  }
}


sub add_dbsnp {
  my ($row, $bi_dbsnp) = @_;
  my $hits = $bi_dbsnp->find(
			     "-sj" => $row,
			    );
  if (@{$hits}) {
    my @rs = join ",", sort map {$_->{name}} @{$hits};
    $row->{$FIELD_DBSNP} = join ",", @rs;
  }

  $row->{$FIELD_DBSNP} = "" unless $row->{dbSNP};

}

sub get_vm_committee {
  my (%options) = @_;
  my $vm = new VariantMatcher();
  $vm->enable_literal_aa_lookup(1);
  # indicate that the annotations we're adding here are in the same
  # format we'll be lookup up later (i.e. SJ annotations in this database
  # and during lookup)
  $vm->enable_literal_variant_lookup(1);
  # we need this so we can be very strict/specific about genomic indel lookups

  my $variants = 0;
  #  foreach my $fn ($medals_file) {

  my $csith_file;
  if ($CSITH_MEDALS_OVERRIDE_LEGACY_MEDALS) {
    my $genome = $FLAGS{genome} || die "-genome";
    my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
    $csith_file = $config_genome->{"CLINCLS_COMMITTEE_MEDALS_FILE_CSITH"};
    unless ($FLAGS{"committee-debug"}) {
      die "first medal file must be csith file" unless $COMMITTEE_GL_MEDALS[0] eq $csith_file;
    }
  }

  die "-gl-current-committee-medals" unless @COMMITTEE_GL_MEDALS;

  my %saw;

  my @usable_rows;

  if (my $f_new = $FLAGS{"hack-committee"}) {
    $COMMITTEE_GL_MEDALS[0] = $csith_file = $f_new;
    printf STDERR "DEBUG: hack committee medals file %s\n", $f_new;
  }

  foreach my $fn (@COMMITTEE_GL_MEDALS) {
    printf STDERR "loading committee calls (%s)...\n", $fn;
    my $df = new DelimitedFile(
			       "-file" => $fn,
			       "-headers" => 1,
			      );
    my $variants_this_file = 0;
    while (my $row = $df->get_hash()) {
#      dump_die($row, "debug", 1);
      die unless exists $row->{$FIELD_PANEL_DECISION};
      my $panel_call = $row->{$FIELD_PANEL_DECISION} || next;
      #    dump_die($row, "Debug", 1);
      $variants++;
      $variants_this_file++;
      my $chrom = $row->{Chr} || dump_die($row, "no Chr field");
      $chrom =~ s/^chr//;

      my $pos = get_sj_pos($row) || die;
      my $base_reference = $row->{ReferenceAllele} || die;
      my $base_variant = $row->{MutantAllele} || die;
      my $gene = $row->{$FIELD_GENE} || die;
      my $aa = $row->{$FIELD_AACHANGE} || die;
      dump_die($row, "missing field $FIELD_REFSEQ") unless exists $row->{$FIELD_REFSEQ};
      # needed during lookup for handshaking

#      printf STDERR "cvariant %s\n", join " ", $chrom, $pos, $base_reference, $base_variant, $gene, $aa, $panel_call;

      my $usable = 1;
      my $debug = 0;
      my $key = join ".", $chrom, $pos, $base_reference, $base_variant;
      printf STDERR "load %s: %s %s\n", $fn, $key, $panel_call if $debug;

      if ($CSITH_MEDALS_OVERRIDE_LEGACY_MEDALS) {
	if ($fn ne $csith_file) {
	  # in this mode, consider csith calls definitive:
	  # they are loaded first and will contain the newest decisions,
	  # which may sometimes change (e.g. 4.41747904.C.T was LB
	  # and is now U).  So, skip any variant from an older
	  # spreadsheet if csith has already loaded it.
	  if ($saw{$key}) {
	    printf STDERR "skip %s\n", $key if $debug;
	    $usable = 0;
	  }
	}
	$saw{$key} = 1;
      }

      if ($usable) {
	push @usable_rows, $row;
	if ($base_reference eq "-") {
	  $vm->add_insertion(
			     "-row" => $row,
			     "-sj" => 1,
			    );
	} elsif ($base_variant eq "-") {
	  $vm->add_deletion(
			    "-row" => $row,
			    "-sj" => 1,
			   );
	} else {
#	  foreach my $b ($base_reference, $base_variant) {
#	    dump_die($row, "FIX ME: dinucleotide? in SJ committee report") unless length($b) == 1;
#	  }
	  # this is OK now

#	  printf STDERR "adding committee SNV %s: %s\n", $key, $row->{$FIELD_PANEL_DECISION};

	  $vm->add_snv(
		       "-row" => $row,
		       "-sj" => 1
		      );
	}

	$vm->add_aa(
		    "-gene" => $gene,
		    "-aa" => $aa,
		    "-row" => $row
		   );
      }
    }
    printf STDERR "   loaded %d\n", $variants_this_file;
    die "no variants loaded from $fn" unless $variants_this_file;
  }
  printf STDERR "committee calls: %d\n", $variants;

  #
  #  sanity-check commitee calls:
  #
  my $db = $vm->db_snv();
  my $busted;
  foreach my $key (keys %{$db}) {
    my $rows = $db->{$key};
    if (@{$rows} > 1) {
      my %calls;
      foreach my $r (@{$rows}) {
	die "missing panel call" unless $r->{$FIELD_PANEL_DECISION};
	my $call = $r->{$FIELD_PANEL_DECISION};
	$calls{$call} = 1;
      }
      if (scalar keys %calls > 1) {
	$busted = 1;
	printf STDERR "ERROR: inconsistent committee call for %s: %s\n", $key, join ",", sort keys %calls;
      }
    }
  }
  if ($busted) {
    my $msg = "ERROR: inconsistent committee calls!!";
    if ($COMMITTEE_CONFLICT_WARN or $FLAGS{"committee-conflict-warn"}) {
      printf STDERR "%s\n", $msg;
    } else {
      die $msg;
    }
  }

  return $options{"-rows"} ? \@usable_rows : $vm;
}

sub tp53_germline_population_callback {
  my (%options) = @_;
  my $reasons = $options{"-reasons"} || die;
  my $hits = $options{"-hits"};
  my $label = $options{"-label"};

  if (0) {
    printf STDERR "hits=%d\n", scalar @{$hits};
    my $count = 0;
    foreach my $hit (@{$hits}) {
      $count++;
      dump_die($hit, "hit $count", 1);
    }
  }

  my %generations_found;
  foreach my $hit (@{$hits}) {
    my $gen = $hit->{Generation};
    unless ($gen) {
      dump_die($hit, "WARNING: no TP53 Generation field, assuming gen 1", 1);
      $gen = 1;
      # hack
    }
    $generations_found{$gen} = 1;
  }
  # there is a "Generations_analyzed", but count manually
  # just in case e.g. variant is not found in a later generation

  my $max_generation = (sort {$b <=> $a} keys %generations_found)[0];

  my %individuals = map {$_->{Individual_ID}, 1} @{$hits};

  push @{$reasons}, sprintf '%s_max_generation=%d', $label, $max_generation;
  push @{$reasons}, sprintf '%s_individuals=%d', $label, scalar keys %individuals;


}

sub get_nhlbi_tabix {
  return new NHLBI_tabix(
			 "-file" => $FLAGS{"tabix-nhlbi"} || die "-tabix-nhlbi"
			);
}

sub add_pcgp_somatic_population {
  my ($row, $vm) = @_;
  my $value = "";
  my $hits = $vm->find_snv("-sj" => $row);
  $hits = $vm->find_literal_variant("-sj" => $row) unless $hits;
  # for SJ-to-SJ lookup, require indel position annotations to match perfectly
  $value = scalar @{$hits} if $hits;
  $row->{$FIELD_PCGP_COUNT_SOMATIC} = $value;
  $row->{$FIELD_PCGP_DENOMINATOR_SOMATIC} = $vm->{$FIELD_PCGP_DENOMINATOR_SOMATIC} || die;
}

sub add_pcgp_germline_population {
  my ($row, $vm) = @_;
  my $value = "";
  my $hits = $vm->find_snv("-sj" => $row);
  # 3/2014: the germline list currently contains only SNVs.
  $value = scalar @{$hits} if $hits;
  $row->{$FIELD_PCGP_COUNT_GERMLINE} = $value;
  $row->{$FIELD_PCGP_DENOMINATOR_GERMLINE} = $vm->{$FIELD_PCGP_DENOMINATOR_GERMLINE} || die;
}

sub add_nhlbi_frequency {
  my ($row, $dbs_nhlbi) = @_;
  my $hits = $dbs_nhlbi->find("-sj" => $row);
  my $freq = "";
  if (@{$hits}) {
    # there may be multiple hits for different isoforms,
    # however underlying population data is the same
    my $ma = $row->{MutantAllele} || die;
    $freq = $dbs_nhlbi->get_parser()->get_genotype_frequency($hits->[0], $ma);
  }
  $row->{$FIELD_NHLBI_FREQ} = $freq;
}

sub add_nhlbi_frequency_tabix {
  my ($row, $tabix_nhlbi) = @_;
  my $hits = $tabix_nhlbi->find_one("-sj-post" => $row);
  my $freq = "";
  if (@{$hits}) {
    # there may be multiple hits for different isoforms,
    # however underlying population data is the same
    my $ma = $row->{MutantAllele} || die;
    $freq = $tabix_nhlbi->get_parser()->get_genotype_frequency($hits->[0], $ma);
  }
  $row->{$FIELD_NHLBI_FREQ} = $freq;
}

sub get_vm_pcgp_somatic {
  my $fn = $FLAGS{"gedi-recurrent"} || die "-gedi-recurrent";
  my $wanted = get_gedi_recurrent_restrict();
  my $gsm = new_gsm_lite();
  foreach my $g (keys %{$wanted}) {
    # safer than hash lookup: can adapt to gene symbol changes
    $gsm->add_gene("-gene" => $g);
  }

  my $df = new DelimitedFile("-file" => $fn,
			     "-headers" => 1,
			    );
  my $vm = get_new_vm();
  $vm->enable_literal_variant_lookup(1);

  my %samples;
  while (my $row = $df->get_hash()) {
    $samples{$row->{Sample} || die} = 1;

    #    dump_die($row);
    my $gene = $row->{Gene};
#    next unless $wanted->{$gene};
    next unless $gsm->find($gene);

    my $chr = $row->{Chr} || die;
    my $pos = $row->{Pos} || die;
    my $ref_allele = $row->{RefAllele} || die;
    my $var_allele = $row->{MutAllele} || die;

    my @opts = (
		"-row" => $row,
		"-reference" => $chr,
		"-reference-base" => $ref_allele,
		"-variant-base" => $var_allele
	       );

    if ($ref_allele eq "-") {
      # insertion
      $vm->add_insertion(
			 @opts,
			 "-start" => $pos,
			 "-end" => $pos,
			 "-count" => length($var_allele)
			);
    } elsif ($var_allele eq "-") {
      # deletion
      $vm->add_deletion(
			@opts,
			"-start" => $pos,
			"-end" => ($pos + length($ref_allele) - 1),
		       );
    } elsif ($ref_allele =~ /\-/ or $var_allele =~ /\-/) {
      # complex indel: literal lookup only!
      #      printf STDERR "complex indel: ref=%s var=%s\n", $ref_allele, $var_allele;
      $vm->add_literal_variant(
			       @opts,
			       "-start" => $pos,
			      );
    } else {
      # substitution
      if (length($ref_allele) != length($var_allele)) {
	dump_die($row, "ERROR: possible ref/var formatting problem ref=$ref_allele var=$var_allele", 1);
	$vm->add_literal_variant(
				 @opts,
				 "-start" => $pos,
				 #				 "-base-number" => $pos,
				);
      } else {
	$vm->add_snv(
		     @opts,
		     "-reference" => $chr,
		     "-base-number" => $pos,
		    );
      }
    }
  }

  $vm->{$FIELD_PCGP_DENOMINATOR_SOMATIC} = scalar keys %samples;
  return $vm;
}

sub get_vm_pcgp_germline {
  #
  # PCGP recurrent germline SNVs
  #
  my $fn = $FLAGS{"gl-pcgp-population"} || die "-gl-pcgp-population";
  printf STDERR "parsing PCGP germline (%s)...", $fn;
  my $df = new DelimitedFile(
			     "-file" => $fn,
			     "-headers" => 1,
			    );
  my $vm = new VariantMatcher();
  my %samples;
  while (my $row = $df->get_hash()) {
    my $sample = $row->{Sample} || die;
    $samples{$sample}++;
    if ($row->{ReferenceAllele} eq "-" or
	$row->{MutantAllele} eq "-") {
      # not seeing any yet
      die "PCGP germline indel, not implemented";
    } else {
      $vm->add_snv("-row" => $row,
		   "-sj" => 1);

      my $gene = $row->{GeneName} || die;
      my $aa = $row->{AAChange} || die;
      $vm->add_aa(
		  "-gene" => $gene,
		  "-aa" => $aa,
		  "-row" => $row
		 );
    }
  }
  $vm->{$FIELD_PCGP_DENOMINATOR_GERMLINE} = scalar keys %samples;
  # hacktacular

  print STDERR "done\n";
  return $vm;
}

sub add_foldx {
  my ($row) = @_;
  $row->{$FIELD_FOLDX} = "" unless $row->{$FIELD_FOLDX};
  # if already provided (e.g. manually by Gang), leave it
}

sub ts_truncating_medal {
  my ($truncation_class) = @_;
  my $medal;
  if (
      $truncation_class == $CLASS_TRUNCATION or
      $truncation_class == $CLASS_SPLICE
     ) {
    $medal = CLASS_GOLD;
  } elsif ($truncation_class == $CLASS_SPLICE_INTRON) {
    $medal = CLASS_SILVER;
  } else {
    die "invalid truncation class $truncation_class";
  }
  return $medal;
}

sub add_cosmic {
  my ($row, $vm) = @_;

  my $hits = $vm->find_snv("-sj" => $row);
  $hits = indel_search($vm, $row) unless $hits;
  #  my $value = "-";
  my $value = "";
  if ($hits and @{$hits}) {
    #    dump_die($row, "hey now cosmic hit", 1);
    #    dump_die($hits->[0], "COSMIC:");
    $value = "Present";
  }
  $row->{$FIELD_COSMIC} = $value;
}

sub parse_nhgri {
  my (%options) = @_;
  my $gene = $options{"-gene"} || die "-gene";

  my $REFSEQ_CHECK = 0;
  my $vm = get_new_vm($REFSEQ_CHECK);
  my $flag = sprintf 'nhgri-%s', lc($gene);
  my $infile = $FLAGS{$flag} || die "-" . $flag;
  die unless -s $infile;

  printf STDERR "parsing NHGRI (%s)...\n", $infile;

  my $chrom;
  my $nm;
  if ($gene eq "BRCA1") {
    $chrom = 17;
    $nm = $NHGRI_BRCA1_NM;
  } elsif ($gene eq "BRCA2") {
    $chrom = 13;
    $nm = $NHGRI_BRCA2_NM;
  } else {
    die;
  }

  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			    );
  my $nsp = new NucleotideSubstitutionParser();
  my %imported;
  while (my $row = $df->get_hash()) {
    #    dump_die($row, "debug", 1);

    $row->{$FIELD_REFSEQ} = $nm;
    # record for transcript handshaking during later lookup

    #
    #  import variant genomically:
    #
    my $gpos = $row->{"HGVS Genomic (hg19)"};
    if (exists $row->{sj_pos}) {
      # refactor + liftover for GRCh38+
      my $pos = $row->{sj_pos};
      if ($pos) {
	my $ra = $row->{sj_ref_allele} || dump_die($row, "no ref allele");
	my $va = $row->{sj_alt_allele} || dump_die($row, "no variant allele");
	my $v = new Variant();
#	printf STDERR "import NHGRI %s\n", join ".", $chrom, $pos, $ra, $va;
	$v->import_generic(
		     "-reference-name" => $chrom,
		     "-base-number" => $pos,
		     "-reference-allele" => $ra,
		     "-variant-allele" => $va,
		    );
	$vm->add_variant($v, "-row" => $row);
      }
    } elsif ($gpos ne "-") {
      # hg19
      $nsp->parse($gpos) || die "can't parse $gpos";
      my $start = $nsp->start;
      my $end = $nsp->end;
      if ($nsp->is_substitution()) {
	my $ref_base = $nsp->reference_sequence();
	my $var_base = $nsp->variant_sequence();
	my $ok = $vm->add_snv(
			      "-row" => $row,
			      "-reference" => $chrom,
			      "-base-number" => $start,
			      "-reference-base" => $ref_base,
			      "-variant-base" => $var_base
			     );
	die unless $ok;
      } elsif ($nsp->is_deletion()) {
	$vm->add_deletion(
			  "-row" => $row,
			  "-reference" => $chrom,
			  "-start" => $start,
			  "-end" => $end
			 );
      } elsif ($nsp->is_insertion()) {
	# simple insertion
	$vm->add_insertion(
			   "-reference" => $chrom,
			   "-start" => $start,
			   "-end" => $end,
			   "-row" => $row,
			   "-count" => ($nsp->event_length() || die "no ins length"),
			  );
      } elsif ($nsp->is_complex_indel()) {
	$vm->add_deletion(
			  "-reference" => $chrom,
			  "-start" => $start,
			  "-end" => $end,
			  "-row" => $row
			 );
	$vm->add_insertion(
			   "-reference" => $chrom,
			   "-start" => $start,
			   "-end" => $end,
			   "-row" => $row,
			   "-count" => ($nsp->complex_insertion_length() || die)
			  );
      } else {
	die sprintf "FIX ME: gpos $gpos unhandled";
      }
    }

    #
    #  import variant by AA:
    #
    my $aa = $row->{"HGVS Protein"} || die;
    $aa = trim_flanking_whitespace($aa);
    # ugh, e.g. "p.Gln148del "

    if ($aa ne "-") {
      if ($vm->add_aa(
		      "-gene" => $gene,
		      "-aa" => $aa,
		      "-row" => $row
		     )) {

      }
    }
  }
  return $vm;
}

sub nhgri_info_callback {
  my (%options) = @_;
  my $reasons = $options{"-reasons"} || die;
  my $hits = $options{"-hits"} || die;
  my $label = $options{"-label"} || die;

  my @fields = (
		"Clinically Important",
		"Depositor",
		#		"Number Reported",
		# initially thought this was # of patients, but field is sparsely
		# populated, row count is a better proxy for that
		"Ethnicity",
	       );

  my %synonyms = (
		  "Western Europe" => "Western European",
		  "Western-European" => "Western European",

		  "Latin-American/Caribbean" => "Latin American/Caribbean",

		  "Central/Eastern Europe" =>  "Central/Eastern European",
		  "Central/Eastern-European" =>  "Central/Eastern European",
		 );

  my %null = map {$_, 1} (
			  "-",
			  "None Specified",
			  "None-Specified",
			  "Not Specified",
			 );

  push @{$reasons}, sprintf '%s_patients=%d', $label, scalar @{$hits};

  foreach my $field (@fields) {
    # summarize hit fields of interest
    my %values;
    foreach my $hit (@{$hits}) {
      #      dump_die($hit, "debug NHGRI hit", 1);

      my $thing = $hit->{$field};
      $thing =~ s/^\"//;
      $thing =~ s/\"$//;

      $thing = lc($thing) if $field =~ /clinically/i;

      my @things = split /,\s*/, $thing;
      # split lists, e.g.:
      # Native American, Central/Eastern European

      foreach my $v (@things) {
	if (my $to = $synonyms{$v}) {
	  $v = $to;
	}
	next if $null{$v};
	$values{$v} = 1;
      }

    }
    my $reason_label = sprintf '%s_%s', $label, $field;
    $reason_label =~ s/\s+/_/g;
    push @{$reasons}, sprintf '%s=%s', $reason_label, join ",", sort keys %values;
  }

}

sub hack_reannotate {
  my ($infile) = @_;
  my $dbnsfp = get_nsfp();
  my $sjpi = get_sjpi();

  my $outfile = sprintf "%s.reannotated.tab", basename($infile);

  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			    );
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   qw(
					       SJ_mRNA_acc
					       SJ_mRNA_acc_difference
					       SJ_AA_NSFP
					       SJ_AA_difference
					    )
					  ]
			     );

  while (my $row = $df->get_hash()) {
    my $nm_preferred = $sjpi->get_preferred_isoform($row->{$FIELD_GENE} || die) || die;
    my $nm_orig = $row->{$FIELD_REFSEQ} || die;
    $row->{SJ_mRNA_acc} = $nm_preferred;
    my $rna_difference = $nm_orig eq $nm_preferred ? 0 : 1;
    $row->{SJ_mRNA_acc_difference} = $rna_difference;
    my $aa_prev = $row->{$FIELD_AACHANGE} || die;

    my $hits = $dbnsfp->find(
			     "-chr" => ($row->{Chr} || die),
			     "-pos" => (get_sj_pos($row) || die),
			     "-reference-base" => ($row->{ReferenceAllele} || die),
			     "-variant-base" => ($row->{MutantAllele} || die),
			     "-nm" => $nm_preferred
			    );
    my $aa;
    if ($hits and @{$hits}) {
      $aa = $dbnsfp->get_uniprot_aachange(
					  "-row" => $hits->[0],
					  "-strict" => 1,
					  "-relax-if-nonsense" => $aa_prev
					 );
    }

    $aa = "" unless $aa;
    if (not($aa) and not($rna_difference)) {
      # nothing from dbNSFP
      # (stops may not have uniprot accessions, etc.)
      # defer to raw call if mRNA hasn't changed
      $aa = $aa_prev;
    }

    $row->{SJ_AA_NSFP} = $aa;
    $row->{SJ_AA_difference} = $aa eq $aa_prev ? 0 : 1;
    $rpt->end_row($row);
  }
  $rpt->finish;

}

sub umd_info_callback {
  my (%options) = @_;
  my $reasons = $options{"-reasons"} || die;
  my $hits = $options{"-hits"};
  my $label = $options{"-label"};

  my %sig;
  foreach my $hit (@{$hits}) {
    my $sig = $hit->{Biological_significance_cleaned} || next;
    $sig{$sig} = 1;
  }
  my @sig = sort keys %sig;

  push @{$reasons}, sprintf "%s_significance=%s",
    $label, join ",", @sig if @sig;
}

sub get_vm_generic_aa {
  #
  # import a flatfile using standardized column names, AA annotations only
  #
  my (%options) = @_;
  my $label = $options{"-label"} || die "-label";
  my $flag = $options{"-flag"} || die "-flag";

  my $infile = $FLAGS{$flag} || die "-" . $flag;
  printf STDERR "parsing %s (%s)...\n", $label, $infile;

  my $vm = get_new_vm();
  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			    );
  my %count;
  while (my $row = $df->get_hash()) {
    die unless exists $row->{$FIELD_AACHANGE};
    my $aa = $row->{$FIELD_AACHANGE};
    if ($aa =~ /\w/) {
      if ($vm->add_aa(
		      "-gene" => ($row->{$FIELD_GENE} || die),
		      "-aa" => $aa,
		      "-row" => $row
		     )) {
	$count{ok}++;
      } else {
	printf STDERR "can't parse %s AA %s\n", $label, $aa if $VERBOSE_AA_PARSE_WARNING;
	$count{unparsable}++;
      }
    } else {
      $count{blank}++;
    }
  }
  printf STDERR "%s summary:\n", $label;
  foreach (sort keys %count) {
    printf STDERR "  %s: %d\n", $_, $count{$_};
  }

  return $vm;
}

sub lovd_germline_callback {
  my (%options) = @_;
  my $reasons = $options{"-reasons"} || die;
  my $hits = $options{"-hits"};
  my $label = $options{"-label"};

  my %sig;
  foreach my $hit (@{$hits}) {
    my $sig = $hit->{Pathogenicity_concluded_cleaned} || next;
    $sig{$sig} = 1;
  }
  my @sig = sort keys %sig;

  push @{$reasons}, sprintf "%s_significance=%s",
    $label, join ",", @sig if @sig;
}

sub nhgri_medal_callback {
  my ($hits, $medal, $perfect_match) = @_;
  my $medal_orig = $medal;

  my $YES = 1;
  my $NO = 0;
  my $UNKNOWN = -1;

  my %saw;
  foreach my $hit (@{$hits}) {
    #    dump_die($hit, "NHGRI match debug", 1);
    my $code;
    my $raw = lc($hit->{"Clinically Important"});
    if ($raw eq "yes") {
      $code = $YES;
    } elsif ($raw eq "no") {
      $code = $NO;
    } elsif ($raw eq "unknown" or $raw eq "-") {
      $code = $UNKNOWN;
    } else {
      die "can't parse importance code $raw";
    }
    $saw{$code} = 1;
  }

  if ($saw{$YES}) {
    $medal = $perfect_match ? CLASS_GOLD : CLASS_SILVER;
  } elsif ($saw{$UNKNOWN}) {
    $medal = $perfect_match ? CLASS_SILVER : CLASS_BRONZE;
  } elsif ($saw{$NO}) {
    $medal = CLASS_BRONZE;
  } else {
    die;
  }

  printf STDERR "NHGRI hit: perfect=%d codes=%s medal=%s\n",
    ($perfect_match ? 1 : 0),
      join(",", sort {$a <=> $b} keys %saw),
	$medal;

  return $medal;
}

sub show_isoforms_used {
  my ($meta) = @_;

  my $df;
  my %gene2nm;

  foreach my $fn (qw(
		      /nfs_exports/genomes/1/projects/ClinicalSeq/germline/LOVD_APC/export_APC_chromium.liacs.nl_APC131203.tab
		      /nfs_exports/genomes/1/projects/ClinicalSeq/germline/LOVD_MSH2/export_MSH2_chromium.liacs.nl_MSH2_131203.tab
		   )
		 ) {
    # HACK: not yet in meta report
    my @f = split /_/, basename($fn);
    my $db = sprintf "LOVD_%s", $f[1];

    $df = new DelimitedFile(
			    "-file" => $fn,
			    "-headers" => 1,
			   );
    my $row = $df->get_hash();
    my $gene = $row->{GeneName} || die;
    my $nm = $row->{mRNA_acc} || die;
    $nm =~ s/\.\d+$//;
    $gene2nm{$gene}{$nm}{$db} = 1;
  }

  #
  #  meta db:
  #
  $df = new DelimitedFile(
			  "-file" => $meta,
			  "-headers" => 1,
			 );
  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    my $gene = $row->{Gene} || die;
    my $db = $row->{database} || die;
    my $nm = $row->{mRNA_acc};
    if ($nm) {
      $nm =~ s/\.\d+$//;
      $gene2nm{$gene}{$nm}{$db} = 1;
    }
  }

  foreach my $gene (keys %gene2nm) {
    #
    #  also include current SJ preferred isoform
    #
    my $sjpi_old = new SJPreferredIsoform("-file" =>
					  "/nfs_exports/apps/gnu-apps/NextGen/nextgensupport2/NCBI/RefSeq/gene_transcript_matrix.withNM.mod"
					 );

    my $sjpi_new = new SJPreferredIsoform("-file" =>
#					  "/nfs_exports/genomes/1/PCGP/BucketIntermediate/PN/gene_transcript_matrix.withNM.forJinghui",
#					  "/nfs_exports/apps/gnu-apps/NextGen/nextgensupport2/NCBI/RefSeq/gene_transcript_matrix.withNM.mod"
					  "sjpi.tab",
					  # temp
					 );
    my $nm_old = $sjpi_old->get_preferred_isoform($gene) || die;
    my $nm_new = $sjpi_new->get_preferred_isoform($gene) || die;
    $gene2nm{$gene}{$nm_old}{SJ_old} = 1;
    $gene2nm{$gene}{$nm_new}{SJ_new} = 1;
  }

  my $rpt = new Reporter(
			 "-file" => "isoforms_used.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   GeneName
					   info
					   consistent
					)
				      ],
			 "-auto_qc" => 1,
			);

  foreach my $gene (sort keys %gene2nm) {
    my @things;
    foreach my $nm (sort keys %{$gene2nm{$gene}}) {
      my $dbs = $gene2nm{$gene}{$nm};
      push @things, sprintf '%s=%s', $nm, join ",", sort keys %{$dbs};
    }
    my %r;
    $r{GeneName} = $gene;
    $r{info} = join " ", @things;
    $r{consistent} = @things == 1 ? 1 : 0;
    $rpt->end_row(\%r);
  }

  print STDERR "NOTE: NHGRI BRCA1 is a placeholder, not 100%% compatible!\n";
  $rpt->finish();
}


sub committee_hack {
  my $vm_committee = get_vm_committee();
  die "parse done.\n";
  foreach my $infile (@INPUT_GL_FILES) {
    my $df = new DelimitedFile(
			       "-file" => $infile,
			       "-headers" => 1,
			      );
    while (my $row = $df->get_hash()) {
      my $medal = "lead";
      my @reasons;
      my $is_novel_these_db = 0;
      my $medal_call_protected_from_nhlbi = 0;

      dump_die($row, "start lookup", 1);

      my ($hits, $perfect) = assign_aa_match(
					     "-vm" => $vm_committee,
					     "-row" => $row,
					     "-transcript-handshake" => $FIELD_REFSEQ,
					     "-current-medal" => $medal,
					     "-reasons" => \@reasons,
					     "-label" => "committee_reviewed",
					     "-novel-flag" => \$is_novel_these_db,
					     "-no-medal" => 1,
					     "-protect" => \$medal_call_protected_from_nhlbi,
					     "-genomic-lookup" => 1,
					     "-return-hits-and-perfect" => 1
					    );

      printf STDERR "Reasons: %s\n", join ",", @reasons;
      printf STDERR "perfect=%d\n", $perfect || 0;

      foreach my $hit (@{$hits}) {
	dump_die($hit, "reported hit:", 1);
      }
      print STDERR "\n\n";

    }
  }

}

sub grf_hack {
  my $gl_gold_ranges_file = $FLAGS{"gl-gold-cancer-ranges"} || die "-gl-gold-cancer-ranges";

  my $grf_gold = new GenomicRangeFinder();
  my $gl_gold_genes = {};

  open(GLTMP, $gl_gold_ranges_file) || die;
  while (<GLTMP>) {
    chomp;
    my @f = split /\t/, $_;
    die unless @f == 2;
    my ($gene, $loc) = @f;
    $gl_gold_genes->{$gene} = 1
      unless $gene eq "TP53" and $FLAGS{"gl-grf-debug"};

    $grf_gold->add(
		   "-range" => $loc,
		   "-value" => {
				"gene" => $gene,
				"location" => $loc
			       },
		  );
  }

  my $hits = $grf_gold->find(
  			     "-chrom" => "chr5",
  			     "-base" => 653248
   			    );
  die scalar @{$hits};

}

sub umd_hack {
  my $vm_umd = get_vm_generic_aa(
				 "-label" => "UMD",
				 "-flag" => "gl-umd-flatfile",
				);
  single_vm_hack(
		 "-vm" => $vm_umd,
		 "-label" => "UMD"
		);
}

sub single_vm_hack {
  my (%options) = @_;
  my $vm = $options{"-vm"} || die;
  my $label = $options{"-label"} || die;
  foreach my $glf (@INPUT_GL_FILES) {
    printf STDERR "processing %s...\n", $glf;
    my $df = new DelimitedFile(
			       "-file" => $glf,
			       "-headers" => 1
			      );

    while (my $row = $df->get_hash()) {
      #      dump_die($row, "start variant lookup", 1);

      my $variant_class = lc($row->{Class} || die "no Class field");
      my $is_silent_intron_utr = 0;
      if ($variant_class eq "silent" or
	  $variant_class eq "intron" or
	  $variant_class eq "utr_5" or
	  $variant_class eq "utr_3") {
	$is_silent_intron_utr = 1;
      }

      my @reasons;
      my $medal;

      my $is_novel_these_db;

      $medal = assign_aa_match(
			       "-vm" => $vm,
			       "-row" => $row,
			       "-protect" => 0,
			       "-is-silent" => $is_silent_intron_utr,
			       "-current-medal" => $medal,
			       "-transcript-handshake" => $FIELD_REFSEQ,
			       "-reasons" => \@reasons,
			       "-label" => $label,
			       "-novel-flag" => \$is_novel_these_db,
			       "-medal-perfect" => CLASS_SILVER,
			       "-medal-codon" => CLASS_SILVER,
			       "-medal-site" => CLASS_SILVER,
			       "-genomic-lookup" => 1,
			      );
      die $medal;

      @reasons = "" unless @reasons;

      printf "%s\n",
	join ",", $medal, @reasons;
    }
  }
}

sub pcgp_somatic_hack {
  my $vm = get_vm_pcgp_somatic() || die;
  die;
}

sub column_qc_hack {
  die unless @INPUT_SNV_INDEL_FILES;
  foreach my $f (@INPUT_SNV_INDEL_FILES) {
    line_delimiter_qc("-file" => $f, "-delimiter" => "\t");
  }
  die "done";
}

sub get_sj_pos {
  my ($row) = @_;
  die unless $row;
  my $pos = $row->{$FIELD_VPOS} || $row->{$FIELD_VPOS_ALT};
  dump_die($row, "can't find base number") unless $pos;
  return $pos;
}

sub file_sanity_check {
  my ($list) = @_;
  my $fail;
  foreach my $fn_raw (@{$list}) {
    my ($fn, $outfile) = get_infile_and_outfile($fn_raw);
    next unless $fn =~ /\w/;
    # in case blank line in input list
    if (not(-e $fn)) {
      printf STDERR "ERROR: %s does not exist\n", $fn;
      $fail = 1;
    } elsif (not(-f $fn)) {
      printf STDERR "ERROR: %s is not a plain file\n", $fn;
      $fail = 1;
      #    } elsif (not(-r $fn)) {
      #      printf STDERR "ERROR: %s is not readable\n", $fn;
      #      $fail = 1;
      #
      # this does NOT seem to work on some filesystems (?),
      # file supposedly not readable but it actually is!
      #
    } elsif (-s $fn) {
      line_delimiter_qc("-file" => $fn, "-delimiter" => "\t");
    }
  }
  die "fatal error: invalid or inaccessible input file" if $fail;
}

sub get_sv_genes {
  my $svc = $FLAGS{"sv"} || die "specify -sv";
  my $svc_manual = $FLAGS{"sv-manual"} || die "-sv-manual";
  my $sv_silver = $FLAGS{"sv-silver-genes-fb"} || die "-sv-silver-genes-fb";

  printf STDERR "SV configs: main=%s manual=%s silver=%s\n", $svc, $svc_manual, $sv_silver;

  my %genes;

  my $lines = read_simple_file($svc);
  # primary SVs
  foreach my $line (@{$lines}) {
    my ($g1, $g2, $desc) = split /\t/, $line;
    foreach my $g ($g1, $g2) {
      $genes{$g} = 1 if $g;
    }
  }

  # manually-curated SVs:
  my $df = new DelimitedFile(
			     "-file" => $svc_manual,
			     "-headers" => 1,
			    );
  while (my $row = $df->get_hash()) {
    my ($gene1, $gene2, $submitter) = @{$row}{qw(pair1 pair2 contact)};
    foreach my $g ($gene1, $gene2) {
      $genes{$g} = 1 if $g;
    }
  }

  # other SV genes which should simply receive a silver:
  $df = new DelimitedFile(
			  "-file" => $sv_silver,
			  "-headers" => 1,
			 );
  while (my $row = $df->get_hash()) {
    my $g = $row->{gene} || die;
    $genes{$g} = 1;
  }

  return \%genes;
}

sub sv_gene_check_old {

  die "old version";
  my ($rf) = @_;

  my $hgnc = new HGNCParser("-file" => $FLAGS{hgnc} || die "-hgnc");

  my $gss = new GeneSymbolStandardizer("-filename" => $FLAGS{"gene-info"} || die "-gene-info");

  my $gsm = new GeneSymbolMapper("-hgnc_file" => $FLAGS{hgnc} || die "-hgnc");

  # gene symbols reported by FusionBuilder:
  open(IN, $rf) || die;
  while (<IN>) {
    chomp;
    my @f = split /\t/, $_;
    my $gene = clean_sharp_gene_symbol($f[0]);
    my $nm = $f[1];
    $gsm->add_gene(
		   "-gene" => $gene,
		   "-refseq" => $nm
		  );
  }

  # gene symbols used by SV annotations
  my %sv_genes;
  foreach my $flag (qw(sv sv-manual)) {
    my $f = $FLAGS{$flag} || die $flag;
    open(IN, $f) || die;
    while (<IN>) {
      chomp;
      my @f = split /\t/, $_;
      foreach my $g ($f[0], $f[1]) {
	next unless $g =~ /\w/;
	next if $g eq "pair1" or $g eq "pair2";
	$sv_genes{$g} = 1;
      }
    }
  }

  #
  # compare:
  #
  my %repair;
  $repair{"6-Sep"} = "SEPT6";
  # excel damage

  foreach my $gene (sort keys %sv_genes) {
    if ($gsm->contains($gene)) {
      printf STDERR "hit for %s (ok)\n", $gene;
    } else {
      my $alt;

      my $alts = $hgnc->find("-symbol" => $gene);
      my @alts;
      if ($alts) {
	@alts = map {$_->{"Approved Symbol"}} @{$alts};
	if (@{$alts} == 1) {
	  # unambiguous
	  #	  dump_die($alts->[0], "hit for $gene");
	  ($alt) = @alts;
	} else {
	  printf STDERR "disqualifying map for %s, %d matches\n", $gene, scalar @alts;
	}
      }
      my $alt_found = $alt ? ($gsm->contains($alt) || 0) : 0;

      my $alt2 = "";
      my @alt2;
      unless ($alt_found) {
	# only populate if HUGO lookup doesn't work
	my $syns = $gss->find($gene);
	if ($syns and ref $syns) {
	  if (@{$syns} == 1) {
	    $alt2 = $syns->[0];
	  }
	  @alt2 = @{$syns};
	}
      }
      my $alt2_found = $alt2 ? ($gsm->contains($alt2) || 0) : 0;

      printf STDERR "miss for %s, approved=%d HUGO=%s HUGO_usable=%s EG=%s EG_usable=%d\n",
	$gene,
	  ($hgnc->is_approved($gene) ? 1 : 0),
	    (join(",", @alts)), $alt_found,
	      join(",", @alt2), $alt2_found;

      $repair{$gene} = $alt if $alt_found;
      #      printf STDERR "alt2 repair: $gene => $alt2\n" if $alt2_found;
      $repair{$gene} = $alt2 if $alt2_found;
    }
  }

  my $infile = $FLAGS{sv} || die;
  my $outfile = sprintf '%s.repaired.txt', basename($infile);
  open(IN, $infile) || die;
  open(OUT, ">" . $outfile) || die;
  while (<IN>) {
    chomp;
    my $raw = $_;
    my @f = split /\t/, $_;
    my ($g1, $g2, $desc) = @f;
    my $repaired = 0;
    foreach my $g ($g1, $g2) {
      if (my $to = $repair{$g}) {
	$repaired++;
	$g = $to;
      }
    }
    my $new = sprintf '%s', join "\t", $g1, $g2, $desc;
    printf OUT "%s\n", $new;
  }
  close IN;
  close OUT;
}


sub clean_sharp_gene_symbol {
  # in Michael Rusch's "sharp" versions of refFlat,
  # ambiguous loci assigned "_loc" suffix, e.g. DUX4_locF.
  # These will also be reported by FusionBuilder with this suffix,
  # so need to remove it to get the raw gene symbol.
  my ($gene) = @_;
  if ($gene =~ /_loc/) {
    my $before = $gene;
    $gene =~ s/_loc.*$//;
#    printf STDERR "stripping %s => %s\n", $before, $gene;
  }
  return $gene;
}

sub get_gsm_fb {
  # gene symbols reported by FusionBuilder
  # config REFSEQ_NM_REFFLAT, e.g.
  # clinical:
  #   - pre-CG2.0.0:
  #     /clingen/dev/tartan/runs/genome/ReJNnZE8/output/refFlat-sharp.txt
  #   - CG2.0.0 test version (contains updates):
  #     /clingen/test/tartan/index/reference/Homo_sapiens/GRCh37-lite/mRNA/RefSeq/refFlat-sharp.txt
  # research:
  # /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/mRNA/RefSeq/refFlat-sharp.txt
  # NOTE: CLINICAL VERSION NEWER THAN RESEARCH
  my $rf = $FLAGS{"refflat-fb"};
  unless ($rf) {
    printf STDERR "-refflat-fb not specified, using config refflat\n";
    my $genome = $FLAGS{genome} || die "-genome";
    my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
    $rf = $config_genome->{REFSEQ_NM_REFFLAT} || die "no config REFSEQ_NM_REFFLAT";
  }

  my $gsm = new_gsm();

  open(IN, $rf) || die;
  while (<IN>) {
    chomp;
    my @f = split /\t/, $_;
    my $gene = clean_sharp_gene_symbol($f[0]);
    if ($gene) {
      my $nm = $f[1];
      $gsm->add_gene(
		     "-gene" => $gene,
		     "-refseq" => $nm
		    );
    } else {
      printf STDERR "no gene symbol in refFlat entry %s\n", $_;
      # e.g. _locPar NR_001526-locPar
    }
  }
  close IN;

  return $gsm;
}

sub sv_gene_check {
  my $gsm = get_gsm_fb();
  # FusionBuilder genes

  my $sv_genes = get_sv_genes();
  # gene symbols used by SV annotations

  #
  # compare:
  #
  my %repair;
  $repair{"6-Sep"} = "SEPT6";
  # excel damage

  foreach my $gene (sort keys %{$sv_genes}) {
    if ($gsm->contains($gene)) {
      printf STDERR "hit for %s (ok)\n", $gene;
    } else {
      my $result = $gsm->resolve("-symbol" => $gene);
      printf STDERR "miss for %s, map=%s\n", $gene, ($result || "n/a");
      $repair{$gene} = $result if $result;
    }
  }

  my $infile = $FLAGS{sv} || die;
  my $outfile = sprintf '%s.repaired.txt', basename($infile);
  open(IN, $infile) || die;
  open(OUT, ">" . $outfile) || die;

  print STDERR "FIX ME: add sanity check?: confirm each symbol in raw list\n";
  printf STDERR "TO DO: report whether any entries changed\n";
  printf STDERR "TO DO: process both -sv and -sv-manual\n";

  while (<IN>) {
    chomp;
    my $raw = $_;
    my @f = split /\t/, $_;
    my ($g1, $g2, $desc) = @f;
    my $repaired = 0;
    foreach my $g ($g1, $g2) {
      if (my $to = $repair{$g}) {
	$repaired++;
	$g = $to;
      }
    }
    my $new = sprintf '%s', join "\t", $g1, $g2, $desc;
    printf OUT "%s\n", $new;
  }
  close IN;
  close OUT;
}

sub cnv_sv_gene_patch {
  #
  # CNV classification also uses SV genes.
  # patch the SV gene list into GENE_EXON_REGION gene name space.
  #

  printf STDERR "SV config file: %s\n", $FLAGS{sv};

  printf "**** is this the latest/finalized SV file? [y/n]: ";
  my $resp = <STDIN>;
  chomp $resp;
  die "quitting\n" unless (lc($resp) eq "y");

  #    $FLAGS{sv} = "SVCheck.txt.repaired.txt";
  # hack: use the latest corrected version

  my $gsm = get_gsm_ger();

  # gene symbols used by SV annotations
  my $sv_genes = get_sv_genes();

  #
  # compare:
  #
  my $wf = new WorkingFile("SV_genes_patched_to_GENE_EXON_REGION.txt");
  my $fh = $wf->output_filehandle();

  foreach my $gene (sort keys %{$sv_genes}) {
    my $updated;
    unless ($gsm->contains($gene)) {
      my $result = $gsm->resolve("-symbol" => $gene);
      printf STDERR "miss for %s, map=%s\n", $gene, ($result || "n/a");
      $updated = $result if $result;
    }
    printf $fh "%s\n", $updated || $gene;
  }

  $wf->finish();
}

sub dump_cnv_genes {
  my $cnv_genes = get_cnv_genes();
  foreach (sort keys %{$cnv_genes}) {
    printf "%s\n", $_;
  }
}

sub dump_sv_genes {
  my $sv_genes = get_sv_genes();
  foreach (sort keys %{$sv_genes}) {
    printf "%s\n", $_;
  }
}

sub get_cnv_genes {
  my $cnv = $FLAGS{cnv} || die "-cnv";
  my $cnv_annot = $FLAGS{"cnv-annotations"} || die "-cnv-annotations";
  printf STDERR "files: %s\n", join " ", $cnv, $cnv_annot;

  open(IN, $cnv) || die;
  my %genes;
  while (<IN>) {
    chomp;
    my @f = split /\t/, $_;
    $genes{$f[0]} = 1;
  }

  my $df = new DelimitedFile(
			     "-file" => $cnv_annot,
			     "-headers" => 1,
			    );
  while (my $row = $df->get_hash()) {
    $genes{$row->{"Gene Symbol"}} = 1;
  }

  return \%genes;
}

sub verify_cnv_genes {
  my $annot = $FLAGS{"gene-exon-region-dir"} || die "specify -gene-exon-region-dir";
  # need genomically-sorted genes in interval for analysis
  my $ga = new GeneAnnotation(
			      "-style" => "gene_exon_region",
			      "-gene_exon_region_dir" => $annot,
			      "-ignore_non_coding" => 0
			     );

  my $cnv_genes = get_cnv_genes();
  my $sv_genes = get_sv_genes();

  foreach my $ref (
		   [ $cnv_genes, "CNV" ],
		   [ $sv_genes, "SV" ]
		  ) {
    my ($hash, $label) = @{$ref};

    foreach my $gene (sort keys %{$hash}) {
      my $valid = $ga->is_valid_gene("-gene" => $gene) || 0;
      printf "%s\n", join " ", $label, $gene, $valid;
    }
  }
}

sub new_gsm {

  unless ($FLAGS{hgnc} and $FLAGS{"gene-info"}) {
    my $config_species = TdtConfig::readConfig('species', "Homo_sapiens") || die "can't find hs species config";
    unless ($FLAGS{hgnc}) {
      $FLAGS{hgnc} = $config_species->{HGNC} || die "no HGNC config";
    }

    unless ($FLAGS{"gene-info"}) {
      $FLAGS{"gene-info"} = $config_species->{ENTREZ_GENE_GENEINFO} || die "no config ENTREZ_GENE_GENEINFO";
    }
  }


  return new GeneSymbolMapper(
			      "-hgnc_file" => ($FLAGS{hgnc} || die "-hgnc"),
			      "-eg_file" => ($FLAGS{"gene-info"} || die "-gene-info")
			     );
}

sub get_gsm_ger {
  my ($no_nm) = @_;
  my $gsm = new_gsm();

  #
  # gene symbols reported by GENE_EXON_REGION:
  #
  my $ger_dir = $FLAGS{"gene-exon-region-dir"} || die "-gene-exon-region-dir";
  my @files = glob("$ger_dir/*txt");
  die unless @files;
  foreach my $f (@files) {
    open(IN, $f) || die;
    while (<IN>) {
      my @f = split /\t/, $_;
      my ($gene, $nm) = @f[0,2];
      $nm = "" if $no_nm;
      $gsm->add_gene(
		     "-gene" => $gene,
		     "-refseq" => $nm
		    );
    }
  }

  return $gsm;
}


sub cnv_gene_patch {
  my $gsm = get_gsm_ger();
  my $cnv_check = $FLAGS{cnv} || die "-cnv";
  my $cnv_annot = $FLAGS{"cnv-annotations"} || die;

  #
  #  patch CNVCheck.txt:
  #
  my $outfile = sprintf "%s.patched_to_GENE_EXON_REGION.txt", basename($cnv_check);
  my $wf = new WorkingFile($outfile);
  my $fh = $wf->output_filehandle();
  open(IN, $cnv_check) || die;
  while (<IN>) {
    chomp;
    my @f = split /\t/, $_;
    die unless @f == 3;
    my ($gene, $type, $desc) = @f;
    my $updated;
    unless ($gsm->contains($gene)) {
      my $result = $gsm->resolve("-symbol" => $gene);
      printf STDERR "miss for %s, map=%s\n", $gene, ($result || "n/a");
      $updated = $result if $result;
    }

    printf $fh "%s\n", join "\t", ($updated || $gene), $type, $desc;
  }
  $wf->finish();

  #
  #  patch Deletion.txt:
  #
  $outfile = sprintf "%s.patched_to_GENE_EXON_REGION.txt", basename($cnv_annot);
  my $df = new DelimitedFile(
			     "-file" => $cnv_annot,
			     "-headers" => 1,
			    );
  my $rpt = $df->get_reporter("-file" => $outfile);

  my $gf = "Gene Symbol";
  while (my $row = $df->get_hash()) {
    my $gene = $row->{$gf} || die;
    my $updated;
    unless ($gsm->contains($gene)) {
      my $result = $gsm->resolve("-symbol" => $gene);
      printf STDERR "miss for %s, map=%s\n", $gene, ($result || "n/a");
      $updated = $result if $result;
    }
    $row->{$gf} = $updated || $gene;
    $rpt->end_row($row);
  }
  $rpt->finish();
}

sub map_snv_gold_genes_to_sv {
  #
  # the SNV/indel gold gene list is also used in SV processing.
  # check compatibility of these symbols with FusionBuilder gene annotations.
  #

  my $rf = $FLAGS{"refflat-fb"} || die "-refflat-fb";
  printf "**** is %s the same refFlat file used by CLINICAL SV? [y/n]: ", $rf;
  if (-t STDOUT) {
    my $resp = <STDIN>;
    chomp $resp;
    die "quitting\n" unless (lc($resp) eq "y");
  } else {
    print STDERR "...it better be!\n";
  }

  my $gsm = get_gsm_fb();
  # FusionBuilder genes

  my $gold_genes_file = $FLAGS{gold} || die "-gold";
  my $outfile = sprintf "%s.patched_for_FB.txt", basename($gold_genes_file);

  my $wf = new WorkingFile($outfile);
  my $fh = $wf->output_filehandle();

  my $gold_genes = read_simple_file($gold_genes_file);

  my $manual = $FLAGS{"genes-manual"} || die "-genes-manual";
  my $df = new DelimitedFile("-file" => $manual,
			     "-headers" => 1,
			    );
  while (my $row = $df->get_hash()) {
    my $type = $row->{Type} || die;
    if ($type eq "Gold") {
      # gold gene
      my $gene = $row->{Gene} || die;
      push @{$gold_genes}, $gene;
    } else {
      die "unhandled gene Type $type";
    }
  }

  my %repair;
  foreach my $gene (@{$gold_genes}) {
    my $repaired;
    if ($gsm->contains($gene)) {
      printf STDERR "hit for %s (ok)\n", $gene;
    } else {
      $repaired = $gsm->resolve("-symbol" => $gene);
      printf STDERR "miss for %s (mapped to %s)\n", $gene, $repaired || "n/a";
    }
    printf $fh "%s\n", $repaired || $gene;
  }
  $wf->finish();

  printf STDERR "diff %s %s\n", $gold_genes_file, $outfile;

}

  sub hack_gsm {
    my $gsm = get_gsm_ger();
    #  die $gsm->resolve("-symbol" => "CMKOR1");

    my ($resolved, $hugo) = translate_symbol("-symbol" => "STL",
					     "-gsm" => $gsm);

    printf "resolved:%s hugo:%s\n", 
      ($resolved || "null"),
	($hugo || "null");
    die;
  }

sub hack_gsm_bare {
  my $gsm = get_gsm_ger(1);
  die $gsm->resolve("-symbol" => "CMKOR1");
  die;
}

sub erin_genes_sv {
  #
  # report genes used in SV processing (for Erin)
  #
  my $sv_main = get_sv_genes();

  my $gold_genes_patched = read_simple_file($FLAGS{"gold-genes-mapped-to-fb"} || die "-gold-genes-mapped-to-fb");

  my %all = %{$sv_main};
  foreach my $g (@{$gold_genes_patched}) {
    $all{$g} = 1;
  }

  foreach (sort keys %all) {
    printf "%s\n", $_;
  }
}

sub generate_germline_cnv_config {
  #
  # new version using TSOncoDB
  #
  my (%options) = @_;
  my $ram_mode = $options{"-ram"};

  my $map_to_ger = 0;
  # map symbols to GENE_EXON_REGION space.  Historically this has been
  # required, however we want to shift to HGNC dynamic lookups instead.

  printf STDERR "mapping to GENE_EXON_REGION: %s\n", $map_to_ger ? "y" : "n";

  my $gold_genes = read_simple_file($FLAGS{"gl-reportable-genes"} || die "-gl-reportable-genes");

  my @opts = (
	      "-gsm" => new_gsm_lite()
	     );
  if (0) {
    # original behavior based on SJ germline config file
    printf STDERR "restricting to SJ germline config only\n";
    push @opts, "-restrict_source" => "SJ_germline";
  } else {
    # more comprehensive annotations from various sources.
    # pro: might be more future-proof as new genes are added
    # con: old annotations might introduce inconsistencies,
    # e.g. APC has a recurrent/oncogen annotation in JZ_cancer annot
    push @opts, ("-exclude_source" => "SJ_JZ_cancer");
    # has some problematic annotations, some of these recurrent
    # events are truncations, e.g. in APC tumor suppressor,
    # so are not reliable oncogene annotations for CNV
    printf STDERR "using TSOncoDB for LoF/GoF annotations\n";
  }

  my $ts_onco_db = new TSOncoDB(@opts);

  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  my $outfile = sprintf "CNVCheck_germline_%d_%02d_%02d.txt", 1900+$year, $mon + 1, $mday;

  my ($wf, $fh, $gsm_ger);
  unless ($ram_mode) {
    $wf = new WorkingFile($outfile);
    $fh = $wf->output_filehandle();
  }

  $gsm_ger = get_gsm_ger() if $map_to_ger;

  my $errors = 0;

  my @out;
  foreach my $gene (@{$gold_genes}) {
    if ($map_to_ger) {
      if ($gsm_ger->contains($gene)) {
	# OK
      } elsif (my $alt = $gsm_ger->resolve("-symbol" => $gene)) {
	printf STDERR "mapping %s to %s\n", $gene, $alt;
	$gene = $alt;
	# annotation file should have been built with the equivalent
	# symbol as well
      } else {
	printf STDERR "ERROR: symbol %s is not in GENE_EXON_REGION\n", $gene;
	$errors++;
      }
    }

    if ($ts_onco_db->find_row("-gene" => $gene)) {
      my $is_ts = $ts_onco_db->is_lof("-gene" => $gene);
      my $is_onco = $ts_onco_db->is_gof("-gene" => $gene);
      my $annot = "";
      push @out, sprintf "%s", join "\t", $gene, "Del", $annot if $is_ts;
      push @out, sprintf "%s", join "\t", $gene, "Amp", $annot if $is_onco;
    } else {
      printf STDERR "ERROR: no TS/oncogene status for %s\n", $gene;
      $errors++;
    }
  }
  if ($errors) {
    die "errors, can't build CNV config file";
  } elsif (!$ram_mode) {
    foreach (@out) {
      printf $fh "%s\n", $_;
    }
    $wf->finish();
  }
  return \@out;
}

sub map_gl_extended_to_ger {
  #
  #  map germline reviewable gene symbols to SNV/indel/CNV namespace
  #  install under CLINCLS_GOLD_CANCER_RANGES_GER variable
  #
  # TO DO: add OPTIONS to:
  #    1. tartan import result file to dev
  #    2. tartan ?put? from dev to prod
  my $gsm = get_gsm_ger();

  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  my $outfile = sprintf "gl_gold_cancer_ranges_mapped_to_GER_%d_%02d_%02d.txt", 1900+$year, $mon + 1, $mday;

  my $wf = new WorkingFile($outfile);
  my $fh = $wf->output_filehandle();

  my $infile = $FLAGS{"gl-gold-cancer-ranges"} || die "-gl-gold-cancer-ranges";
  printf "infile: %s\n", $infile;
  open(IN, $infile) || die;

  my $count_same = 0;
  my $count_updated = 0;
  my $count_missing = 0;
  while (<IN>) {
    chomp;
    my @f = split /\t/, $_;
    die unless @f == 2;
    my $gene = $f[0];
    my $updated = $gene;
    if ($gsm->contains($gene)) {
      $count_same++;
    } else {
      my $result = $gsm->resolve("-symbol" => $gene);
      if ($result) {
	printf STDERR "mapping %s to %s\n", $gene, $result;
	$updated = $result;
	$count_updated++;
      } else {
	printf STDERR "missing %s\n", $gene;
	$count_missing++;
      }
    }
    printf $fh "%s\n", join "\t", $updated, $f[1];
    # same format with updated symbol
  }
  $wf->finish();

  printf "unchanged:%d  updated:%d  missing:%d\n",
    $count_same, $count_updated, $count_missing;
}

sub sanity_check_cosmic_somatic_genes_in_main_lists {
  my $gold_genes = flatfile_to_hash($FLAGS{"gold"} || die "-gold");
  my $silver_genes = flatfile_to_hash($FLAGS{"silver"} || die "-silver");
  # somatic SNV/indel classifier: main silver/gold gene lists

  my $cosmic_gold = $FLAGS{"cosmic-gold"} || die "-cosmic-gold";
  my $cosmic_silver = $FLAGS{"cosmic-silver"} || die "-cosmic-silver";

  my %cosmic;
  foreach my $fn ($cosmic_gold, $cosmic_silver) {
    open(IN, $fn) || die;
    while (<IN>) {
      chomp;
      my @f = split /\t/, $_;
      $cosmic{$f[0]} = 1;
    }
  }

  foreach my $gene (sort keys %cosmic) {
    next if $gold_genes->{$gene} or $silver_genes->{$gene};
    die "ERROR $gene not found in gold/silver gene list!";
    # doesn't happen
  }

  print STDERR "OK\n";

}

sub run_cnvs_integrated_old {
  # DELETE ME
  #
  # CNV classification for either somatic or germline.
  #
  my (%options) = @_;
  my $SOMATIC = 1;
  my $GERMLINE = 2;
  my $mode;

  my @cnv_configs;
  my $min_abs_copy_number;

  if ($options{"-somatic"}) {
    printf STDERR "starting SOMATIC CNV classification...\n";
    $mode = $SOMATIC;
    push @cnv_configs, $FLAGS{"cnv"} || die "specify -cnv";
    push @cnv_configs, $FLAGS{"gl-cnv"} || die "-gl-cnv";
    # 3/2016: somatic CNV should also report on germline reportable genes.
    # since germline CNV config only contains reportable genes, we
    # can just use that config as well.
  } elsif ($options{"-germline"}) {
    printf STDERR "starting GERMLINE CNV classification...\n";
    $mode = $GERMLINE;
    push @cnv_configs, $FLAGS{"gl-cnv"} || die "-gl-cnv";
  } else {
    die "must specify either -somatic or -germline";
  }
  die unless @cnv_configs;

  my (%amp_genes, %del_genes, %annotations);
  my %all_genes;
  my %sv_genes;

  #
  #  load target amplifications/deletions:
  #
  foreach my $cnv_check (@cnv_configs) {
    my $lines = read_simple_file($cnv_check);
    foreach my $line (@{$lines}) {
      my ($gene, $type, $annot) = split /\t/, $line;
      if ($type eq "Amp") {
	$amp_genes{$gene} = 1;
	$all_genes{$gene} = 1;
      } elsif ($type eq "Del") {
	$del_genes{$gene} = 1;
	$all_genes{$gene} = 1;
      } else {
	die;
      }

      if ($annot) {
	die "duplicate annotation for $gene, split by type??" if $annotations{$gene} and $annotations{$gene} ne $annot;
	$annotations{$gene} = $annot;
      }
    }
  }

  my %supplemental;
  if ($mode == $GERMLINE) {
    my $gl_reviewable = read_simple_file(($FLAGS{"gl-reviewable-genes-ger"} || die "-gl-reviewable-genes-ger"), "-hash1" => 1);
    # 3/2016: use mapping directly from Dale's master list which contains
    # both the HUGO and an alternate symbol.  This cuts out the middleman
    # vs. the germline intervals list originally created by Gang,
    # just in case there are issues with older reviewable symbols.
    foreach my $gene (keys %{$gl_reviewable}) {
      $all_genes{$gene} = 1;
      # if found in a CNV, these genes will get a silver
      # (we don't have tumor suppressor/oncogene annotations so
      # amplification/deletion related checks won't be performed)
    }
  } elsif ($mode == $SOMATIC) {
    # supplemental annotations:
    my $cnv_annot = $FLAGS{"cnv-annotations"} || die "specify -cnv-annotations";
    my $df = new DelimitedFile("-file" => $cnv_annot,
			       "-headers" => 1,
			      );
    # while (my $row = $df->next("-ref" => 1)) {  # headerless
    while (my $row = $df->get_hash()) {
      my $gene = $row->{"Gene Symbol"} || die;
      $supplemental{$gene} = $row->{"Tumour Types(Germline)"};
    }
  } else {
    die;
  }

  my $annot = $FLAGS{"gene-exon-region-dir"} || die "specify -gene-exon-region-dir";
  # need genomically-sorted genes in interval for analysis
  printf STDERR "loading %s...", $annot;
  my $ga = new GeneAnnotation(
			      "-style" => "gene_exon_region",
			      "-gene_exon_region_dir" => $annot,
			      "-ignore_non_coding" => 0
			     );
  print STDERR "done\n";

  my $sv_genes;
  # my $sv_genes = get_sv_genes();
  # FAIL:
  # - SV events have now been mapped to the gene symbol space used by
  #   FusionBuilder ("sharp" refFlat).
  # - this is NOT the same symbol set used by GENE_EXON_REGION.
  # - need a custom list of SV genes mapped to GENE_EXON_REGION space
  if ($mode == $SOMATIC) {
    $sv_genes = flatfile_to_hash($FLAGS{"sv-genes-mapped-to-ger"} || die "-sv-genes-mapped-to-ger");
  } elsif ($mode == $GERMLINE) {
    $sv_genes = {};
    # not applicable to germline mode since we don't do SV analysis there
  }

  foreach my $infile (@INPUT_CNV_FILES) {
    printf STDERR "processing %s...\n", $infile;

    my $df = new DelimitedFile(
			       "-file" => $infile,
			       "-headers" => 1,
			      );
    my @out_rows;
    while (my $row = $df->get_hash()) {
      #      my $genes = $row->{Genes};
      #      foreach (sort keys %{$row}) {
      #	printf "%s: %s\n", $_, $row->{$_};
      #      }
      my $log;
      if ($mode == $SOMATIC) {
	$log = $row->{LogRatio};
	die "ERROR: LogRatio field missing in input file" unless defined $log;
	# LogRatio only appears in paired output (somatic)
      } elsif ($mode == $GERMLINE) {
	$log = $row->{"seg.mean"};
	# unpaired (germline) runs, see
	# http://hc-wiki.stjude.org/display/compbio/Clinical+Genomics+CNV+Manual+Analysis+SOP
	# note this is NOT on log scale so shouldn't be used the same
	# way as in somatic classification.  In germline classification
	# we just use it to detect amplification/deletion status,
	# for high-level amplification the separate "seg.mean" check
	# still applies.
	die "ERROR: seg.mean field missing in input file" unless defined $log;
      } else {
	die;
      }

      my $medal;
      $row->{sort_score} = 0;

      my $chrom = $row->{chrom} || die;
      my $start = $row->{"loc.start"} || dump_die($row, "no loc.start");
      my $end = $row->{"loc.end"} || dump_die($row, "no loc.end");

      if ($start > $end) {
	# some records seem to have start and end swapped,
	# /cgs01/clingen/prod/tartan/runs/classify-cnv-somatic/SCKIrJDS/input/SJOS030272_D1_G1/WHOLE_GENOME/conserting-crest/SJOS030272_D1_G1_CONSERTING_Mapability_100.txt.QualityMerge
	# /home/medmonso/work/dale/2016_12_19_somatic_CNV_crash
	printf STDERR "WARNING: start > end for %s.%s.%s, swapping\n",
	  $chrom, $start, $end;
	my $temp = $start;
	$start = $end;
	$end = $temp;
      }

      #	die keys %{$row};
      #	die "fix me: reannotate @g $chrom $start $end";
      # /nfs_exports/genomes/1/projects/ClinicalSeq/medal_config/2013_09_20/medal_config_2013_09_20.txt

      $ga->find(
		"-reference" => $chrom,
		"-start" => $start,
		"-end" => $end
	       );
      #      my @genes_report = split /,/, $genes;
      my $genes_genomic = $ga->results_genes_genomic_order();
      # genes in interval, sorted genomically
      my $gene_count = scalar @{$genes_genomic};
      #      printf STDERR "genes: %d\n", $gene_count;
      my $is_intronic = $ga->is_intronic();

      #      printf STDERR "gene counts:%s report:%s us:%s\n",
      #	@genes_report == @{$genes_genomic} ? "same" : "different",
      #	  join(",", @genes_report), join(",", @{$genes_genomic});

      #      if ($genes eq "> 10 genes" or $genes eq "NA") {

      if (my $test_case = $FLAGS{"cnv-test"}) {
	printf STDERR "DEBUG, test case %d\n", $test_case;
	die if $FLAGS{"outfile-fq"};
	if ($test_case == 1) {
	  # focal deletion
	  $genes_genomic = ["PAX5"];
	  $log = -1;
	} elsif ($test_case == 2) {
	  # focal amplification, high-level amplification
	  $genes_genomic = ["EGFR"];
	  $log = CNV_MIN_LOG2_FOR_HLA;
	} elsif ($test_case == 3) {
	  # focal amplification, low-level amplification
	  $genes_genomic = ["EGFR"];
	  $log = CNV_MIN_LOG2_FOR_HLA - 1;
	} elsif ($test_case == 4) {
	  # reversed direction, single gene
	  $genes_genomic = ["EGFR"];
	  $log = -1;
	} elsif ($test_case == 5) {
	  # breakpoint gene, perfect match
	  #	  $genes_genomic = [1 .. 7, "EGFR"];
	  $genes_genomic = ["EGFR", 1 .. 7];
	  $log = CNV_MIN_LOG2_FOR_HLA;
	} elsif ($test_case == 6) {
	  # breakpoint gene, opposite match
	  #	  $genes_genomic = [1 .. 7, "EGFR"];
	  $genes_genomic = ["EGFR", 1 .. 7];
	  $log = - CNV_MIN_LOG2_FOR_HLA;
	} else {
	  die "unimplemented test case";
	}
	$gene_count = @{$genes_genomic};
      }

      # 10/16/2013:
      # add gene annotation columns, replacing now-retired
      # earlier annotation step
      $row->{"Number ofGenes"} = $gene_count;
      # (sic)
      $row->{Genes} = $gene_count > $CNV_MAX_GENES_ANNOTATE ? sprintf("> %d genes", $CNV_MAX_GENES_ANNOTATE) : join(",", @{$genes_genomic});

      my @medal_genes;
      my @reasons;
      my $is_focal = $gene_count <= CNV_MAX_GENES_FOCAL;

      #      $row->{AmpOrDel} = $log < 0 ? "Del" : "Amp";

      my $is_deletion;
      my $hash;
      if ($log < 0) {
	$is_deletion = 1;
	$hash = \%del_genes;
      } else {
	$hash = \%amp_genes;
      }
      die "ERROR, missing seg.mean" unless defined $row->{"seg.mean"};
      my $copy_number = $row->{"seg.mean"} * 2;

      if (CNV_HIGH_LEVEL_ENABLE) {
	# if a CNV is highly amplified or deleted
	if ($copy_number >= CNV_MIN_COPY_NUMBER_GAIN_FOR_GOLD or
	    $copy_number <= CNV_HIGH_LEVEL_DELETION_COPIES) {
	  # CNV is a high-level amplification or deletion
	  # search ALL contained genes vs. the watch list

	  printf "high-level event: %s, genes=%s\n", $copy_number, join ",", @{$genes_genomic};

	  foreach my $gene (@{$genes_genomic}) {
	    if ($hash->{$gene}) {
	      push @medal_genes, $gene;
	    }
	  }
	  if (@medal_genes) {
	    push @reasons, general_label($is_focal, $is_deletion);
	    # downstream users might want/expect this basic info,
	    # so leave it, however it seems a little misleading
	    # because these details are unrelated to the logic behind the
	    # event.
	    push @reasons, sprintf "high_level_%s", $is_deletion ? "deletion" : "amplification";
	    $medal = CLASS_GOLD;
	    $row->{Genes} = join(",", @{$genes_genomic});
	    # if medaled, annotate genes regardless of CNV size
	    # (remove usual cap)
	  }
	}
      }

      my $mode_germline_any = ($mode == $GERMLINE and $CNV_GERMLINE_MEDAL_ANY_SIZE);

      if ($medal) {
	# already handled
      } elsif ($mode_germline_any ? 0 : $gene_count > $CNV_MAX_GENES_MEDAL) {
	$medal = CLASS_UNKNOWN;
      } else {
	#
	#  focal event:
	#
	if ($is_focal) {
	  foreach my $g (@{$genes_genomic}) {
	    push @medal_genes, $g if $hash->{$g};
	  }
	  if (@medal_genes) {
	    push @reasons, general_label($is_focal, $is_deletion);
	    if ($is_deletion) {
	      $medal = CLASS_GOLD;
	    } else {
	      # amplification
	      if ($mode == $SOMATIC) {
		if ($log >= CNV_MIN_LOG2_FOR_HLA) {
		  $medal = CLASS_GOLD;
		  push @reasons, "high_level_amplification";
		} else {
		  $medal = CLASS_SILVER;
		  #		push @reasons, "low_level_amplification";
		  # remove for simplicity: only use "high_level_amplification"
		}
	      } elsif ($mode == $GERMLINE) {
		$medal = CLASS_GOLD;
		# JZ 5/23/2014:
		# I just had a brief discussion with Zhaojie. I think
		# that his concern is valid. We do need to implement CNV
		# for germline somewhat differently.  One big difference
		# is that for any CNVs that affect the coding exons of
		# the 31 genes currently used for reporting, we should
		# mark them as Gold. We can also be a bit more
		# sophisticated to make the tumor suppressor genes
		# “Gold” with deletion and oncogene “Gold” with
		# amplification.  We actually never implemented Germline
		# CNV or SV medal classification. It will be good if we
		# can throw something quick and simple together.
		#
		# JZ 5/27/2014:
		# For germline CNV, we do not need to distinguish
		# high-amplification or low amplification. Just treat
		# every germline amplification as if they were somatic
		# high-amplification.
		#
		# MNE notes:
		# - the genes in germline CNV config file currently
		#   perfectly match the the reportable gene list.
		#   So, this logic will assign gold if the amp/del status
		#   matches the config file (i.e. oncogene/tumor suppressor),
		#   and silver later if only the gene symbol matches.
	      } else {
		die;
	      }
	    }
	  }
	}

	if (not($medal) and
	    ($gene_count == 1 or $gene_count > CNV_MAX_GENES_FOCAL)) {
	  # - broad events
	  # - single-gene focal events that did not already receive a medal
	  my @search_genes;
	  if ($gene_count == 1) {
	    @search_genes = @{$genes_genomic};
	  } elsif ($mode_germline_any) {
	    @search_genes = @{$genes_genomic};
	  } else {
	    @search_genes = ($genes_genomic->[0],
			     @{$genes_genomic}[$#$genes_genomic]);
	  }

	  # first check: look for perfect matches,
	  # (will only ever happen for the non-focal search type)
	  foreach my $g (@search_genes) {
	    push @medal_genes, $g if $hash->{$g};
	  }
	  if (@medal_genes) {
	    die if $is_focal;
	    # "that's unpossible!"

	    if ($is_deletion) {
	      $medal = CLASS_GOLD;
	    } else {
	      # amplification
	      if ($mode == $SOMATIC) {
		if ($log >= CNV_MIN_LOG2_FOR_HLA) {
		  push @reasons, "high_level_amplification";
		  $medal = CLASS_GOLD;
		  # JZ 9/19/2013:
		  # For amplification, I am thinking about marking
		  # amplification >4x as Gold. This will represent those
		  # that have log2ratio >=2. I looked over the list of
		  # amplification and felt that most of them are well
		  # characterized oncogenes expected to have very
		  # high-level of amplification. We can mark anything
		  # with log2ration <=2 as silver.
		  #
		  # 10/9/2013: in light of the below should this be silver??
		} else {
		  # JZ 10/9/2013:
		  # I thought that we can assign bronze to the low-level
		  # amplification (>5 genes)?
		  $medal = CLASS_BRONZE;
		}
	      } elsif ($mode == $GERMLINE) {
		$medal = CLASS_GOLD;
		# see comments in focal code above
	      } else {
		die;
	      }
	    }
	  } else {
	    # second check: look for matches to all genes regardless
	    # of amplification/deletion status
	    foreach my $g (@search_genes) {
	      push @medal_genes, $g if $all_genes{$g};
	    }
	    if (@medal_genes) {
	      $medal = CLASS_SILVER;
	      push @reasons, "global_gene_match";
	      # match to global gene list regardless of amp/del status
	    }
	  }
	  if ($medal) {
	    if ($medal and @search_genes > 1) {
	      my $is_breakpoint;
	      foreach my $g (@medal_genes) {
		$is_breakpoint = 1 if $genes_genomic->[0] eq $g or
		  $genes_genomic->[$#$genes_genomic] eq $g;
	      }
	      unshift @reasons, "breakpoint_gene" if $is_breakpoint;
	    }
	    unshift @reasons, general_label($is_focal, $is_deletion);
	  }
	}			# broad

	if (not($medal)) {
	  #
	  #  silver for any gene in SV list:
	  #
	  my @search;
	  if ($is_focal or $gene_count == 1) {
	    @search = @{$genes_genomic};
	  } elsif ($mode_germline_any) {
	    @search = @{$genes_genomic};
	  } else {
	    @search = ($genes_genomic->[0],
		       @{$genes_genomic}[$#$genes_genomic]);
	  }

	  my $hit;
	  foreach (@search) {
	    if ($sv_genes->{$_}) {
	      push @medal_genes, $_;
	      $hit = 1;
	    }
	  }

	  if ($hit) {
	    $medal = CLASS_SILVER;
	    push @reasons, "sv_gene_match";
	    my $is_breakpoint;
	    foreach my $g (@medal_genes) {
	      $is_breakpoint = 1 if $genes_genomic->[0] eq $g or
		$genes_genomic->[$#$genes_genomic] eq $g;
	    }
	    unshift @reasons, "breakpoint_gene" if $is_breakpoint;
	    unshift @reasons, general_label($is_focal, $is_deletion);
	  }
	}
      }

      $row->{copy_number_gain} = $copy_number;
      if ($copy_number > 0 and $medal and $medal ne CLASS_UNKNOWN) {
	#
	#  supplemental/override logic for amplifications, per JZ 10/8/2013.
	#  only applies if a medal has been assigned (i.e. some kind of match)
	#
	if ($copy_number >= CNV_MIN_COPY_NUMBER_GAIN_FOR_GOLD) {
	  # a) For copy number gain >10, report amplification as gold regardless of num. of genes included in the interval.
	  # [presumably >= 10]
	  $medal = CLASS_GOLD;
	  #	  push @reasons, "high_level_amplification";
	  add_reason(\@reasons, "high_level_amplification");
	  # hopefully not a duplicate
	} elsif ($is_focal) {
	  $medal = CLASS_GOLD;
	  # b)For copy number gain <10, report gold for focal (<=5 genes)
	} else {
	  $medal = CLASS_SILVER unless $medal eq CLASS_GOLD;
	  # - silver for large segment (>5 genes)
	  # - don't override a better medal if assigned earlier
	}
      }

      push @reasons, sprintf "genes=%s", join "|", @medal_genes if @medal_genes;

      $medal = CLASS_UNKNOWN unless $medal;
      $row->{intronic_only} = $is_intronic;

      if ($mode == $GERMLINE and $is_intronic and $medal ne CLASS_UNKNOWN) {
	# would have received a medal, but intronic
	$medal = CLASS_UNKNOWN;
	unshift @reasons, "intronic";
      }

      if ($mode == $GERMLINE and
	  $medal ne CLASS_UNKNOWN and
	  defined $CNV_GERMLINE_MIN_COPY_DELTA) {
	unless (abs($copy_number) >= $CNV_GERMLINE_MIN_COPY_DELTA) {
	  push @reasons, "copy_number_below_threshold";
	  $medal = CLASS_UNKNOWN;
	}
      }

      $row->{GSBClass} = $medal;
      $row->{Reason} = join ",", @reasons;
      $row->{sort_score} = $MEDAL_SORT_WEIGHT{$medal} || die "no weight for $medal";

      my @annotations;
      foreach my $gene (@medal_genes) {
	if (my $thing = $annotations{$gene}) {
	  # won't be present if medal is based on a SV gene only
	  my $extra = $thing eq "COSMIC" ? $supplemental{$gene} : "";
	  if ($extra) {
	    printf STDERR "NOTE: patching supplemental info for %s: %s\n", $gene, $extra;
	    $thing .= ";" . $extra;
	  }
	  push @annotations, $thing;
	}
      }
      $row->{Evidence} = join "|", @annotations;

      push @out_rows, $row;
    }

    my @sorted = sort {$b->{sort_score} <=> $a->{sort_score}} @out_rows;
    my $rpt = $df->get_reporter(
				"-file" => get_outfile($infile),
				"-extra" => [
					     (
					      "Number ofGenes",
					      "Genes",
					      # 10/16/2013:
					      # we are now doing the
					      # gene annotation
					      "copy_number_gain",
					      "intronic_only",
					      "GSBClass",
					      "Reason",
					      "Evidence",
					     )
					    ]
			       );
    $rpt->write_headers();
    # force header line even if empty
    foreach my $row (@sorted) {
      $rpt->end_row($row);
    }
    $rpt->finish();
  }
  printf STDERR "CNVs done\n";
}

sub get_germline_reportable_genes {
  return read_simple_file($FLAGS{"gl-reportable-genes"} || die "-gl-reportable-genes");
}

sub add_reason {
  my ($reasons, $entry) = @_;
  push @{$reasons}, $entry unless grep {$_ eq $entry} @{$reasons};
}

sub meta_classifier_gene_list {
  #
  #  report union of all gene symbols used in somatic analysis
  #  (SNV/indel, SV, CNV).  Report each translated to
  #  SNV/indel, CNV (GENE_EXON_REGION) and SV (refFlat "sharp").
  #
  #  Use the raw/source reports rather than any translated versions.
  #
  #  see
  #  http://hc-wiki.stjude.org/display/compbio/Variant+Classification+Gene+Lists
  #

  my $report_type = $FLAGS{"meta-classifier-gene-list"} || die;
  my $use_somatic;
  my $use_germline;
  my $outfile;
  if ($report_type eq "somatic") {
    $use_somatic = 1;
    $use_germline = 0;
    $outfile = "somatic_gene_symbols.tab";
  } elsif ($report_type eq "germline") {
    $use_somatic = 0;
    $use_germline = 1;
    $outfile = "germline_gene_symbols.tab";
  } elsif ($report_type eq "all") {
    $use_somatic = 1;
    $use_germline = 1;
    $outfile = "germline_and_somatic_gene_symbols.tab";
  } elsif ($report_type eq "custom") {
    $use_somatic = 1;
    $use_germline = 1;
    $outfile = "custom_gene_symbols.tab";
  } else {
    die "must be somatic/germline/all/custom";
  }

  my %headered_file_flags;
  # delimited files, one or more headers
  my %column_flags;
  # delimited files, one or more indexes

  if ($use_somatic) {
    $column_flags{gold} = [ 0 ];
    $column_flags{silver} = [ 0 ];
    $column_flags{cnv} = [ 0 ];
    $column_flags{sv} = [ 0, 1 ];

    $column_flags{"gl-reportable-genes-ger"} = [ 0 ];
    $column_flags{"gl-reportable-genes-refflat"} = [ 0 ];
    # germline reportable genes now used in somatic

    $headered_file_flags{"genes-manual"} = [ "Gene" ];
    $headered_file_flags{"sv-silver-genes-fb"} = [ "gene" ];
    $headered_file_flags{"sv-manual"} = [ qw(pair1 pair2) ];
  }

  if ($use_germline) {
    $column_flags{"gl-cnv"} = [ 0 ];
    $column_flags{"gl-gold-cancer-ranges"} = [ 0 ];
    $column_flags{"gl-reportable-genes"} = [ 0 ];
  }

  #
  #  associate each source with somatic or germline:
  #
  my %flag2classifier;
  foreach my $f (qw(
		     gold
		     silver
		     cnv
		     sv
		     genes-manual
		     sv-silver-genes-fb
		     sv-manual
                     gl-reportable-genes-ger
                     gl-reportable-genes-refflat
		  )) {
    $flag2classifier{$f} = "somatic";
  }

  foreach my $f (qw(
		     gl-gold-cancer-ranges
		     gl-reportable-genes
		     gl-cnv
		  )) {
    $flag2classifier{$f} = "germline";
  }

  my %flag2type;
  #
  # associate each source with classifier type (SNV/indel, CNV, and/or SV).
  # See:
  #
  # http://hc-wiki.stjude.org/display/compbio/Variant+Classification+Gene+Lists
  #
  $flag2type{gold}{somatic_snv} = 1;
  $flag2type{gold}{somatic_sv} = 1;
  # this report deals with primary lists only since we are doing the
  # symbol remapping dynamically.  e.g. the "gold" gene list is used
  # by both SNV and SV, even though SV classifier uses the remapped
  # equivalent flag "-gold-genes-mapped-to-fb".
  $flag2type{silver}{somatic_snv} = 1;
  $flag2type{"genes-manual"}{somatic_snv} = 1;
  $flag2type{"genes-manual"}{somatic_sv} = 1;
  # SV: code doesn't actually use this (BUG)

  $flag2type{cnv}{somatic_cnv} = 1;

  foreach my $f (qw(sv sv-manual sv-silver-genes-fb)) {
    $flag2type{$f}{somatic_sv} = 1;
    $flag2type{$f}{somatic_cnv} = 1;
    # SV genes used by CNV classifier as well, e.g. EBF1
    # BUG: -sv-silver-genes-fb not seen by CNV classifier
  }

  foreach my $f (qw(gl-gold-cancer-ranges gl-reportable-genes)) {
    # -gl-cnv not included
    $flag2type{$f}{germline_snv} = 1;
    $flag2type{$f}{germline_cnv} = 1;
  }
  $flag2type{"gl-cnv"}{"germline_cnv"} = 1;

  $flag2type{"gl-reportable-genes-ger"}{"somatic_snv"} = 1;
  $flag2type{"gl-reportable-genes-ger"}{"somatic_cnv"} = 1;
  $flag2type{"gl-reportable-genes-refflat"}{"somatic_sv"} = 1;

  if ($report_type eq "custom") {
    foreach my $hash (\%column_flags, \%headered_file_flags) {
      foreach my $src (sort keys %{$hash}) {
	printf "include %s (%s: %s)? ",
	  $src,
	    $flag2classifier{$src},
	      join(",", sort keys %{$flag2type{$src}});
	my $response = <STDIN>;
	chomp $response;
	$response = lc($response);
	if ($response eq "y") {
	  # keep
	} elsif ($response eq "n") {
	  delete $hash->{$src};
	} else {
	  die "specify y or n";
	}
      }
    }
  }


  #
  #  load gene symbols:
  #
  my %sym2source;
  my %sym2classifier;

  foreach my $flag (sort keys %column_flags) {
    my $switch = "-" . $flag;
    my $fn = $FLAGS{$flag} || die $switch;
    my $indexes = $column_flags{$flag};
    my $classifiers = $flag2type{$flag} || die "no classifier associated with flag $flag";

    open(TMPIN, $fn) || die;
    while (<TMPIN>) {
      chomp;
      my @f = split /\t/, $_;
      foreach my $index (@{$indexes}) {
	if (my $sym = $f[$index]) {
	  # might not be defined, e.g. SV
	  $sym2source{$sym}{$flag} = 1;
	  foreach my $c (keys %{$classifiers}) {
	    $sym2classifier{$sym}{$c} = 1;
	  }
	}
      }
    }
  }

  foreach my $flag (sort keys %headered_file_flags) {
    my $fields = $headered_file_flags{$flag} || die;
    my $classifiers = $flag2type{$flag} || die;
    my $switch = "-" . $flag;
    my $fn = $FLAGS{$flag} || die $switch;
    my $df = new DelimitedFile(
			       "-file" => $fn,
			       "-headers" => 1,
			      );
    while (my $row = $df->get_hash()) {
      foreach my $field (@{$fields}) {
	die "no entry for $field" unless exists $row->{$field};
	if (my $sym = $row->{$field}) {
	  # might not be defined (SV)
	  $sym2source{$sym}{$flag} = 1;
	  foreach my $c (keys %{$classifiers}) {
	    $sym2classifier{$sym}{$c} = 1;
	  }
	}
      }
    }
  }


  my @labels = qw(
		   symbol_raw
		   sources
		   classifiers
		   symbol_HUGO
		   symbol_CNV
		   symbol_SV
		   lookup_ok
		   lookup_identical
		   missing_any
		   missing_both
		);

  push @labels, qw(
		    somatic_snv
		    somatic_cnv
		    somatic_sv
		 ) unless $report_type eq "germline";

  push @labels, qw(
		    germline_snv
		    germline_cnv
		 ) unless $report_type eq "somatic";

  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => \@labels,
			 "-auto_qc" => 1,
			);

  my $gsm_cnv = get_gsm_ger();
  my $gsm_sv = get_gsm_fb();

  if (0) {
    #    my $restrict = "SSX9";
    #    my $restrict = "SIL";
    my $restrict = "C15orf65";
    printf STDERR "DEBUG: restrict to %s only!\n", $restrict;
    foreach (keys %sym2source) {
      delete $sym2source{$_} unless $_ eq $restrict;
    }
    die unless %sym2source;
  }

  foreach my $sym (sort keys %sym2source) {
    my %r;

    my ($symbol_cnv, $symbol_hugo) = translate_symbol(
						      "-symbol" => $sym,
						      "-gsm" => $gsm_cnv
						     );

    my $symbol_sv = translate_symbol(
				     "-symbol" => $sym,
				     "-gsm" => $gsm_sv
				    );

    my @sources = sort keys %{$sym2source{$sym}};
    my %broad_classifiers = map {($flag2classifier{$_} || die "no classifier for $_") => 1} @sources;
    die "WTF: no classifiers for $sym" unless %broad_classifiers;

    $r{symbol_raw} = $sym;
    $r{sources} = join ",", @sources;
    $r{classifiers} = join ",", sort keys %broad_classifiers;

    my $used = $sym2classifier{$sym} || die;
    my $somewhere;
    foreach my $f (qw(
		       somatic_snv
		       somatic_cnv
		       somatic_sv
		       germline_snv
		       germline_cnv
		    )) {
      my $value = $used->{$f} || 0;
      $r{$f} = $value;
      $somewhere = 1 if $value;
    }
    die unless $somewhere;

    $r{symbol_HUGO} = $symbol_hugo;
    $r{symbol_CNV} = $symbol_cnv;
    $r{symbol_SV} = $symbol_sv;
    $r{lookup_ok} = ($symbol_cnv and $symbol_sv) ? 1 : 0;
    $r{lookup_identical} = (($symbol_cnv || "") eq $sym and
			    ($symbol_sv || "") eq $sym) ? 1 : 0;
    $r{missing_any} = (not($symbol_cnv) or not($symbol_sv)) ? 1 : 0;
    $r{missing_both} = (not($symbol_cnv) and not($symbol_sv)) ? 1 : 0;

    $rpt->end_row(\%r);
  }

  $rpt->finish();
}

sub translate_symbol {
  my (%options) = @_;
  my $gsm = $options{"-gsm"} || die;
  my $sym = $options{"-symbol"} || die;

  my $mapped = "";
  if ($gsm->contains($sym)) {
    # target database already contains this symbol.
    $mapped = $sym;
    $gsm->resolve("-symbol" => $sym);
    # run through lookup code anyway to determine HUGO
  } else {
    $mapped = $gsm->resolve("-symbol" => $sym);
  }
  return wantarray ? ($mapped, $gsm->approved_symbol) : $mapped;
}

sub sheila_sv_to_manual {
  my ($fn) = @_;

  my $sv_manual = $FLAGS{"sv-manual"} || die "-sv-manual";

  my $outfile = "sv_manual_appended.txt";

  my $df = new DelimitedFile(
			     "-file" => $sv_manual,
			     "-headers" => 1,
			    );
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			     );

  while (my $row = $df->get_hash()) {
    # copy original SV manual file
    $rpt->end_row($row);
  }

  my $df_new = new DelimitedFile(
				 "-file" => $fn,
				 "-headers" => 1,
				);
  my $fusion_col = "Gene Fusion";
  my %skip;
  $skip{">50 partners"} = 1;
  # already handled elsewhere

  my @pairs;
  while (my $row = $df_new->get_hash()) {
    die unless exists $row->{$fusion_col};
    my $fusion = $row->{$fusion_col};
    if ($fusion and !$skip{$fusion}) {
      # some rows blank (sections)
      $fusion =~ s/ fusions?$//;
      die $fusion if $fusion =~ /\s/;

      if ($fusion =~ /\//) {
	# a list of genes on LHS, each of which pairs with gene on RHS,
	# e.g. EAP/MDS1/EVI-RUNX1
	my @p = split /\-/, $fusion;
	die unless @p == 2;
	my @l = split /\//, $p[0];
	die unless @l > 1;
	foreach my $l (@l) {
	  push @pairs, [$l, $p[1]];
	}
      } elsif ($fusion =~ /\-/) {
	# single pair
	my @p = split /\-/, $fusion;
	die unless @p == 2;
	push @pairs, \@p;
      } else {
	# single gene
	push @pairs, [ $fusion, "" ];
      }
    }
  }

  foreach my $p (@pairs) {
    # pair1   pair2   contact date

    my @g = @{$p};
    foreach (@g) {
      $_ = uc($_) if /^ig\w$/i;
      # IgH, igK, etc.
    }

    my %r;
    @r{qw(pair1 pair2)} = @g;
    $r{contact} = "Sheila Shurtleff";
    $r{date} = "1/8/2015";
    # hack
    $rpt->end_row(\%r);
  }


  $rpt->finish();
}

sub meta_qc_checks {
  #
  #  various gene list QC/sanity checks
  #
  #  TO DO:
  #  - check for incompatible gene symbols in committee paneldecision files

  #
  #  are all germline reportable gene symbols present in the
  #  reviewable gene list?:
  #
  my $gl_reportable = get_germline_reportable_genes();

  my $gl_reviewable = $FLAGS{"gl-gold-cancer-ranges"} || die "-gl-gold-cancer-ranges";
  #  die $gl_reviewable;

  open(IN, $gl_reviewable) || die;
  my %gl_reviewable = map {(split /\t/, $_)[0] => 1} <IN>;
  foreach my $g (sort keys %gl_reviewable) {
    printf STDERR "ERROR: whitespace in reviewable gene symbol %s\n", $g if $g =~ /\s/;
  }

  foreach my $gl (@{$gl_reportable}) {
    printf STDERR "ERROR: reportable gene %s not in reviewable list\n", $gl
      unless $gl_reviewable{$gl};
    printf STDERR "ERROR: whitespace in reportable gene symbol %s\n", $gl if $gl =~ /\s/;
  }


  #
  #  do all germline reportable genes have an entry in TS/onco list?:
  #
  my $ts_onco_db = new TSOncoDB(
				"-gsm" => new_gsm_lite()
			       );

  foreach my $gene (@{$gl_reportable}) {
    my $found = $ts_onco_db->find_row("-gene" => $gene);
    printf STDERR "ERROR: no TSOncoDB entry for %s\n", $gene unless $found;
  }

  # TO DO:
  # - are all CNVCheck.txt gene symbols in GENE_EXON_REGION?
  # - GermlineReportableGeneAnnotation.txt: are all these in GENE_EXON_REGION?

}

sub export_panel_decisions {
  #
  #  export of GERMLINE panel decisions.
  #
  my $dbi = get_dbi_clinical();

  die "-gl-current-committee-medals" unless @COMMITTEE_GL_MEDALS;
  my %old;
  foreach my $orig_list (@COMMITTEE_GL_MEDALS) {
    # old committee calls for QC purposes
    my $df = new DelimitedFile("-file" => $orig_list,
			       "-headers" => 1,
			      );
    # while (my $row = $df->next("-ref" => 1)) {  # headerless
    my @f = qw(Chr WU_HG19_Pos ReferenceAllele MutantAllele);
    while (my $row = $df->get_hash()) {
      my @v = @{$row}{@f};
      foreach (@v) {
	die "undef field" unless defined $_;
      }
      $v[0] =~ s/^chr//;
      my $key = join "_", @v;
      $old{$key} = 1;
    }
  }

  my $query = "
WITH grouped_vc AS
  (SELECT variant_id,
          max(version_nbr) AS max_version_nbr
   FROM sith.variant_classifications
   GROUP BY variant_id)
SELECT vc.committee_classification,
       vc.origin,
       v.chr,
       v.wu_hg19_pos,
       v.referenceallele,
      v.mutantallele,
       v.gene_symbol,
       v.cache->'mc'->'class' as class,
       v.mrna_acc,
       v.aachange
FROM sith.variant_classifications vc
INNER JOIN grouped_vc ON (vc.variant_id = grouped_vc.variant_id
                          AND vc.version_nbr = grouped_vc.max_version_nbr)
JOIN sith.variants v ON (vc.variant_id=v.id)
WHERE vc.origin='GERMLINE'

";
  # Aman: TEMPORARY query, needs intervention from him to update

  my $rows = selectall_hashref($dbi, $query);
  $dbi->disconnect();

  #
  #  write output in the format used by the SNV/indel classifier
  #  so we can run the file through for QC
  #
  # http://hc-wiki.stjude.org/display/compbio/Variant+Classification+for+Clinical+Genomics+Details
  #
  # GeneName, Chr, WU_HG19_Pos, ReferenceAllele, MutantAllele, Class, AAChange, and mRNA_acc

  my $rpt = new Reporter(
			 "-file" => "panel_decisions_csith_gedi.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   GeneName
					   Chr
					   WU_HG19_Pos
					   ReferenceAllele
					   MutantAllele
					   Class
					   AAChange
					   mRNA_acc
					   paneldecision

					   origin
					   found_in_old
					)
				      ],
			 "-auto_qc" => 1,
			);

  my %long2short;
  $long2short{"PATHOGENIC"} = "P";
  $long2short{"LIKELY_PATHOGENIC"} = "LP";
  $long2short{"UNCERTAIN_SIGNIFICANCE"} = "U";
  $long2short{"LIKELY_BENIGN"} = "LB";
  $long2short{"BENIGN"} = "B";

  foreach my $row (@{$rows}) {
    #    dump_die($row, "debug", 1);

    my $origin = $row->{origin} || die;
    if ($origin eq "SOMATIC") {
      # this list is germline only,
      # somatic classifier has never tracked committee decisions
      next;
    } elsif ($origin ne "GERMLINE") {
      die "unknown origin $origin";
    }

    my %r;
    $r{origin} = $origin;
    $r{GeneName} = $row->{gene_symbol} || dump_die($row, "gene_symbol");
    $r{Chr} = "chr" . ($row->{chr} || die);
    $r{WU_HG19_Pos} = $row->{wu_hg19_pos} || die;
    $r{ReferenceAllele} = $row->{referenceallele} || die;
    $r{MutantAllele} = $row->{mutantallele} || die;
    my $class = $row->{class} || die;
    $class =~ s/^\"//;
    $class =~ s/\"$//;
    $r{Class} = $class;

    $r{AAChange} = $row->{aachange} || die;
    $r{mRNA_acc} = $row->{mrna_acc} || die;

    die unless exists $row->{committee_classification};
    my $decision = $row->{committee_classification} || next;
    # skip rows without a classification

    $r{paneldecision} = $long2short{$decision} || dump_die($row, "no short code for $decision");

    my @f = @r{qw(Chr WU_HG19_Pos ReferenceAllele MutantAllele)};
    $f[0] =~ s/^chr//;
    my $key = join "_", @f;
    my $found_in_old = $old{$key} ? 1 : 0;
    $old{$key} = 2 if $found_in_old;
    # track which old variants were found in the new set
    $r{found_in_old} = $found_in_old;

    #    $decision = "broken" unless $decision;

    if ($decision) {
      # don't include unless a paneldecision is present
      $rpt->end_row(\%r);
    } else {
      print "SKIPPING entry with no decision\n";
    }
  }
  $rpt->finish();

  my $missing = 0;
  my $found = 0;
  foreach my $k (keys %old) {
    my $status = $old{$k};
    if ($status == 1) {
      $missing++;
    } elsif ($status == 2) {
      $found++;
    } else {
      die;
    }
  }
  printf STDERR "old committee calls new: found=%d missing=%d\n", $found, $missing;

}

sub export_germline_reportable_genes {
  my $dbi = get_dbi_clinical();
  my $query = 'select * from rep.germline_reportable_genes order by gene_symbol;';
  export_query_to_flatfile(
			   "-dbi" => $dbi,
			   "-sql" => $query,
			   "-outfile" => "germline_reportable_genes.lst",
			   "-header" => 0,
			  );
  $dbi->disconnect();

}

sub get_dbi_clinical {
  my $db_user = $FLAGS{"gedi-user"} || getlogin() || getpwuid($<) || die;
  my $db_pw = $FLAGS{"gedi-password"} || die "-gedi-password";
  # TO DO:
  # this will eventually be replaced by config entries
  # (read-only user/pw?)

  my $dbi = get_dbi_gedi(
			 "-type" => "clinical",
			 "-password" => $db_pw,
			 "-dev" => $FLAGS{"gedi-dev"},
			) || die;
  # BLEH: get rid of this module?  (external distribution)
  return $dbi;
}

sub export_for_clinical {
  my $param = $FLAGS{"export-for-clinical"};
  my $clr = new CommandLineRebuilder("-parameters" => \@COMMAND_OPTIONS,
				     "-flags" => \%FLAGS);
  my $values = $clr->get_values_for("-param" => $param) || die;

  foreach my $f (@{$values}) {
    printf STDERR "%s\n", $f;
    die unless -s $f;
    my $out = basename($f);
    copy($f, $out) || die;
  }
}

sub pcgp_indel_hack {
  my $vm = get_vm_pcgp_somatic();
  my %row;
  $row{Chr} = 5;
  $row{WU_HG19_Pos} = 35874572;
  $row{ReferenceAllele} = "------TTAC";
  $row{MutantAllele} = "TCGCCCTGCA";

  my $hits_literal = $vm->find_literal_variant("-sj" => \%row);

  printf "hits: %d\n", $hits_literal ? scalar @{$hits_literal} : 0;
}

sub hack_exac {
  my $exac_fn = $FLAGS{"exac"} || die "-exac";
  my $tf_exac = new TabixFile("-file" => $exac_fn);
  $tf_exac->wiggle_bases(INDEL_NUCLEOTIDE_NT_FUZZY_MATCH);
  # pull in variants in the area

  my $row = {};
  if (0) {
    # apparently invalid AN (all 0)
    $row->{Chr} = 11;
    $row->{WU_HG19_Pos} = 2906165;
    $row->{ReferenceAllele} = "A";
    $row->{MutantAllele} = "G";
  } elsif (0) {
    # multiple ALT
    $row->{Chr} = 1;
    $row->{WU_HG19_Pos} = 13528;
    $row->{ReferenceAllele} = "C";
    $row->{MutantAllele} = "T";
  } elsif (0) {
    # deletion
    $row->{Chr} = 1;
    $row->{WU_HG19_Pos} = 144873963;
    $row->{ReferenceAllele} = "T";
    $row->{MutantAllele} = "-";
  } elsif (0) {
    # insertion test: works ok
    $row->{Chr} = 1;
#    $row->{WU_HG19_Pos} = 13417;
    $row->{WU_HG19_Pos} = 13418;
#    $row->{WU_HG19_Pos} = 13419;
    $row->{ReferenceAllele} = "-";
    $row->{MutantAllele} = "GAGA";
  } elsif (0) {
    # insertion test #2: not found??
    # it is findable, but not a gold gene
    $row->{Chr} = "chr14";
    $row->{WU_HG19_Pos} = 24470691;
    $row->{ReferenceAllele} = "-";
    $row->{MutantAllele} = "A";
  } elsif (1) {
    $row->{Chr} = "chr3";
    $row->{WU_HG19_Pos} = 38627165;
    $row->{ReferenceAllele} = "-";
    $row->{MutantAllele} = "TGTGTGTGTGTGTGTGTGTGTGTGTG";
    # 26 nt
  }
  $row->{GeneName} = "dummy";
  $row->{AAChange} = "dummy";

  my @reasons;

  my $vm = get_vm_exac(
		       "-row" => $row,
		       "-exac" => $tf_exac
		      );
  # build a VariantMatcher including only variants in region

  my $hits;
  if ($hits = $vm->find_snv(
			    "-sj" => $row
			   )) {
    die "ambiguous SNV hits" unless @{$hits} == 1;
  } elsif ($hits = indel_search($vm, $row)) {
    if (@{$hits} > 1) {
      # problematic: what to do here?
      push @reasons, "exac_multiple_indel_matches";
      $hits = [];
    }
  }

  if ($hits and @{$hits} == 1) {
    # usable match
    populate_exac(
		  "-reasons" => \@reasons,
		  "-hit" => $hits->[0],
		 );
  }

  die join ";", @reasons;

}

sub populate_exac {
  #
  #  populate Reason field info tags for ExAC match
  #
  my %options = @_;
  my $reasons = $options{"-reasons"} || die;
  my $hit = $options{"-hit"} || die;
  my $info_split = $hit->{INFO} || dump_die($hit, "no ExAc INFO");
  my $info_main = $hit->{vcf_hash}{INFO} || die;

  my %fmt = (
#	     "AF" => "%f",
	     "AF" => "%s",
	     "AC" => "%d",
	     "AN" => "%d",
	    );


  my @passthrough_tags = qw(AF AC AN);
  my @population_tags = qw(
			    AFR
			    AMR
			    EAS
			    FIN
			    NFE
			    SAS
			    OTH
			 );

  foreach my $tag (@passthrough_tags) {
    my $v = $info_split->{$tag};
    $v = $info_main->{$tag} unless defined $v;
    my $tag_print = sprintf "ExAC_%s", $tag;
    my $fmt = $fmt{$tag} || die "no format for $tag";
    push @{$reasons}, join "=", $tag_print, sprintf($fmt, $v) if defined $v;
  }

  my @pops;
  foreach my $pop (@population_tags) {
    my $tag_allele = sprintf 'AC_%s', $pop;
    my $tag_chroms = sprintf 'AN_%s', $pop;

    my $count_allele = $info_split->{$tag_allele};
    if (defined $count_allele) {
      my @out = $count_allele;

#      push @pop_ac, sprintf '%s:%d', $pop, $count_allele;

      my $count_chrom = $info_main->{$tag_chroms};
      if (defined($count_chrom)) {
#	my $freq = $count_allele / $count_chrom;
	push @out, $count_chrom;
	my $freq = $count_chrom ? ($count_allele / $count_chrom) : 0;
	# might be 0, e.g. 11.2906165.A.G
	push @out, $freq > 0 ? sprintf("%4.3e", $freq) : $freq;
      }

      push @pops, join ":", $pop, @out;
    }
  }
  push @{$reasons}, sprintf "ExAC_pop=%s", join ",", @pops if @pops;

}


sub get_vm_exac {
  my %options = @_;
  my $row = $options{"-row"} || die "-row";
  my $tf_exac = $options{"-exac"} || die "-exac";

  my $chr = $row->{Chr} || die "Chr";
  my $pos = get_sj_pos($row);

  my $vm = new VariantMatcher();

  my $vcf = $tf_exac->query(
                       "-reference" => $chr,
                       "-pos" => $pos,
                       "-vcf-parser" => 1,
                      );
  if ($vcf) {
    my $vcfu = new VCFUtils("-vcf" => $vcf);

    while (my $row_main = $vcf->next_data_hash()) {
#      dump_die($row_main, "exac debug: main line", 1);
      my $rows = $vcfu->get_alt_rows("-hash" => $row_main);
      # split row into one row for each alternate allele
      # with appropriately parsed-out alt info

      foreach my $r (@{$rows}) {
#	dump_die($r, "exac import", 1);
	my $chrom = $r->{CHROM} || die;
	my $pos = $r->{POS} || die;
	my $ref_base = $r->{REF} || die;
	my $var_base = $r->{ALT} || die;

	if (length($ref_base) == length($var_base)) {
	  # substitution (SNV, MNV, etc.)
	  $vm->add_snv(
		       "-row" => $r,
		       "-reference" => $chrom,
		       "-base-number" => $pos,
		       "-reference-base" => $ref_base,
		       "-variant-base" => $var_base
		      );
	} elsif (length($ref_base) == 1 and
		 length($ref_base) < length($var_base)) {
	  # simple insertion
	  die unless substr($ref_base, 0, 1) eq substr($var_base, 0, 1);
	  $vm->add_insertion(
			     "-row" => $r,
			     "-reference" => $chrom,
			     "-start" => $pos,
			     "-end" => $pos + 1,
			     "-count" => length($var_base) - 1
			    );
	  # hack (both before and after bases touched)
	} elsif (length($var_base) == 1 and
		 length($ref_base) > length($var_base)) {
	  # simple deletion
	  dump_die($r, "WARNING: VCF deletion formatting problem?", 1) unless substr($ref_base, 0, 1) eq substr($var_base, 0, 1);
	  # unexpected formatting e.g.
	  # CHROM: 20
	  # POS: 62550824
	  # REF: GGGACAGCGCCACGGAAGAGGACGCACCCGGCTGTGTGCACATGTGCCCAGGGCCCGGGACAGCGCCACGGAAGAGGACGCACCCGGCTGTGTGCACATGTGCCCAGGGCCCGGGACAGCGCCACGGAAGAGGACGCACCCGGCTGTGTGCACATGTGCCCAGGGCCCGGGACAGCGCCACGGAAGAGGACGCACCCGGCTGTGTGGACATGTGCCCAGGGCCCGGGACAGCGCCACGGAAGAGGACGCACA
	  # ALT: A
	  # vcf_alt_entry: 2

	  my $del_start = $pos + 1;
	  my $del_length = length($ref_base) - length($var_base);
	  my $del_end = $del_start + $del_length - 1;
#	  dump_die($r, "exac deletion $chrom $del_start $del_end", 1);
#	  dump_die($r->{INFO}, "INFO", 1);
	  $vm->add_deletion(
			    "-row" => $r,
			    "-reference" => $chrom,
			    "-start" => $del_start,
			    "-end" => $del_end,
			   );
	} else {
	  printf STDERR "ERROR: unhandled complex indel: chrom=$chrom pos=$pos ref=$ref_base var=$var_base\n";
	}
      }
    }
  }
  return $vm;
}

sub check_gold_isoforms {
  my $rpt_type = $FLAGS{"check-gold-isoforms"} || die;

  my $genes;
  my $outfile;
  if ($rpt_type eq "reportable") {
    my $f_rpt = $FLAGS{"gl-reportable-genes"} || die;
    $genes = read_simple_file($f_rpt);
    $outfile = "reportable_isoform_check.tab";
  } elsif ($rpt_type eq "reviewable") {
    my $ranges = $FLAGS{"gl-gold-cancer-ranges"} || die;
    open(IN, $ranges) || die;
    my %genes;
    while (<IN>) {
      chomp;
      my ($gene) = (split /\t/, $_)[0];
      $genes{$gene} = 1;
    }
    $genes = [ sort keys %genes ];
    $outfile = "reviewable_isoform_check.tab";
  } else {
    die "specify reviewable/reportable";
  }

  my $hgmd = $FLAGS{"hgmd-dm-mp-all"} || die;

  my $df = new DelimitedFile("-file" => $hgmd,
			     "-headers" => 1,
			     );
  printf STDERR "hgmd: %s\n", $hgmd;
  my %used;
  while (my $row = $df->get_hash()) {
    my $gene = $row->{gene};
    my $nm = $row->{refcore};
    $used{$gene}{HGMD} = $nm if $nm;
  }

  my $sjpi = get_sjpi();

  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   gene
					   SJ_preferred
					   HGMD
					   consistent
					)
				      ],
			 "-auto_qc" => 1,
			);

  foreach my $gene (sort @{$genes}) {
    my $pref = $sjpi->get_preferred_isoform($gene) || "sj_NA";
    my $hgmd = $used{$gene}{HGMD} || "hgmd_NA";

    my %r;
    $r{gene} = $gene;
    $r{HGMD} = $hgmd;
    $r{SJ_preferred} = $pref;
    $r{consistent} = $pref eq $hgmd ? 1 : 0;
    $rpt->end_row(\%r);
  }
  $rpt->finish();
}


sub add_batch_nhlbi_frequency {
  my %options = @_;

  my $rows = $options{"-rows"} || die;
  my $start_time = time;
  log_msg("batch NHLBI start");
  my $tabix = get_tabix_nhlbi();
  my $f_row_key = "_user_row";
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $twobit = $config_genome->{TWOBIT} || die "no TWOBIT";

  my @query_variants;
  foreach my $r (@{$rows}) {
    my $v = get_variant_from_row($r);
    $v->{$f_row_key} = $r;
    push @query_variants, $v;
  }

  my %map;
  $map{$FIELD_NHLBI_FREQ} = "genotype_frequency";

  my $tba = new TabixBatchAnnotation(
				     "-tabix" => $tabix,
				     "-twobit" => $twobit,
				     "-split_count" => TABIX_BATCH_SIZE_NHLBI,
				     "-f_tabix_chr" => "Chr",
				     "-f_tabix_pos" => "Pos",
				     "-f_tabix_ref_allele" => "ReferenceAllele",
				     "-f_tabix_var_allele" => "MutantAllele",

				     "-user_row_key" => $f_row_key,
				     "-annotation_map" => \%map,
				    );

  $tba->query(
	      "-query" => \@query_variants,
	     );

  log_msg(sprintf "batch NHLBI annotation for %d rows took %d",
	  scalar(@{$rows}), time - $start_time);
}

sub add_batch_cosmic {
  my %options = @_;

  if (not(%main::COSMIC_PUBMED) and my $cps = $FLAGS{"cosmic-pubmed-summary"}) {
    ram_debug("COSMIC summary start");
    printf STDERR "loading COSMIC PubMed/recurrence summary (%s)...\n", $cps;
    my $df = new DelimitedFile("-file" => $cps,
			       "-headers" => 1,
			     );
    my $db = \%main::COSMIC_PUBMED;
    my $db_aa = \%main::COSMIC_PUBMED_AA;
    # you don't see this

    my %pm;
    my %saw_pm;
    my %pm_aa;
    my %saw_pm_aa;

    while (my $row = $df->get_hash()) {
      my $key = get_cosmic_snv4($row);
      # safe to use uncooked chrom since this file was derived
      # from cooked cosmic file
      my $key_aa = get_cosmic_key($row);

      if ($row->{TotalVerifiedSample} >= $COSMIC_MIN_VALIDATED_SAMPLES_FOR_RECURRENT) {
	$db->{$key}{recurrent} = 1;
	$db_aa->{$key_aa}{recurrent} = 1;
      }

#      printf STDERR "current %s %s recurrent %d\n", $key, $aa_key, $db->{$key}{recurrent} || 0;
      # check this on a per-row rather than total bases:
      # db contains duplicate entries by AA because alternate transcripts
      # may be reported for the same sample

      if (my $pms = $row->{PubmedInfo}) {
	foreach my $pm (split /,/, $pms) {
	  next if $pm =~ /^0;/;
	  # observed but without a PMID
	  push @{$pm{$key}}, $pm unless $saw_pm{$key}{$pm};
	  $saw_pm{$key}{$pm} = 1;

	  push @{$pm_aa{$key_aa}}, $pm unless $saw_pm_aa{$key_aa}{$pm};
	  $saw_pm_aa{$key_aa}{$pm} = 1;

#	  printf STDERR "save pm %s %s\n", $key, $pm;
	}
      }
    }

    foreach my $key (keys %pm) {
      $db->{$key}{PubmedInfo} = join ",", @{$pm{$key}};
#      printf STDERR "final pm %s %s\n", $key, $db->{$key}{PubmedInfo};
    }

    foreach my $key_aa (keys %pm_aa) {
      $db_aa->{$key_aa}{PubmedInfo} = join ",", @{$pm_aa{$key_aa}};
    }

    ram_debug("COSMIC summary start");
  }

  my $germline_mode = $options{"-germline"};
  my $rows = $options{"-rows"} || die;
  my $start_time = time;
  log_msg("batch COSMIC start");
  my $tabix = get_tabix_cosmic();
  my $f_row_key = "_user_row";
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $twobit = $config_genome->{TWOBIT} || die "no TWOBIT";

  my @query_variants;
  foreach my $r (@{$rows}) {
    my $v = get_variant_from_row($r);
    $v->{$f_row_key} = $r;
    push @query_variants, $v;
  }

  my %map;
  $map{$FIELD_COSMIC} = "WU_HG19_Pos";

  my $tba = new TabixBatchAnnotation(
				     "-tabix" => $tabix,
				     "-twobit" => $twobit,
				     "-split_count" => TABIX_BATCH_SIZE_NHLBI,
				     "-f_tabix_chr" => "Chr",
				     "-f_tabix_pos" => "WU_HG19_Pos",
				     "-f_tabix_ref_allele" => "ReferenceAllele",
				     "-f_tabix_var_allele" => "MutantAllele",

				     "-user_row_key" => $f_row_key,
				     "-store_hits" => $FIELD_COSMIC_TABIX,
				    );

  $tba->query(
	      "-query" => \@query_variants,
	     );

  if ($germline_mode) {
    foreach my $row (@{$rows}) {
      if (my $hits = $row->{$FIELD_COSMIC_TABIX}) {
	$row->{$FIELD_COSMIC_TRUNCATING_FLAG} = 1 if cosmic_truncating_check($hits);
	delete $row->{$FIELD_COSMIC_TABIX};
	# free raw matched rows
      }
    }
  } else {
    foreach my $row (@{$rows}) {
      my $is_recurrent;
      my @pmi;
      my %saw_pmi;
      my $v_cosmic = "";
      my @xrefs;
      if (my $hits = $row->{$FIELD_COSMIC_TABIX}) {
	$v_cosmic = "Present";
	# set status here because hit is not guaranteed to have
	# pubmed info or be recurrent

	my @reasons;
	add_equivalent_tag(
			   "-row" => $row,
			   "-label" => "COSMIC",
			   "-tags" => \@reasons
			  );
	push @{$row->{$INTERNAL_FIELD_QUEUE_REASONS}}, @reasons if @reasons;

	foreach my $hit (@{$hits}) {
	  my $snv4 = get_cosmic_snv4($hit);
	  if (%main::COSMIC_PUBMED and my $xref = $main::COSMIC_PUBMED{$snv4}) {
	    push @xrefs, $xref;
	  }
	}

	delete $row->{$FIELD_COSMIC_TABIX};
	# free raw matched rows
      } elsif (%main::COSMIC_PUBMED_AA) {
	my $gene = $row->{$FIELD_GENE};
	my $aa = $row->{$FIELD_AACHANGE};
	my $key = sprintf "%s_p.%s", $gene, $aa;
	# this is such hot garbage
	if (my $xref = $main::COSMIC_PUBMED_AA{$key}) {
	  $v_cosmic = "Present";
	  push @xrefs, $xref;
	}
      }

      if (@xrefs) {
	foreach my $xref (@xrefs) {
	  $is_recurrent = 1 if $xref->{recurrent};
	  my $pmi = $xref->{PubmedInfo};
	  push @pmi, $pmi unless $saw_pmi{$pmi};
	  $saw_pmi{$pmi} = 1;
	  # hack: needs work
	}
      }

      # TO DO:
      # if specific tabix matching fails, try coarse gene/AA matching

      $row->{$FIELD_COSMIC} = $v_cosmic;
      $row->{$FIELD_COSMIC_RECURRENT} = $is_recurrent;
      $row->{$FIELD_COSMIC_PUBMED} = join ",", @pmi if @pmi;
    }
  }

  log_msg(sprintf "batch COSMIC annotation for %d rows took %d",
	  scalar(@{$rows}), time - $start_time);
}

sub add_batch_cosmic_old {
  # obsolete: this version uses the raw COSMIC data in tabix format
  my %options = @_;
  my $cosmic = $options{"-cosmic"} || die;
  my $rows_all = $options{"-rows"} || die;
  my $start_time = time;
  log_msg("batch COSMIC start");

  foreach my $rows (split_list($rows_all, TABIX_BATCH_SIZE_COSMIC)) {
    my @query_variants;
    foreach my $row (@{$rows}) {
      my $v = get_variant_from_row($row);
      push @query_variants, $v if $v->is_substitution;
    }

    my $t_rows = $cosmic->query(
				"-variants" => \@query_variants,
				"-hash" => 1
			       );

    my $vm_cosmic_somatic;
    if ($t_rows) {
      my $tfw = new TemporaryFileWrangler();
      my $tempfile = $tfw->get_tempfile("-append" => ".cosmic_tmp");
      #    printf STDERR "COSMIC tempfile: %s\n", $tempfile;

      my $rpt = new Reporter(
			     "-file" => $tempfile,
			     "-delimiter" => "\t",
			     "-labels" => [ sort keys %{$t_rows->[0]} ],
			     "-auto_qc" => 1,
			    );
      foreach my $r (@{$t_rows}) {
	$rpt->end_row($r);
      }
      $rpt->finish();
      # hack but easier to use a temporary file in this case
      # since code expects to iterate through a (very large)
      # COSMIC source file

      $vm_cosmic_somatic = parse_cosmic(
					"-file" => $tempfile,
					"-all-genes" => 1,
					"-all-variants" => 1
				       );
    } else {
      $vm_cosmic_somatic = new VariantMatcher();
    }

    foreach my $row (@{$rows}) {
      add_cosmic($row, $vm_cosmic_somatic);
    }
  }

  log_msg(sprintf "batch COSMIC annotation for %d rows took %d",
	  scalar(@{$rows_all}), time - $start_time);

}

sub add_batch_nsfp {
  my %options = @_;
  my $nsfp = get_tabix_nsfp();
  my $rows_all = $options{"-rows"} || die "rows";
  my $sjpi = $options{"-sjpi"} || die "sjpi";
  my $f_pos = $FLAGS{"tabix-dbnsfp-pos"} || die "-tabix-dbnsfp-pos";

  my $start_time = time;
  log_msg("batch dbNSFP start");

  foreach my $rows (split_list($rows_all, TABIX_BATCH_SIZE_DBNSFP)) {
    log_msg("batch of " . scalar(@{$rows}));
    my @query_variants;
    foreach my $row (@{$rows}) {
      foreach my $f (@NSFP_FIELDS) {
	# init blank entries in case not found
	$row->{$f} = "";
      }
      my $v = get_variant_from_row($row);
      push @query_variants, $v if $v->is_substitution;
    }

    my $t_rows = $nsfp->query(
			      "-variants" => \@query_variants,
			      "-hash" => 1
			     );
    my $vm;
    if ($t_rows) {
      $vm = new VariantMatcher();
      foreach my $r (@{$t_rows}) {
	my $v = new Variant();
	$v->import_dbnsfp_row("-row" => $r, "-f-pos" => $f_pos);
	$v->{dbnsfp} = $r;
	$vm->add_snv(
		     "-row" => $v,
		     "-variant" => 1,
		    );
      }
    } else {
      $vm = new VariantMatcher();
    }

    my $o_dbnsfp = get_nsfp();
    foreach my $row (@{$rows}) {
      # find any raw dbNSFP matches for this site
      my $v = get_variant_from_row($row);
      if (my $hits = $vm->find_snv("-variant" => $v)) {
	my @hit_rows = map {$_->{dbnsfp}} @{$hits};

	add_nsfp(
		 "-row" => $row,
		 "-dbnsfp" => $o_dbnsfp,
		 "-sj" => $row,
		 "-aa" => ($row->{$FIELD_AACHANGE} || ""),
		 "-nm" => $row->{$FIELD_REFSEQ},
		 "-preferred" => $sjpi,
		 "-hits" => \@hit_rows,
		 # pre-cooked hit rows for this site
		);

      }
    }
  }

  log_msg(sprintf "batch dbNSFP annotation for %d rows took %d",
	  scalar(@{$rows_all}), time - $start_time);

}


sub get_dbsnp_list {
  my ($hits) = @_;
  my %dbsnp;
  foreach my $h (@{$hits}) {
    my $dbsnp = $h->{dbsnp} || dump_die($h, "no dbsnp hash");
    $dbsnp{$dbsnp->{name} || die} = 1;
  }
  # we may actually match multiple SNP records, e.g.
  # chr19.40389657.T.G hits rs3746010 (-) and rs148187888 (+)
  return join ",", sort keys %dbsnp;
}

sub get_variant_from_row {
  my ($row) = @_;
  my $v = new Variant();
  $v->import_bambino_row("-row" => $row, "-postprocessed" => 1);
  return $v;
}

sub load_genome_config {
  my (%options) = @_;
  my $include_basename = $options{"-checksum"};
  my $checksum_mode = defined $include_basename;
  my $dump_mode = $FLAGS{"config-dump"};
  my $genome = $FLAGS{genome};
  %FLAGS = () if $checksum_mode or $dump_mode;

  my %var_string = map {$_, 1} qw(
				   DBNSFP_TABIX_POS_FIELD
				);

  my %var_simple = map {$_, 1} qw(
				   CLINCLS_SV_CHECK_MANUAL_FILE
				   CLINCLS_NHGRI_BRCA1_FILE
				   CLINCLS_GEDI_MUTATION_FILE
				   CLINCLS_RB1_FLAT_FILE
				   CLINCLS_GL_ARUP_FILE
				   CLINCLS_HGMD_CLASSIFICATION_EVERYTHING_FILE
				   CLINCLS_COSMIC_MUTANT_FILE
				   CLINCLS_GL_BAD_SNV_LIST_FILE
				   CLINCLS_COSMIC_SNV_INDEL_PUBMED_FILE
				   CLINCLS_CANCER_RELATED_GENES_FILE
				   CLINCLS_NHGRI_BRCA2_FILE
				   CLINCLS_GERMLINE_REPORTABLE_GENES
				   CLINCLS_CANCER_GENE_CENSUS_DELETION_FILE
				   CLINCLS_COSMIC_HOTSPOTS_FILE
				   CLINCLS_ASU_TERT_FILE
				   CLINCLS_GL_UMD_FLAT_FILE
				   CLINCLS_CNV_CHECK_FILE
				   CLINCLS_IARC_TP53_SOMATIC_FILE
				   CLINCLS_IARC_TP53_GERMLINE_FILE
				   GENE_TRANSCRIPT_MATRIX
				   CLINCLS_GL_PCGP_POPULATION_FILE
				   CLINCLS_GOLD_GENE_LIST_FILE
				   GL_APC_FLAT_FILE
				   CLINCLS_UNIPROT_ID_MAPPING_GZIP_FILE
				   GL_MSH2_FLAT_FILE
				   CLINCLS_GOLD_SNV_INDEL_FILE
				   CLINCLS_SV_CHECK_FILE
CLINCLS_SV_SILVER_GENES_FB
CLINCLS_GENES_MANUAL_FILE
CLINCLS_SILVER_SNV_INDEL_FILE
CLINCLS_SILVER_GENE_LIST_FILE

CLINCLS_CANCER_GENE_SJ_MUTATION_FILE
CLINCLS_NON_CANCER_GENE_SJ_MUTATION_FILE
CLINCLS_COMMITTEE_MEDALS_FILE
CLINCLS_COMMITTEE_MEDALS_FILE_SET2
CLINCLS_COMMITTEE_MEDALS_FILE_CSITH
CLINCLS_PRIVATE_VARIANTS_FILE

CLINCLS_COSMIC_PUBMED_SUMMARY
CLINCLS_PROMOTER_REGIONS
CLINCLS_PROMOTER_SITES
CLINCLS_MEDAL_CEREMONY_DB
				);

  # files to MD5 contents of

  my %var_size = map {$_, 1} (
			      "TWOBIT",
			      );
  # large files, check size only

  my %var_glob_size = (
		       "GENE_EXON_REGION_DIR" => "chr*region.txt",
		       "FASTA_CHR" => "*.fa",
		      );
  # sizes of files hit by glob patterns

  my %var_ignore;

  if ($genome) {
    my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";

    my %single;
    # a single config variable corresponds to a single parameter
    $single{gold} = "CLINCLS_GOLD_GENE_LIST_FILE";
    $single{silver} = "CLINCLS_SILVER_GENE_LIST_FILE";

    $single{cosmic} = "CLINCLS_COSMIC_SNV_INDEL_PUBMED_FILE";
    $single{"cosmic-gold"} = "CLINCLS_GOLD_SNV_INDEL_FILE";
    $single{"cosmic-silver"} = "CLINCLS_SILVER_SNV_INDEL_FILE";
    $single{"cosmic-cleaned"} = "CLINCLS_COSMIC_MUTANT_FILE";
    $single{"cosmic-hotspots"} = "CLINCLS_COSMIC_HOTSPOTS_FILE";

    $single{"gedi-recurrent"} = "CLINCLS_GEDI_MUTATION_FILE";

    $single{"cnv"} = "CLINCLS_CNV_CHECK_FILE";
    $single{"cnv-annotations"} = "CLINCLS_CANCER_GENE_CENSUS_DELETION_FILE";

    $single{sv} = "CLINCLS_SV_CHECK_FILE";
    $single{"sv-manual"} = "CLINCLS_SV_CHECK_MANUAL_FILE";

    $single{"gene-exon-region-dir"} = "GENE_EXON_REGION_DIR";
    $single{"genes-manual"} = "CLINCLS_GENES_MANUAL_FILE";

    $single{"asu-tert"} = "CLINCLS_ASU_TERT_FILE";

    $single{"iarc-tp53-germline"} = "CLINCLS_IARC_TP53_GERMLINE_FILE";
    $single{"iarc-tp53-somatic"} = "CLINCLS_IARC_TP53_SOMATIC_FILE";

    $single{"gl-gold-cancer-ranges"} = "CLINCLS_CANCER_RELATED_GENES_FILE";
    $single{"gl-arup"} = "CLINCLS_GL_ARUP_FILE";

    $single{"uniprot-idmapping"} = "CLINCLS_UNIPROT_ID_MAPPING_GZIP_FILE";

    $single{"preferred-isoforms"} = "GENE_TRANSCRIPT_MATRIX";
    # CLINCLS_TRANSCRIPT_MATRIX_FILE is obsolete: this is the
    # same preferred isoforms file used everywhere else

    $single{"gl-bad-snv"} = "CLINCLS_GL_BAD_SNV_LIST_FILE";
    $single{"rb1-flatfile"} = "CLINCLS_RB1_FLAT_FILE";

    $single{"gl-pcgp-population"} = "CLINCLS_GL_PCGP_POPULATION_FILE";

    $single{"nhgri-brca1"} = "CLINCLS_NHGRI_BRCA1_FILE";
    $single{"nhgri-brca2"} = "CLINCLS_NHGRI_BRCA2_FILE";
    $single{"gl-umd-flatfile"} = "CLINCLS_GL_UMD_FLAT_FILE";
    $single{"gl-apc-flatfile"} = "GL_APC_FLAT_FILE";
    $single{"gl-msh2-flatfile"} = "GL_MSH2_FLAT_FILE";
    $single{"hgmd-dm-mp-all"} = "CLINCLS_HGMD_CLASSIFICATION_EVERYTHING_FILE";

#    $single{"clinvar-flatfile"} = "CLINCLS_CLINVAR_VARIANTS_FLAT_FILE";
    $single{"sv-silver-genes-fb"} = "CLINCLS_SV_SILVER_GENES_FB";
    $single{"gl-reportable-genes"} = "CLINCLS_GERMLINE_REPORTABLE_GENES";

    $single{"fasta-dir"} = "FASTA_CHR";
    $single{"private-variants"} = "CLINCLS_PRIVATE_VARIANTS_FILE";
    $single{"tabix-dbnsfp-pos"} = "DBNSFP_TABIX_POS_FIELD";

    $single{"known-promoter-intervals"} = "CLINCLS_PROMOTER_REGIONS";
    $single{"known-promoter-sites"} = "CLINCLS_PROMOTER_SITES";
    $single{"2bit"} = "TWOBIT";
    $single{"cosmic-pubmed-summary"} = "CLINCLS_COSMIC_PUBMED_SUMMARY";
    $single{"sqlite"} = "CLINCLS_MEDAL_CEREMONY_DB";

    foreach my $p (keys %single) {
      unless (defined $FLAGS{$p}) {
	# command-line overrides genome config
	my $config_key = $single{$p};

	my $cv = $config_genome->{$config_key};
	unless ($cv) {
	  printf STDERR "ERROR: no value for config key $config_key\n" unless $var_ignore{$config_key};
	  next;
	}

	$FLAGS{$p} = $cv;
	if ($checksum_mode or $dump_mode) {
	  my $thing = $config_genome->{$config_key} || die;
	  if ($var_simple{$config_key}) {
	    # single variable and file
	    if ($dump_mode) {
	      dump_config_entry("-key" => $config_key);
	    } else {
	      my $md5;
	      if (-f $thing) {
		$md5 = md5_file($thing);
	      } else {
		$md5 = "file_missing";
	      }
	      my @things = ($config_key, $md5);
	      push @things, $include_basename == 2 ? $thing : basename($thing) if $include_basename;
	      printf "%s\n", join "\t", @things;
	    }
	  } elsif ($var_size{$config_key}) {
	    # very large single file
	    if ($dump_mode) {
	      dump_config_entry("-key" => $config_key);
	    } else {
	      my $md5 = -s $thing;
	      # simple proxy for very large files
	      my @things = ($config_key, $md5);
	      push @things, $include_basename == 2 ? $thing : basename($thing) if $include_basename;
	      printf "%s\n", join "\t", @things;
	    }
	  } elsif (my $pattern = $var_glob_size{$config_key}) {
	    my $dir = $config_genome->{$config_key} || die;
	    if ($dump_mode) {
	      dump_config_entry("-key" => $config_key);
	    } else {
	      my @files = glob($dir . "/" . $pattern);
	      if (@files) {
		foreach my $f (@files) {
		  my @things = $config_key . "." . basename($f);
		  push @things, -s $f;
		  push @things, $include_basename == 2 ? $thing : basename($thing) if $include_basename;
		  printf "%s\n", join "\t", @things;
		}
	      } else {
		my @things = ($config_key, "broken_no_files", "broken");
		printf "%s\n", join "\t", @things;
	      }
	    }
	  } elsif ($var_string{$config_key}) {
	    if ($dump_mode) {
	      dump_config_entry("-key" => $config_key);
	    } else {
	      my @things = ($config_key, $thing);
	      push @things, $include_basename == 2 ? $thing : basename($thing) if $include_basename;
	      printf "%s\n", join "\t", @things;
	    }
	  } elsif ($var_ignore{$config_key}) {
	    unless ($dump_mode) {
	      my @things = ($config_key, "ignored");
	      push @things, $include_basename == 2 ? $thing : basename($thing) if $include_basename;
	      printf "%s\n", join "\t", @things;
	    }
	  } else {
	    die "unhandled config key $config_key $thing";
	  }
	}
      }
    }

    #
    #  tabix directories:
    #  a TARTAn directory containing a single tabix .gz file
    #
    my %config2flag;
    $config2flag{THOUSAND_GENOMES_AF_DIR} = "tabix-thousand-genomes";
    $config2flag{EXAC_COOKED_DIR} = "exac-vcf2tab";
    $config2flag{EXAC_COVERAGE_DIR} = "exac-coverage";
    $config2flag{DBSNP_TABIX_DIR} = "tabix-dbsnp";
    $config2flag{NHLBI_TABIX_DIR} = "tabix-nhlbi";
    $config2flag{CLINCLS_COSMIC_TABIX_DIR} = "tabix-cosmic" if $ENABLE_TABIX_COSMIC;
    $config2flag{DBNSFP_TABIX_DIR} = "tabix-dbnsfp";
    $config2flag{CLINCLS_GEDI_TABIX_DIR} = "tabix-gedi";
    $config2flag{CLINVAR_TABIX_DIR} = "tabix-clinvar";
    $config2flag{CURATED_PROMOTER_REGIONS} = "known-promoter-intervals" if $ENABLE_KNOWN_PROMOTER_REGIONS;

    foreach my $cv (sort keys %config2flag) {
      if (my $tf = config2tabix($config_genome, $cv)) {
	my $flag = $config2flag{$cv};
	if ($FLAGS{$flag}) {
	  printf STDERR "tabix config: skipping %s, already specified on command line\n", $flag;
	} else {
	  $FLAGS{$flag} = $tf;
	  printf STDERR "tabix config: %s => %s\n", $flag, $tf;
	}

	if ($checksum_mode) {
	  my $md5;
	  if (1) {
	    $md5 = -s $tf;
	  } else {
	    # SLOW
	    $md5 = md5_file($tf);
	  }
	  my @things = ($cv, $md5);
	  push @things, $include_basename == 2 ? $tf : basename($tf) if $include_basename;
	  printf "%s\n", join "\t", @things;
	} elsif ($dump_mode) {
	  dump_config_entry("-key" => $cv);
	}
      }
    }

    #
    # params where each may be specified multiple times to build up a list:
    #
    my %set;

    $set{"gl-current-committee-medals"} = [
					   \@COMMITTEE_GL_MEDALS,
					   "CLINCLS_COMMITTEE_MEDALS_FILE_CSITH",
					   "CLINCLS_COMMITTEE_MEDALS_FILE",
					   "CLINCLS_COMMITTEE_MEDALS_FILE_SET2",
					  ];

    $set{"gedi-recurrent-restrict"} = [
				       \@GEDI_RECURRENT_RESTRICT,
				       "CLINCLS_CANCER_GENE_SJ_MUTATION_FILE",
				       "CLINCLS_NON_CANCER_GENE_SJ_MUTATION_FILE",
				      ];





    foreach my $p (keys %set) {
      my ($listref, @keys) = @{$set{$p}};

      # die "list already populated for $p, how to handle?" if @{$listref} and !$checksum_mode;
      # debug, e.g. test new -gl-current-committee-medals

      foreach my $key (@keys) {
	my $value = $config_genome->{$key};
	unless ($value) {
	  printf STDERR "ERROR: no value for config key %s\n", $key unless $var_ignore{$key};
	  next;
	}

	push @{$listref}, $value;
	if ($checksum_mode) {
	  if ($var_simple{$key}) {
	    my @things;
	    if (-f $value) {
	      my $md5 = md5_file($value);
	      @things = ($key, $md5);
	      push @things, $include_basename == 2 ? $value : basename($value) if $include_basename;
	    } else {
	      @things = ($key, "broken_no_file", "broken");
	    }
	    printf "%s\n", join "\t", @things;
	  } else {
	    die "unhandled $key $value";
	  }
	} elsif ($dump_mode) {
	  dump_config_entry("-key" => $key);
	}
      }
    }

    if ($checksum_mode) {
      my %ignore_params = map {$_, 1} (
				       "-config",
				       "-fb-column",
				       "-rsc",
				       "-indel-match-size",
				       "-config-checksum",
				       "-genome",
				       "-in",
				       "-outfile-fq",
				       "-hgmd-table",
				       "-nhlbi-x",
				       "-clinvar-manual",
				       "-one-list-to-rule-them-all",
				       "-collapse-meta-list",
				       "-germline-gene-survey",
				       "-parse-rb1",
				       "-map-rb1",
				       "-meta-reannotate",
				       "-mini-gl-run",
				       "-seatbelts-are-for-sissies",
				       "-rb1-type",
				       "-rb1-fasta",
				       "-show-isoforms-used",
				       "-committee-conflict-warn",
				       "-sv-gene-check",
				       "-cnv-gene-patch",
				       "-cnv-sv-gene-patch",
				       "-verify-cnv-genes",
				       "-map-snv-gold-genes-to-sv",
				       "-refflat-fb",
				       "-erin-genes-cnv",
				       "-erin-genes-sv",
				       "-cnv-germline",
				       "-map-gl-extended-to-ger",
				       "-sanity-somatic-cosmic-genes-in-main-lists",
				       "-meta-classifier-gene-list",
				       "-sheila-sv-to-manual",
				       "-meta-qc",
				       "-gedi-user",
				       "-gedi-password",
				       "-gedi-dev",
				       "-committee-export",
				       "-config-checksum-old",
				       "-export-for-clinical",
				       "-check-gold-isoforms",
				       "-gl-custom-genes",
				       "-gl-custom-allow-interval-fail",
				       "-gl-no-reviewable",
				       "-enable-thousand-genomes",
				       "-tabix",
				       "-config-dump",

				       # now obsolete:
				       "-exac",

				       # these might be part of formal
				       # config eventually (?):
				       "-hgnc",
				       "-gene-info",

				       "-enable-taylor-hotspots",
				       "-taylor-hotspots",
				       # dev

				       "-enable-known-promoter-regions",

				       "-cnv-max-genes-annotate",
				       # TO DO: might this move to genome
				       # config at some point?
				       "-cnv-germline-any-size",
				       "-cnv-germline-min-copy-delta",
				       "-config-somatic-cnv-add-manual-genes",
				       "-max-population-frequency",
				       "-count-medalable-nt",

				       "-enable-sv-gsm",
				       "-enable-cnv-gsm",
				       "-gl-custom-genes-lite",
				       "-clinvar-min-gold-stars",
				      );

      $ignore_params{"-tabix-cosmic"} = 1 unless $ENABLE_TABIX_COSMIC;

      my @all_flags = grep {not(ref $_)} @COMMAND_OPTIONS;
      my $broken;

      my %c2f = map {$_, 1} values %config2flag;

      foreach my $flag (@all_flags) {
	$flag =~ s/=.*//;
	next if $ignore_params{$flag};

	next if $flag =~ /hack/;
	next if $flag =~ /\-test/;
	next if $flag =~ /debug/;
	# debug/development flags
	next if $flag =~ /^\-no\-/;
	# logical switches
	next if $flag =~ /^\-in\-/;
	next if $flag =~ /^\-single\-/;
	# input file spec
	next if $flag =~ /^\-dump\-/;
	# data export
	next if $flag =~ /^\-generate\-/;
	# utility
	next if $flag =~ /verbose/;
	next if $flag =~ /^\-enable\-tabix/;
	# tabix enable/disable

	$flag =~ s/^\-//;
	my $status;
	if ($single{$flag}) {
	  $status = "single";
	} elsif ($c2f{$flag}) {
	  $status = "config2flag";
	} elsif ($set{$flag}) {
	  $status = "set";
	}
	unless ($status) {
	  printf STDERR "ERROR: unhandled parameter %s\n", $flag;
	  $broken = 1;
	}
      }
      die "ERROR: unhandled parameter(s)!" if $broken;

      #
      #  also include some species config entries
      #

      my $species = "Homo_sapiens";
      my $config_species = TdtConfig::readConfig('species', $species) || die "can't find species config for $species";

      foreach my $param (qw(
			     HGNC
			     OMIM_MIM2GENE
			     TS_ONCO_DB
			  )) {
	my $value;
	my $file;
	if ($file = $config_species->{$param} and -s $file) {
	  $value = md5_file($file);
	} else {
	  $value = "missing";
	}

	printf "%s\n", join "\t", $param, $value;
      }
    }
  }

  exit(0) if $dump_mode;

}

sub check_required_snv_indel_headers {
  # columns required for SNV/indel classification
  my %options = @_;
  my $headers = $options{"-headers"} || die;
  my $file = $options{"-file"} || die;
  my %headers = map {$_, 1} @{$headers};
  my @required = (
		  $FIELD_GENE,
		  $FIELD_CHR,
		  $FIELD_REF_ALLELE,
		  $FIELD_MUT_ALLELE,
		  $FIELD_CLASS,
		  $FIELD_AACHANGE,
		  $FIELD_REFSEQ
		  );
  my %missing;
  foreach my $h (@required) {
    $missing{$h} = 1 unless $headers{$h};
  }
  if (%missing) {
    die sprintf "ERROR: file %s is missing required column(s) %s", $file, join ",", sort keys %missing;
  }

  my $found_pos;
  foreach my $h ($FIELD_VPOS, $FIELD_VPOS_ALT) {
    $found_pos = 1 if $headers{$h};
  }
  die sprintf "ERROR: file %s is missing required position column (%s or %s)", $file, $FIELD_VPOS, $FIELD_VPOS_ALT unless $found_pos;


}

sub dbnsfp_hack {
  die unless @INPUT_SNV_INDEL_FILES;

  my $dbnsfp = get_nsfp();
  my $sjpi = get_sjpi();

  foreach my $infile (@INPUT_SNV_INDEL_FILES) {
    my $df = new DelimitedFile(
			       "-file" => $infile,
			       "-headers" => 1,
			      );
    while (my $row = $df->get_hash()) {
      add_nsfp(
	       "-row" => $row,
	       "-dbnsfp" => $dbnsfp,
	       "-sj" => $row,
	       "-aa" => ($row->{AAChange} || die),
	       "-nm" => $row->{$FIELD_REFSEQ},
	       "-preferred" => $sjpi,
	      );

      die;
    }
  }
}

sub get_vm_private {
  #
  # private/internal variants
  #
  my $fn = $FLAGS{"private-variants"} || die "-private-variants";

  printf STDERR "parsing private variants (%s)...", $fn;
  my $df = new DelimitedFile(
			     "-file" => $fn,
			     "-headers" => 1,
			    );
  my $vm = new VariantMatcher();
  while (my $row = $df->get_hash()) {
    my $chr = $row->{$FIELD_CHR};
    my $pos = $row->{$FIELD_VPOS};
    my $ra = $row->{$FIELD_REF_ALLELE};
    my $va = $row->{$FIELD_MUT_ALLELE};
    my $gene = $row->{$FIELD_GENE};
    my $nm = $row->{$FIELD_REFSEQ};
    my $aa = $row->{$FIELD_AACHANGE};

    my $imported;

    if ($chr and $pos and $ra and $va) {
      if ($ra eq "-" or $va eq "-") {
	die "indel, not implemented";
      } else {
	# SNV
	$vm->add_snv("-row" => $row,
		     "-sj" => 1);
	$imported = 1;
      }
    }

    if ($gene and $nm and $aa) {
      $vm->add_aa(
		  "-gene" => $gene,
		  "-aa" => $aa,
		  "-row" => $row
		 );
      $imported = 1;
    }

    dump_die($row, "can't import") unless $imported;
  }
  print STDERR "done\n";
  return $vm;
}

sub committee_export {
  my $rows = get_vm_committee("-rows" => 1);
  my $rpt = new Reporter(
			 "-file" => "committee_reclassify.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   Chr
					   Pos
					   Chr_Allele
					   Alternative_Allele
					   paneldecision_orig
					)
				      ],
			 "-auto_qc" => 1,
			);
  # use raw Bambino format so Annovar will adjust indel coordinates

  foreach my $r (@{$rows}) {
    my %r = %{$r};
    $r{Pos} = $r->{WU_HG19_Pos} || die;
    $r{Chr_Allele} = $r->{ReferenceAllele} || die;
    $r{Alternative_Allele} = $r->{MutantAllele} || die;
    $r{paneldecision_orig} = $r->{paneldecision} || die;
    $rpt->end_row(\%r);
  }
  $rpt->finish;
}

sub add_custom_gl_genes {
  my (%options) = @_;
  my $grf_gold = $options{"-grf-gold"} || die;
  my $gl_gold_genes = $options{"-gl-gold-genes"} || die;
  my $gsm_gold_genes = $options{"-gsm-gold-genes"} || die;
  my $f_gl_custom = $FLAGS{"gl-custom-genes"};
  my $f_gl_custom_lite = $FLAGS{"gl-custom-genes-lite"};
  if ($f_gl_custom) {
    #
    # user has additional genes to process (e.g. extensions for St. Jude Life).
    #
    my $genome = $FLAGS{genome} || die "-genome";
    my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
#    my $f_refflat = $config_genome->{REFSEQ_REFFLAT} || die;
    # meh; prefer annovar instead
#    require Bambino2Annovar;
#    my $b2a = new Bambino2Annovar();
#    my $code = 'requre Bogus;';
    my $code = 'require Bambino2Annovar;';
    eval $code;
    if ($@) {
      die "can't load Bambino2Annovar.pl: custom genes code requires clinload snv-annovar";
    }

    my $b2a = new Bambino2Annovar();
    my $f_refflat = $b2a->get_refgene_db("-genome" => $genome);

    my $rff = new RefFlatFile("-canonical_references_only" => 1);
    $rff->parse_file("-refflat" => $f_refflat,
		     "-type" => "refgene");
    my $rows = $rff->get_rows();
    my %g2r;
    foreach my $r (@{$rows}) {
      my $g = $r->{gene};
      push @{$g2r{$g}}, $r if $g;
    }

    #
    # genes user wants to use:
    #
    my $df = new DelimitedFile("-file" => $f_gl_custom,
			       "-headers" => 1,
			      );

    my @fail;
    my %saw_genes;
    while (my $row = $df->get_hash()) {
      my $gene = $row->{gene} || die;
      die "duplicate gene $gene" if $saw_genes{$gene};
      $saw_genes{$gene} = 1;
      if ($gene) {
	$gsm_gold_genes->add_gene("-gene" => $gene);
	if (my $wanted = $g2r{$gene}) {
	  $gl_gold_genes->{$gene} = 1;
	  my $intervals = $rff->find_gene_interval(
						   "-gene" => $gene,
						   "-buffer" => $REVIEWABLE_FLANKING_BUFFER_NT,
						   "-single" => 0
						  );

	  foreach my $interval (@{$intervals}) {
	    # add a separate interval for each transcript:
	    # sometimes there are mapping to multiple locations/strands/chroms
	    my ($chr, $start, $end) = @{$interval};
	    my $loc = sprintf '%s:%d-%d', $chr, $start, $end;
	    printf STDERR "adding custom gene interval %s: %s\n", $gene, $loc;
	    $grf_gold->add(
			   "-range" => $loc,
			   "-value" => {
					"gene" => $gene,
					"location" => $loc
				       },
			  );
	  }
	} else {
	  push @fail, $gene;
	}
      }
    }
    if (@fail) {
      my $gsm = new_gsm();
      foreach my $g (keys %g2r) {
	$gsm->add_gene("-gene" => $g);
      }

      my %alt;
      foreach my $g (sort @fail) {
	$alt{$g} = $gsm->resolve("-symbol" => $g);
      }

      printf STDERR "*** ERROR: can't find region for the following genes:\n";
      foreach my $g (sort @fail) {
	printf STDERR "  %15s: ", $g;
	if (my $alt = $alt{$g}) {
	  printf STDERR "possibly %s?\n", $alt;
	} else {
	  print STDERR "no alternate found\n";
	}
      }
      die "fatal error: unresolved gene symbols.  Specify -gl-custom-allow-interval-fail to drop these genes and continue anyway." unless $FLAGS{"gl-custom-allow-interval-fail"};
    }
  }

  if ($f_gl_custom_lite) {
    # just gene symbols
    my $genes = read_simple_file($f_gl_custom_lite, "-hash1" => 1, "-tokenize" => 1);
    foreach my $gene (keys %{$genes}) {
      $gl_gold_genes->{$gene} = 1;
      $gsm_gold_genes->add_gene("-gene" => $gene);
    }
  }


}

sub add_custom_ts {
  my (%options) = @_;
  my $ts_genes = $options{"-ts-genes"} || die;
  my $f_gl_custom = $FLAGS{"gl-custom-genes"};
  my $ts_onco_db = $options{"-ts-onco-db"};
  if ($f_gl_custom) {
    #
    # user has additional genes to process (e.g. extensions for St. Jude Life).
    #
    confess "-ts-onco-db" if $GERMLINE_TS_ONCO_DB and !$ts_onco_db;
    my $df = new DelimitedFile("-file" => $f_gl_custom,
			       "-headers" => 1,
			      );

    while (my $row = $df->get_hash()) {
      my $gene = $row->{gene} || die;
      my $tgold = $row->{truncation_gold};
      dump_die($row, "no truncation_gold") unless defined $tgold;
      if ($gene) {
	$ts_genes->{$gene} = $tgold;
	if ($GERMLINE_TS_ONCO_DB) {
	  my $is_ts = 0;
	  my $is_onco = 0;
	  if ($tgold == 1) {
	    $is_ts = 1;
	    $is_onco = 0;
	  } elsif ($tgold == 0) {
	    $is_ts = 0;
	    $is_onco = 1;
	  } elsif ($tgold == GERMLINE_ANNOTATION_BOTH_TS_AND_ONCO) {
	    $is_ts = $is_onco = 1;
	  }
	  $ts_onco_db->add_gene_info(
				     "-gene" => $gene,
				     "-is-ts" => $is_ts,
				     "-is-onco" => $is_onco
				    );
	}
      }
    }
  }
}

sub hack_exac2 {
  my $infile = $FLAGS{in} || die "-in";

  my $exac_fn = $FLAGS{"exac"} || die "-exac";
  printf STDERR "ExAC: %s\n", $exac_fn;
  my $tf_exac = new TabixFile("-file" => $exac_fn);
  $tf_exac->wiggle_bases(INDEL_NUCLEOTIDE_NT_FUZZY_MATCH);

  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			     );
  my $rpt = $df->get_reporter("-file" => $infile . ".exac.tab",
			      "-extra" => [
					   "ExAC_AF",
					   "ambiguous",
					  ]
			     );

  while (my $row = $df->get_hash()) {
    my $vm_exac = get_vm_exac(
			      "-row" => $row,
			      "-exac" => $tf_exac
			     );
    my @reasons;
    my $hits = $vm_exac->find_snv(
				  "-sj" => $row
				 );
    $hits = indel_search($vm_exac, $row) unless $hits;

    my $af = "";
    my $ambiguous = 0;
    if ($hits) {
      if (@{$hits} == 1) {
	# usable single match
	#	    dump_die($row, "debug before exac", 1);
	populate_exac(
		      "-hit" => $hits->[0],
		      "-reasons" => \@reasons,
		     );
      } else {
	$ambiguous = 1;
	# multiple matches
      }

      foreach (@reasons) {
	$af = sprintf "%.9f", $1 if /ExAC_AF=(.*)/;
      }
    }
    $row->{ExAC_AF} = $af;
    $row->{ambiguous} = $ambiguous;
    $rpt->end_row($row);
  }
  $rpt->finish();
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

sub hack_refgene_sv_symbols {
  my $sv_genes = get_sv_genes();
  my $gsm = get_gsm_fb();

  foreach my $gene (sort keys %{$sv_genes}) {
    my $status;
    if ($gsm->contains($gene)) {
      $status = "identical";
    } elsif (my $new = $gsm->resolve("-symbol" => $gene)) {
      $status = "resolved to $new";
    } else {
      $status = "missing";
    }
    printf "status for %s: %s\n", $gene, $status;
  }
}

sub debug_genes2ger {
  my ($listfile) = @_;
  my $genes = read_simple_file($listfile);
  my $gsm = get_gsm_ger();

  my %genes = map {$_, 1} @{$genes};

  foreach my $gene (sort keys %genes) {
    my $status;
    if ($gsm->contains($gene)) {
      $status = "identical";
    } elsif (my $new = $gsm->resolve("-symbol" => $gene)) {
      $status = "resolved to $new";
    } else {
      $status = "missing";
    }
    printf "status for %s: %s\n", $gene, $status;
  }
}

sub hack18 {
  my $db = get_vm_tp53_gl();
  my $hits = $db->find_aa_codons(
				 "-gene" => "TP53",
				  "-aa" => "P128L"
				 );
  foreach my $hit (@{$hits}) {
    dump_die($hit, "Debug", 1);
  }
  die;
}

sub add_batch_thousand_genomes {
  my (%options) = @_;
  my $rows = $options{"-rows"} || die "-rows";
  my $start_time = time;
  log_msg("batch 1,000 Genomes start");
  my $tabix = get_tabix_1kg();
  my $f_row_key = "_user_row";
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $twobit = $config_genome->{TWOBIT} || die "no TWOBIT";

  my @query_variants;
  foreach my $r (@{$rows}) {
    my $v = get_variant_from_row($r);
    $v->{$f_row_key} = $r;
    push @query_variants, $v;
  }

  my %map;
  $map{$FIELD_THOUSAND_GENOMES} = "AF";

  my $tba = new TabixBatchAnnotation(
				     "-tabix" => $tabix,
				     "-twobit" => $twobit,
				     "-split_count" => TABIX_BATCH_SIZE_THOUSAND_GENOMES,
				     "-f_tabix_chr" => "Chr",
				     "-f_tabix_pos" => "Pos",
				     "-f_tabix_ref_allele" => "ReferenceAllele",
				     "-f_tabix_var_allele" => "MutantAllele",

				     "-user_row_key" => $f_row_key,
				     "-annotation_map" => \%map,
				    );

  $tba->query(
	      "-query" => \@query_variants,
	     );

  log_msg(sprintf "batch 1,000 Genomes annotation for %d rows took %d",
	  scalar(@{$rows}), time - $start_time);

#  foreach my $r (@{$rows}) {
#    dump_die($r, "debug" ,1);
#  }

}

sub add_batch_dbsnp {
  my (%options) = @_;
  my $rows = $options{"-rows"} || die "-rows";
  my $start_time = time;
  log_msg("batch dbSNP start");
  my $tabix = get_tabix_dbsnp();
  my $f_row_key = "_user_row";
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $twobit = $config_genome->{TWOBIT} || die "no TWOBIT";

  my @query_variants;
  foreach my $r (@{$rows}) {
    my $v = get_variant_from_row($r);
    $v->{$f_row_key} = $r;
    push @query_variants, $v;
  }

  my %map;
  $map{$FIELD_DBSNP} = "rs";

  my $tba = new TabixBatchAnnotation(
				     "-tabix" => $tabix,
				     "-twobit" => $twobit,
				     "-split_count" => TABIX_BATCH_SIZE_DBSNP,
				     "-f_tabix_chr" => "Chr",
				     "-f_tabix_pos" => "WU_HG19_Pos",
				     "-f_tabix_ref_allele" => "ReferenceAllele",
				     "-f_tabix_var_allele" => "MutantAllele",

				     "-user_row_key" => $f_row_key,
				     "-annotation_map" => \%map,
				    );

  $tba->query(
	      "-query" => \@query_variants,
	     );

  foreach my $row (@{$rows}) {
    if ($row->{$FIELD_DBSNP}) {
      my @reasons;
      add_equivalent_tag(
			 "-row" => $row,
			 "-label" => "dbSNP",
			 "-tags" => \@reasons
			);
      push @{$row->{$INTERNAL_FIELD_QUEUE_REASONS}}, @reasons if @reasons;
    }
  }

  log_msg(sprintf "batch dbSNP annotation for %d rows took %d",
	  scalar(@{$rows}), time - $start_time);


}

sub add_batch_gedi {
  my (%options) = @_;
  my $rows = $options{"-rows"} || die "-rows";
  my $start_time = time;
  log_msg("batch GEDI start");
  my $tabix = get_tabix_gedi();
  my $f_row_key = "_user_row";
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $twobit = $config_genome->{TWOBIT} || die "no TWOBIT";

  my @query_variants;
  foreach my $r (@{$rows}) {
    my $v = get_variant_from_row($r);
    $v->{$f_row_key} = $r;
    push @query_variants, $v;
  }

  my $tba = new TabixBatchAnnotation(
				     "-tabix" => $tabix,
				     "-twobit" => $twobit,
				     "-split_count" => TABIX_BATCH_SIZE_GEDI,
				     "-f_tabix_chr" => "chromosome",
				     "-f_tabix_pos" => "pos",
				     "-f_tabix_ref_allele" => "reference_allele",
				     "-f_tabix_var_allele" => "non_reference_allele",

				     "-user_row_key" => $f_row_key,

				     "-store_hits" => $FIELD_GEDI_TABIX,

				    );

  $tba->query(
	      "-query" => \@query_variants,
	     );

  log_msg(sprintf "batch GEDI annotation for %d rows took %d",
	  scalar(@{$rows}), time - $start_time);


}

sub get_tabix_1kg {
  my $f_tabix = $FLAGS{"tabix-thousand-genomes"} || die "-tabix-thousand-genomes";
  printf STDERR "1000 genomes AF: %s\n", $f_tabix;

  return get_batch_tabix($f_tabix);
}

sub get_tabix_dbsnp {
  my $f_tabix = $FLAGS{"tabix-dbsnp"} || die "-tabix-dbsnp";
  printf STDERR "dbSNP: %s\n", $f_tabix;
  return get_batch_tabix($f_tabix);
}

sub get_tabix_gedi {
  my $f_tabix = $FLAGS{"tabix-gedi"} || die "-tabix-gedi";
  printf STDERR "GEDI: %s\n", $f_tabix;
  return get_batch_tabix($f_tabix);
}

sub get_tabix_clinvar {
  my $f_tabix = $FLAGS{"tabix-clinvar"} || die "-tabix-clinvar";
  printf STDERR "ClinVar: %s\n", $f_tabix;
  return get_batch_tabix($f_tabix);
}

sub get_tabix_taylor_hotspots {
  my $f_tabix = $FLAGS{"taylor-hotspots"} || die "-taylor-hotspots";
  printf STDERR "TaylorHotspot: %s\n", $f_tabix;
  return get_batch_tabix($f_tabix);
}

sub get_batch_tabix {
  my ($f_tabix) = @_;
  return new TabixFile(
		       "-file" => $f_tabix,
		       "-indel_wiggle_bases" => $TABIX_INDEL_EQUIVALENCE_DISTANCE
		      );
}

sub get_tabix_nhlbi {
  my $f_tabix = $FLAGS{"tabix-nhlbi"} || die "-tabix-nhlbi";
  printf STDERR "NHLBI: %s\n", $f_tabix;
  return get_batch_tabix($f_tabix);
}

sub get_tabix_cosmic {
  my $f_tabix = $FLAGS{"tabix-cosmic"} || die "-tabix-cosmic";
  printf STDERR "COSMIC: %s\n", $f_tabix;
  return get_batch_tabix($f_tabix);
}

sub thousand_genomes_hack {
  die unless @INPUT_SNV_INDEL_FILES;
  foreach my $f (@INPUT_SNV_INDEL_FILES) {
    printf STDERR "processing %s...\n", $f;
    my $df = new DelimitedFile(
			       "-file" => $f,
			       "-headers" => 1
			      );
    my @rows;
    while (my $row = $df->get_hash()) {
      push @rows, $row;
    }

    add_batch_thousand_genomes(
			       "-rows" => \@rows,
			      );
    foreach my $row (@rows) {
      printf STDERR "%s\n", $row->{$FIELD_THOUSAND_GENOMES};
    }

  }
}

sub hack_exac_vcf2tab {
  my @infiles = @INPUT_GL_FILES;
  die "no input gl files" unless @infiles;
  foreach my $f (@infiles) {
    printf STDERR "processing %s...\n", $f;
    my $df = new DelimitedFile(
			       "-file" => $f,
			       "-headers" => 1
			      );
    my @rows;
    while (my $row = $df->get_hash()) {
      push @rows, $row;
    }

    add_batch_exac (
		    "-rows" => \@rows,
		   );
  }
}

sub hack_exac_coverage {
  my @infiles = @INPUT_GL_FILES;
  die "no input gl files" unless @infiles;
  foreach my $f (@infiles) {
    printf STDERR "processing %s...\n", $f;
    my $df = new DelimitedFile(
			       "-file" => $f,
			       "-headers" => 1
			      );
    my @rows;
    while (my $row = $df->get_hash()) {
      push @rows, $row;
    }

    add_batch_exac_coverage(
			    "-rows" => \@rows,
			   );
  }
}

sub add_batch_exac {
  my (%options) = @_;
  my $rows = $options{"-rows"} || die "-rows";

  my $fn = $FLAGS{"exac-vcf2tab"} || die "-exac-vcf2tab";
  printf STDERR "ExAC AF: %s\n", $fn;

  my $tabix = get_batch_tabix($fn);

  my $f_row_key = "_user_row";
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $twobit = $config_genome->{TWOBIT} || die "no TWOBIT";

  my @tabix_fields = qw(
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
			 PubMed
		      );

  my @query_variants;
  foreach my $r (@{$rows}) {
    delete @{$r}{@tabix_fields};
    # wipe in case present for some reason (shouldn't be)
    my $v = get_variant_from_row($r);
    $v->{$f_row_key} = $r;
    push @query_variants, $v;
  }

  my %map = map {$_, $_} @tabix_fields;
  # to do: add leading _ etc. for key safety?

  #
  #  add raw annotation fields to each row:
  #
  my $tba = new TabixBatchAnnotation(
				     "-tabix" => $tabix,
				     "-twobit" => $twobit,
				     "-split_count" => TABIX_BATCH_SIZE_EXAC,
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

  #
  #  for each row, generate old-style Reason code info:
  #
  my @population_tags = qw(
			    AFR
			    AMR
			    EAS
			    FIN
			    NFE
			    SAS
			    OTH
			 );
  foreach my $r (@{$rows}) {
    if ($r->{AC}) {
      my $ac = $r->{AC};
      my $an = $r->{AN} || die;
      my $ac_adj = $r->{AC_Adj};
      my $an_adj = $r->{AN_Adj};
      my $af = $r->{AF} || dump_die($r, "no AF");
      # AF may be either decimal or exponential notation
      my @things;
      push @things, sprintf 'ExAC_AF=%s', $af;
      push @things, sprintf 'ExAC_AC=%d', $ac;
      push @things, sprintf 'ExAC_AN=%d', $an;
      push @things, sprintf 'ExAC_AC_Adj=%d', $ac_adj;
      push @things, sprintf 'ExAC_AN_Adj=%d', $an_adj;

#      dump_die($r);

      my @pops;
      foreach my $pop (@population_tags) {
	my $tag_allele = sprintf 'AC_%s', $pop;
	my $tag_chroms = sprintf 'AN_%s', $pop;

	my $count_allele = $r->{$tag_allele};
	if (defined $count_allele) {
	  my @out = $count_allele;

	  #      push @pop_ac, sprintf '%s:%d', $pop, $count_allele;

	  my $count_chrom = $r->{$tag_chroms};
	  if (defined($count_chrom)) {
	    #	my $freq = $count_allele / $count_chrom;
	    push @out, $count_chrom;
	    my $freq = $count_chrom ? ($count_allele / $count_chrom) : 0;
	    # might be 0, e.g. 11.2906165.A.G
	    push @out, $freq > 0 ? sprintf("%4.3e", $freq) : $freq;
	  }

	  push @pops, join ":", $pop, @out;
	}
      }
      push @things, sprintf "ExAC_pop=%s", join ",", @pops if @pops;

      add_equivalent_tag(
			 "-row" => $r,
			 "-label" => "ExAC",
			 "-tags" => \@things
			);

      if (my $pms = $r->{PubMed}) {
	push @{$r->{$INTERNAL_FIELD_QUEUE_EVIDENCE}}, split /,/, $pms;
      }

      delete @{$r}{@tabix_fields};
      delete $r->{TABIX_VARIANT()};
      # cleanup of temporary fields (not much help)

      push @{$r->{$INTERNAL_FIELD_QUEUE_REASONS}}, @things;
    }

  }

}

sub add_batch_exac_columnar {
  my (%options) = @_;
  my $rows = $options{"-rows"} || die "-rows";

  my $fn = $FLAGS{"exac-vcf2tab"} || die "-exac-vcf2tab";
  printf STDERR "ExAC AF: %s\n", $fn;

  my $tabix = get_batch_tabix($fn);

  my $f_row_key = "_user_row";
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $twobit = $config_genome->{TWOBIT} || die "no TWOBIT";

  my @tabix_fields = qw(
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
			 PubMed
		      );

  my @query_variants;
  foreach my $r (@{$rows}) {
    delete @{$r}{@tabix_fields};
    # wipe in case present for some reason (shouldn't be)
    my $v = get_variant_from_row($r);
    $v->{$f_row_key} = $r;
    push @query_variants, $v;
  }

  my %map = map {$_, $_} @tabix_fields;
  # to do: add leading _ etc. for key safety?

  #
  #  add raw annotation fields to each row:
  #
  my $tba = new TabixBatchAnnotation(
				     "-tabix" => $tabix,
				     "-twobit" => $twobit,
				     "-split_count" => TABIX_BATCH_SIZE_EXAC,
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

  #
  #  for each row, generate old-style Reason code info:
  #
  my @population_tags = qw(
			    AFR
			    AMR
			    EAS
			    FIN
			    NFE
			    SAS
			    OTH
			 );
  foreach my $r (@{$rows}) {
    foreach my $h (@EXAC_HEADERS) {
      $r->{$h} = "";
    }
    if ($r->{AC}) {
      foreach my $f (@tabix_fields) {
	next if $f eq "PubMed";
	my $v = $r->{$f};
	$v = "" unless defined $v;
	$r->{"ExAC_" . $f} = $v;
      }

      my @reasons;
      add_equivalent_tag(
			 "-row" => $r,
			 "-label" => "ExAC",
			 "-tags" => \@reasons
			);
      push @{$r->{$INTERNAL_FIELD_QUEUE_REASONS}}, @reasons if @reasons;

      if (my $pms = $r->{PubMed}) {
	push @{$r->{$INTERNAL_FIELD_QUEUE_EVIDENCE}}, split /,/, $pms;
      }

      delete @{$r}{@tabix_fields};
      delete $r->{TABIX_VARIANT()};
      # cleanup of temporary fields (not much help)
    }
  }
}


sub add_batch_exac_coverage {
  my (%options) = @_;
  my $rows = $options{"-rows"} || die "-rows";
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $twobit = $config_genome->{TWOBIT} || die "no TWOBIT";

  my $fn = $FLAGS{"exac-coverage"} || die "-exac-coverage";
  printf STDERR "ExAC coverage: %s\n", $fn;
  my $tabix = new TabixFile(
			    "-file" => $fn,
			    "-indel_wiggle_bases" => 0,
			   );
  my $f_row_key = "_user_row";

  my @tabix_fields = qw(
			 mean
			 median
			 1
			 5
			 10
			 15
			 20
			 25
			 30
			 50
			 100
		      );

  my @query_variants;
  foreach my $r (@{$rows}) {
    delete @{$r}{@tabix_fields};
    # wipe in case present for some reason (shouldn't be)
    my $v = get_variant_from_row($r);
    $v->{$f_row_key} = $r;
    push @query_variants, $v;
  }

  my %map = map {$_, $_} @tabix_fields;
  # to do: add leading _ etc. for key safety?

  #
  #  add raw annotation fields to each row:
  #
  my $tba = new TabixBatchAnnotation(
				     "-tabix" => $tabix,
				     "-twobit" => $twobit,
				     "-split_count" => TABIX_BATCH_SIZE_EXAC_COVERAGE,
				     "-f_tabix_chr" => "chrom",
				     "-f_tabix_pos" => "pos",
				     "-user_row_key" => $f_row_key,
				     "-annotation_map" => \%map,
				     "-site_only" => 1,
				    );

  $tba->query(
	      "-query" => \@query_variants,
	     );

  #
  #  for each row, generate old-style Reason code info:
  #
  foreach my $r (@{$rows}) {
    my @things;
    if ($r->{$tabix_fields[0]} =~ /\w/) {
      foreach my $f (@tabix_fields) {
	next if $f =~ /^\d+$/ and $r->{$f} == 0;
	push @things, sprintf 'ExAC_cvg_%s=%s', $f, $r->{$f};
      }
    } else {
      push @things, "ExAC_cvg=none";
    }
    my $formatted = join ";", @things;
    $r->{$FIELD_EXAC_COVERAGE} = $formatted;
#    dump_die($r, "debug", 1);
  }

}


sub add_batch_exac_coverage_columnar {
  my (%options) = @_;
  my $rows = $options{"-rows"} || die "-rows";
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $twobit = $config_genome->{TWOBIT} || die "no TWOBIT";

  my $fn = $FLAGS{"exac-coverage"} || die "-exac-coverage";
  printf STDERR "ExAC coverage: %s\n", $fn;
  my $tabix = new TabixFile(
			    "-file" => $fn,
			    "-indel_wiggle_bases" => 0,
			   );
  my $f_row_key = "_user_row";

  my @query_variants;
  foreach my $r (@{$rows}) {
    my $v = get_variant_from_row($r);
    $v->{$f_row_key} = $r;
    push @query_variants, $v;
  }

  my %map = map {"ExAC_cvg_" . $_, $_} @EXAC_COVERAGE_H_TABIX;

  #
  #  add raw annotation fields to each row:
  #
  my $tba = new TabixBatchAnnotation(
				     "-tabix" => $tabix,
				     "-twobit" => $twobit,
				     "-split_count" => TABIX_BATCH_SIZE_EXAC_COVERAGE,
				     "-f_tabix_chr" => "chrom",
				     "-f_tabix_pos" => "pos",
				     "-user_row_key" => $f_row_key,
				     "-annotation_map" => \%map,
				     "-site_only" => 1,
				    );

  $tba->query(
	      "-query" => \@query_variants,
	     );

  #
  #  for each row, generate old-style Reason code info:
  #
  foreach my $r (@{$rows}) {
    foreach my $h (@EXAC_COVERAGE_H) {
      # make sure each output field has a blank value unless populated
      $r->{$h} = "" unless defined $r->{$h};
    }
  }

}


sub cosmic_truncating_check {
  my ($rows) = @_;

  my $STATUS_NO = 0;
  my $STATUS_MAYBE = 1;
  my $STATUS_YES = 2;

  my %wanted_desc = (
		     "Substitution - Missense" => $STATUS_NO,
		     "Unknown" => $STATUS_MAYBE,
		     "Complex" => $STATUS_MAYBE,
		     # for unknown/complex, AA annotation may
		     # reveal frameshift status

		     "Substitution - coding silent" => $STATUS_NO,

		     "Substitution - Nonsense" => $STATUS_YES,
		     "Deletion - Frameshift" => $STATUS_YES,
		     "Insertion - Frameshift" => $STATUS_YES,
		     "Complex - frameshift" => $STATUS_YES,

		     "Nonstop extension" => $STATUS_YES,
		     # e.g. p.*2844E (run-on)

		     "Whole gene deletion" => $STATUS_NO,

		     "Deletion - In frame" => $STATUS_NO,
		     "Insertion - In frame" => $STATUS_NO,
		     "Complex - deletion inframe" => $STATUS_NO,
		     "Complex - insertion inframe" => $STATUS_NO,
		     # QUESTION: what about inframe indels???
		     "Complex - compound substitution" => $STATUS_NO,
		     # are these always inframe?

		     "No detectable mRNA/protein" => $STATUS_NO

		    );
  # variant categories we definitely want to include

  my $any_usable;

  foreach my $row (@{$rows}) {
    my $aa = $row->{"Mutation AA"};
    my $desc = $row->{"Mutation Description"};

    my $wanted_desc = $wanted_desc{$desc};
    die "unhandled COSMIC Mutation Description \"$desc\" aa=$aa" unless defined $wanted_desc;

    my $usable;
    if ($wanted_desc == $STATUS_YES) {
      $usable = 1;
    } elsif ($wanted_desc == $STATUS_MAYBE) {
      if ($aa =~ /fs/) {
	# desc=Unknown aa=p.L35fs*10
	# desc=Complex aa=p.R158fs
	# desc=Complex aa=p.R158fs
	# desc=Unknown aa=p.L35fs*10
	# desc=Complex aa=p.R158fs
	# desc=Complex aa=p.R158fs
	# desc=Unknown aa=p.L35fs*10
	# desc=Complex aa=p.R26fs
	# desc=Complex aa=p.R26fs
	# desc=Complex aa=p.R65fs
	# desc=Complex aa=p.R65fs
	$usable = 1;
      }
    }

    if ($usable) {
      $any_usable = 1;
      last;
    }
  }

  return $any_usable;
}

sub add_batch_clinvar {
  my (%options) = @_;
  my $rows = $options{"-rows"} || die "-rows";
  my $start_time = time;
  my $label = "ClinVar";

  log_msg("batch $label start");
  my $vm_clinvar_aa = get_vm_clinvar_reviewable_pathogenic_aa();
  if (0) {
    print STDERR "DEBUG: disable AA\n";
    $vm_clinvar_aa = new VariantMatcher();
  }

  my $tabix = get_tabix_clinvar();
  my $f_row_key = "_user_row";
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $twobit = $config_genome->{TWOBIT} || die "no TWOBIT";

  my @query_variants;
  foreach my $r (@{$rows}) {
    my $v = get_variant_from_row($r);
    $v->{$f_row_key} = $r;
    push @query_variants, $v;
  }

  my $tba = new TabixBatchAnnotation(
				     "-tabix" => $tabix,
				     "-twobit" => $twobit,
				     "-split_count" => TABIX_BATCH_SIZE_CLINVAR,
				     "-f_tabix_chr" => "Chr",
				     "-f_tabix_pos" => "WU_HG19_Pos",
				     "-f_tabix_ref_allele" => "ReferenceAllele",
				     "-f_tabix_var_allele" => "MutantAllele",

				     "-user_row_key" => $f_row_key,

				     "-store_hits" => $FIELD_CLINVAR_TABIX,
				     "-store_site" => $FIELD_CLINVAR_TABIX_SITE,
				    );

  $tba->query(
	      "-query" => \@query_variants,
	      "-callback" => [ \&batch_clinvar_medal, "-vm" => $vm_clinvar_aa ],
	     );

  log_msg(sprintf "batch %s annotation for %d rows took %d",
	  $label, scalar(@{$rows}), time - $start_time);

}

sub add_batch_taylor_hotspots {
  my (%options) = @_;
  my $rows = $options{"-rows"} || die "-rows";
  my $start_time = time;
  my $label = "TaylorHotspot";

  log_msg("batch $label start");
  my $tabix = get_tabix_taylor_hotspots();
  my $f_row_key = "_user_row";
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $twobit = $config_genome->{TWOBIT} || die "no TWOBIT";

  my @query_variants;
  foreach my $r (@{$rows}) {
    my $v = get_variant_from_row($r);
    $v->{$f_row_key} = $r;
    push @query_variants, $v;
  }

  my $tba = new TabixBatchAnnotation(
				     "-tabix" => $tabix,
				     "-twobit" => $twobit,
				     "-split_count" => TABIX_BATCH_SIZE_CLINVAR,
				     "-f_tabix_chr" => "Chr",
				     "-f_tabix_pos" => "Pos",
				     "-f_tabix_ref_allele" => "Chr_Allele",
				     "-f_tabix_var_allele" => "Alternative_Allele",

				     "-user_row_key" => $f_row_key,

				     "-store_hits" => $FIELD_TAYLOR_TABIX,
				    );

  $tba->query(
	      "-query" => \@query_variants,
	      "-callback" => \&batch_taylor_hotspot_medal,
	     );

  log_msg(sprintf "batch %s annotation for %d rows took %d",
	  $label, scalar(@{$rows}), time - $start_time);

}

sub get_is_rare {
  my ($row, $return_freq) = @_;
#  my $nhlbi_freq = $row->{$FIELD_NHLBI_FREQ} || 0;
#  my $freq = $row->{$FIELD_NHLBI_FREQ} || 0;
  my $freq = $row->{ExAC_AF};
  # 1/2018: switch from NHLBI to ExAC
  my $dbname = "ExAC_ex_TCGA";

  unless (defined $freq and length($freq)) {
    # non-empty
    $freq = $row->{$FIELD_THOUSAND_GENOMES};
    if ($freq) {
      # in some rare cases a variant will not be present in ExAC,
      # making it seem rare, however it will appear at high
      # frequency in 1,000 genomes, e.g. for hg19:
      #   chr4.1795557.T.C
      #   chr22.19747082.G.C
      # Note both of these are splice_region variants and ExAC might
      # not cover them anyway.
      # In any event 1KG lets us catch a few more high-frequency variants
      # ExAC misses.  May become moot once we switch to gnomad.
      $dbname = "1KG";
    } else {
      $freq = 0;
    }

    dump_die($row, "found in NHLBI but not ExAC", 1) if $row->{$FIELD_NHLBI_FREQ};

  }

  my $is_rare = $freq > $GERMLINE_MAX_NHLBI_FREQ ? 0 : 1;
  if ($return_freq) {
    return ($freq, $dbname);
  } else {
    return $is_rare;
  }
}

sub get_is_silent_intron_utr {
  my ($row) = @_;
  my $variant_class = lc($row->{Class} || "");
  my $is_silent_intron_utr = 0;
  if ($variant_class eq "silent" or
      $variant_class eq "intron" or
      $variant_class eq "utr_5" or
      $variant_class eq "utr_3") {
    $is_silent_intron_utr = 1;
  }
  return $is_silent_intron_utr;
}

sub batch_clinvar_medal {
  my (%options) = @_;
  my $rows = $options{"-rows"} || die "-rows";
  my $vm_aa = $options{"-vm"} || die "-vm";

  my $VERBOSE = $ENV{VERBOSE_CLINVAR};

  foreach my $row (@{$rows}) {
    my @reasons;
    my $medal_expert;

    my $aa = $row->{$FIELD_AACHANGE};
    my $gene = $row->{$FIELD_GENE};
    my $nm = $row->{$FIELD_REFSEQ};

    if (clinvar_usable_hits($row, $FIELD_CLINVAR_TABIX)) {
      printf STDERR "already have clinvar match\n" if $VERBOSE;
    } else {
      # if we don't already have a perfect match, try amino acid lookup
      if ($gene and $aa and $nm) {
	printf STDERR "ClinVar AA lookup %s %s\n", $gene, $aa if $VERBOSE;
	my $hits = $vm_aa->find_aa(
				  "-gene" => $gene,
				  "-aa" => $aa,
				 );
	if ($hits and
	    $vm_aa->match_type == VariantMatcher::MATCH_TYPE_PERFECT) {
	  my $handshake_ok = check_transcript_handshake(
							"-row" => $row,
							"-transcript-handshake" => $FIELD_REFSEQ,
							"-hits" => $hits
						       );
	  if ($handshake_ok) {
	    # replace existing result with AA hit so below logic can
	    # be used unchanged
	    $row->{$FIELD_CLINVAR_TABIX} = $hits;
	  }
	}
      }
    }

    my $found_ACMG_PS1;

    if (clinvar_usable_hits($row, $FIELD_CLINVAR_TABIX)) {
      # perfect match
      push @reasons, "ClinVar";
      foreach my $r (@{$row->{$FIELD_CLINVAR_TABIX}}) {
	if (my $cvid = $r->{ClinVar_Variation_ID}) {
	  push @reasons, sprintf 'ClinVar_Variation_ID=%s', $cvid;
	} else {
	  dump_die($r, "ERROR: where is ClinVar_Variation_ID?", 1);
	  # format change??
	}

	my $is_p;
	foreach my $clnsig (split /,/, $r->{CLNSIG}) {
	  $is_p = 1 if $clnsig == CLINVAR_CLNSIG_P;
	}
	if ($is_p) {
	  push @reasons, generate_acmg_reason(ACMG_2015, "PS1", "ClinVar");
	  $found_ACMG_PS1 = 1;
	  # TO DO:
	  # - specify exact variant?
	  # - only do this if trusted curation, i.e. below?
	}

	# check review curation status:
	my $clinvar_code = $r->{CLNREVSTAT};
	dump_die($r, "where is CLNREVSTAT") unless $clinvar_code;
	my $gold_stars = $CLINVAR_CLNREVSTAT_TO_GOLD_STARS{$clinvar_code};
	if (defined $gold_stars) {
	  if ($gold_stars >= $CLINVAR_MIN_GOLD_STARS_TO_TRUST) {
	    my $clnsigs = $r->{CLNSIG};
	    dump_die($r, "where is CLNSIG") unless defined $clnsigs;
	    push @reasons, "ClinVar_trusted_curation";

	    my %all_clnsig = map {$_, 1} (split /,/, $clnsigs);

	    if ($all_clnsig{CLINVAR_CLNSIG_P()} or
		$all_clnsig{CLINVAR_CLNSIG_LP()}) {
	      $medal_expert = CLASS_GOLD;
	      $row->{$INTERNAL_FIELD_NO_POPULATION_FILTER} = 1;
	      # exempt trusted ClinVar P/LP variants from population filtering.
	      # some founder variants are at higher frequencies in the
	      # population, yet can be pathogenic if monoallelic.
	      # e.g. chr1.45797228.C.T see threads w/Zhaoming and Scott,
	      # 9/2018.
	    } elsif ($all_clnsig{CLINVAR_CLNSIG_U()}) {
	      $medal_expert = CLASS_SILVER;
	      # lean towards keeping population filtering policy in place;
	      # unknown status doesn't seem like a compelling enough
	      # reason to lift it.
	    } elsif ($all_clnsig{CLINVAR_CLNSIG_LB()} or
		     $all_clnsig{CLINVAR_CLNSIG_B()}) {
	      # note that B/LB calls are still subject to population
	      # frequency filtering, e.g. LIG4.A3V (ClinVar ID 7676)
	      # is a B/LB high-frequency variant.
	      $medal_expert = CLASS_BRONZE;
	    }
	    # miscellaneous types (e.g. drug response) will get a
	    # standard silver (below)
	  }
	} else {
	  # shouldn't happen; could big a major problem like a big file
	  # format change (which ClinVar has done before)
	  dump_die($row, "ERROR: unknown ClinVar CLNREVSTAT value $clinvar_code, format change in VCF distribution??");
	}

	if (my $ids = $r->{PubMed}) {
	  # pass through PubMed IDs if present
	  # TO DO: might also be links for PubMedCentral and NCBIBookShelf
	  push @{$row->{$INTERNAL_FIELD_QUEUE_EVIDENCE}}, split /,/, $ids;
	}
      }
    } elsif (clinvar_usable_hits($row, $FIELD_CLINVAR_TABIX_SITE)) {
      # site-only match
      push @reasons, "ClinVar_site_only";
    }

    if (not($found_ACMG_PS1) and
	lc($row->{$FIELD_CLASS}) eq "missense" and $aa) {
      # check for matches to a different AA change in the same codon
      # that's pathogenic.  Does't quite fit into code above though
      # since sometimes we have a perfect ClinVar match to a VUS
      # (e.g. TP53 R337G), however this is also the site of a pathogenic
      # variant (R337H).
      my $hits = $vm_aa->find_aa_specific_codon(
						"-gene" => $gene,
						"-aa" => $aa
					       );
      my @cvid;
      foreach my $h (@{$hits}) {
	next unless lc($h->{$FIELD_CLASS}) eq "missense";
	# "A novel missense amino acid change occurring at the same
	# position as another pathogenic missense change (e.g.,
	# Trp38Ser and Trp38Leu) is considered moderate evidence"
	#
	# => want to match to missense specifically, ignoring e.g.
	#    frameshifts which are inherently pathogenic
	push @cvid, $h->{ClinVar_Variation_ID};
#	dump_die($h, "debug", 1);
      }

      push @reasons, generate_acmg_reason(ACMG_2015, "PM5", "ClinVar", @cvid) if @cvid;
    }

    if (@reasons) {
      my $is_rare = get_is_rare($row);
      my $is_silent_intron_utr = get_is_silent_intron_utr($row);
      my $medal;
      if ($medal_expert) {
	# medal based on expert/trusted ClinVar classification
	$medal = $medal_expert;
      } elsif ($is_rare and not($is_silent_intron_utr)) {
	# this filtering still applies if review status is < 2 stars
	$medal = CLASS_SILVER;
	# both perfect/imperfect get a silver
      }
      push @{$row->{$INTERNAL_FIELD_QUEUE_MEDALS}}, $medal if $medal;
      push @{$row->{$INTERNAL_FIELD_QUEUE_REASONS}}, @reasons;
    }

    delete $row->{$FIELD_CLINVAR_TABIX};
    delete $row->{$FIELD_CLINVAR_TABIX_SITE};
    # release tabix resources

#    dump_die($row);
  }
}

sub batch_taylor_hotspot_medal {
  my (%options) = @_;
  my $rows = $options{"-rows"} || die "-rows";
  foreach my $row (@{$rows}) {
    my ($reason, $medal);

    if (my $hits = $row->{$FIELD_TAYLOR_TABIX}) {
      $medal = CLASS_SILVER;
      $reason = "taylor_hotspot";
    }

    push @{$row->{$INTERNAL_FIELD_QUEUE_MEDALS}}, $medal if $medal;
    push @{$row->{$INTERNAL_FIELD_QUEUE_REASONS}}, $reason if $reason;

    delete $row->{$FIELD_TAYLOR_TABIX};
    # release tabix resources

#    dump_die($row);
  }
}


sub clinvar_usable_hits {
  # clinvar batch tabix: do the result rows have a usable CLNSIG value?
  my ($row, $store_key) = @_;
  die unless defined $store_key;
  my $usable = 0;

  my %wanted = map {$_, 1} GERMLINE_CLINVAR_CLNSIG_TYPES;

  if (my $hits = $row->{$store_key}) {
    foreach my $h (@{$hits}) {
      my $list = $h->{CLNSIG};
      dump_die($h, "WTF: no CLNSIG") unless defined $list;
      foreach my $v (split /,/, $list) {
	$usable = 1 if $wanted{$v};
      }
    }
  }
  return $usable;
}

sub dump_config_entry {
  my (%options) = @_;
  my $config_key = $options{"-key"} || die "-key";
  printf "%s\n", $config_key;
}

sub add_equivalent_tag {
  my (%options) = @_;
  my $r = $options{"-row"} || die "-row";
  my $label = $options{"-label"} || die "-label";
  my $tags = $options{"-tags"} || die "-tags";
  my $query_snv4 = get_variant_from_row($r)->get_snv4();
  my $match_snv4 = $r->{TABIX_VARIANT()}->get_snv4();
  push @{$tags}, sprintf "%s_equivalent=%s", $label, $match_snv4 if $match_snv4 ne $query_snv4;
  # add tag indicating tabix equivalent variant site if query != match
  # (i.e. equivalent indel match)
}

sub ram_debug {
  my ($label) = @_;
  if ($FLAGS{"debug-ram"}) {
    my $cmd = sprintf 'ps u';
    open(PSTMP, sprintf '/bin/ps u %d|', $$) || die;
    my $hl = <PSTMP>;
    chomp $hl;
    my $dl = <PSTMP>;
    close PSTMP;

    my @h = split /\s+/, $hl;
    my @d = split /\s+/, $dl;
    # HACK: breaks for command portion, but should work for earlier fields
    my %info;
    @info{@h} = @d;

    log_msg(sprintf "RAM at %s: RSS:%d VSZ:%d", $label, @info{qw(RSS VSZ)});
  }
}


sub get_cosmic_snv4 {
  my ($row) = @_;
  my $chr = $row->{Chr} || die;
  my $pos = $row->{WU_HG19_Pos} || die;
  my $ra = $row->{ReferenceAllele} || die;
  my $va = $row->{MutantAllele} || die;
  return join ".", $chr, $pos, $ra, $va;
}

sub get_vo_promoters {
  my $infile = $FLAGS{"known-promoter-sites"} || die "-known-promoter-sites";
  my $vo = new VariantOverlap();
  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  while (my $row = $df->get_hash()) {
    my $chr = $row->{Chr} || die;
    my $pos = $row->{Pos} || die;
    $vo->add_site("-reference" => $chr, "-position" => $pos);
  }
  return $vo;
}

sub promoter_check {
  my (%options) = @_;
  my $row = $options{"-row"} || die;
  my $medal_ref = $options{"-medal-ref"} || die;
  my $reasons_ref = $options{"-reasons-ref"} || die;

  if ($ENABLE_KNOWN_PROMOTER_SITES) {
    my $v = get_variant_from_row($row);

    $main::VO_PROMOTERS = get_vo_promoters() unless $main::VO_PROMOTERS;

    if ($main::VO_PROMOTERS->overlaps("-variant" => $v)) {
      $$medal_ref = request_medal($$medal_ref, CLASS_GOLD);
      unshift @{$reasons_ref}, "known_promoter_site";
    }
  }

  if ($ENABLE_KNOWN_PROMOTER_REGIONS) {
    unless ($main::CBM_PR) {
      my $glr = read_simple_file($FLAGS{"gl-reviewable-genes-refflat"} || die);
      die "convert to main reviewable list and add GSM disambiguation";
      # never in production, needs updates so we can obsolete
      # gl-reviewable-genes-refflat
      my $f_gl_custom = $FLAGS{"gl-custom-genes"};
      printf STDERR "WARNING: custom genes not implemented with promoter intervals\n" if $f_gl_custom;

      my $infile = $FLAGS{"known-promoter-intervals"} || die "-known-promoter-intervals";
      my $df = new DelimitedFile("-file" => $infile,
				 "-headers" => 1,
				);
      $main::CBM_PR = new ChrBucketMap(
				       "-f_chr" => "chrom",
				       "-f_start" => "txStart",
				       "-f_end" => "txEnd",
				      );

      my %glr = map {$_, 1} @{$glr};
      my %p_needed = %glr;
      while (my $row = $df->get_hash()) {
	my $gene = $row->{name2} || die;
	if ($glr{$gene}) {
	  delete $p_needed{$gene};
	  $row->{txStart}++;
	  # convert from UCSC to 1-based
	  $main::CBM_PR->add_row("-row" => $row);
	}
      }

      delete $p_needed{STL};
      # won't be present (RNF217-AS1)
      if (%p_needed) {
	foreach my $g (sort keys %p_needed) {
	  printf STDERR "ERROR: can't find promoter intervals for \"%s\"\n", $g;
	}
	die "missing promoter intervals";
      }
    }

    #
    #  check:
    #
    my $in_promoter_region;
    my $v = new Variant();
    $v->import_bambino_row("-row" => $row, "-postprocessed" => 1);
    if (my $hits = $main::CBM_PR->find(
			     "-chr" => $v->reference_name,
			     "-start" => $v->start,
			     "-end" => $v->end,
			    )) {
      $in_promoter_region = $hits->[0]->{name2};
    }

    my $variant_class = lc($row->{Class} || "");
    my $is_rare = get_is_rare($row);

    if ($is_rare and
	($variant_class eq "promoter_region" or $in_promoter_region)) {
      unshift @{$reasons_ref}, sprintf "promoter_region=%s", $in_promoter_region || "unknown";
      $$medal_ref = request_medal($$medal_ref, CLASS_BRONZE);
    }
  }
}


sub config_somatic_cnv_add_manual_genes {
  my $f_gm = $FLAGS{"genes-manual"} || die "-genes-manual";
  my $f_cnv = $FLAGS{"cnv"} || die "-cnv";
  printf STDERR "CNV: %s\n", $f_cnv;

  my $f_cnv_out = basename($f_cnv) . ".patched.txt";
  my $fh_in = new FileHandle();
  $fh_in->open($f_cnv) || die;
  my $wf_out = new WorkingFile($f_cnv_out);
  my $fh_out = $wf_out->output_filehandle();
  while (<$fh_in>) {
    print $fh_out $_;
  }
  my $df = new DelimitedFile(
			  "-file" => $f_gm,
			  "-headers" => 1
			 );
  while (my $row = $df->get_hash()) {
    my $gene = $row->{Gene} || die;
    my $class = $row->{Class} || die;
    my ($is_ts, $is_recurrent) = parse_gene_class($class);
#    die "$gene $class $is_ts $is_recurrent";
    my $annot = "not_available";
    die unless $is_ts or $is_recurrent;

    printf $fh_out "%s\n", join "\t", $gene, "Del", $annot if $is_ts;
    printf $fh_out "%s\n", join "\t", $gene, "Amp", $annot if $is_recurrent;
  }
  $wf_out->finish();
}

sub committee_check {
  my (%options) = @_;
  my $row = $options{"-row"} || die "-row";
  my $is_silent = $options{"-is-silent"};
  die "-is-silent" unless defined $is_silent;

  my $panel_decision = "";
  my $committee_call_medal = "";
  my $request;
  my ($hits, $perfect_hit) = assign_aa_match(
					  %options,
					  "-silent-allowed" => 1,
					  "-protect" => 1,
					  "-transcript-handshake" => $FIELD_REFSEQ,
					  "-label" => "committee_reviewed",
					  "-no-medal" => 1,
					  "-genomic-lookup" => 1,
					  "-return-hits-and-perfect" => 1
					 );

  if (@{$hits}) {
    # - multiple matches are possible, because the variant may
    #   have been observed in multiple samples
    # - only populate committee field for perfect matches
    # - CHECK: INDEL WIGGLE ROOM?
    #	  dump_die($row, "debug, hit to:", 1);
    if ($perfect_hit) {
      my %medals;
      foreach my $hit (@{$hits}) {
	#	      dump_die($hit, "debug hit:", 1);
	my $decision = $hit->{$FIELD_PANEL_DECISION} || die;
	my $m;
	if ($decision =~ /gold/i) {
	  $m = CLASS_GOLD;
	} elsif ($decision =~ /silver/i) {
	  $m = CLASS_SILVER;
	} elsif ($decision =~ /bronze/i) {
	  $m = CLASS_BRONZE;
	} elsif ($decision =~ /^([A-Z]+)/) {
	  my $code = $1;
	  if ($code eq CLASS5_BENIGN or
	      $code eq CLASS5_LIKELY_BENIGN or
	      $code eq CLASS5_UNKNOWN or
	      $code eq CLASS5_LIKELY_PATHOLOGIC or
	      $code eq CLASS5_PATHOLOGIC) {
	    $m = $code;
	  } else {
	    die "can't parse committee code $code";
	  }
	} else {
	  die "can't parse committee call: $decision";
	}
	$medals{$m}++;
      }

      my %weight = %MEDAL_SORT_WEIGHT;
      $weight{CLASS_UNKNOWN()} = 150;
      $weight{CLASS5_UNKNOWN()} = 150;
      # weighting used for output file puts unknown at the bottom.
      # however here unknown should be better than LB but worse than LP.

      my @sorted = sort {($weight{$b} || die) <=> ($weight{$a} || die)} keys %medals;
      dump_die($row, sprintf("WARNING: committee medal conflict at %s: %s",
			     #				   join(".", map {$row->{$_}} qw(Chr WU_HG19_Pos ReferenceAllele MutantAllele)),
			     get_gedi_key($row),
			     join(",", sort keys %medals)), 1) if keys %medals > 1;

      #	    ($panel_decision) = (keys %medals);
      ($panel_decision) = $sorted[0];
      die unless $panel_decision;

      if ($panel_decision eq CLASS5_PATHOLOGIC or
	  $panel_decision eq CLASS5_LIKELY_PATHOLOGIC) {
	$request = CLASS_GOLD;
      } elsif ($panel_decision eq CLASS5_UNKNOWN) {
	$request = CLASS_SILVER;
      } elsif ($panel_decision eq CLASS5_LIKELY_BENIGN or
	       $panel_decision eq CLASS5_BENIGN) {
	$request = CLASS_BRONZE;
      } else {
	die "unknown panel call $panel_decision";
      }
      $committee_call_medal = $request;
    } elsif (not($is_silent)) {
      #
      #  imperfect match to committee variant, and not silent:
      #
      $request = CLASS_SILVER;
    }
  }

  return ($request, $panel_decision, $committee_call_medal);
}

sub generate_acmg_reason {
  # TO DO: discuss/finalize formatting?
  my ($spec_name, $tag_name, @info) = @_;
  return sprintf '%s=%s', TAG_CLASSIFY, join ",", $spec_name, $tag_name, @info;
}

sub get_acmg_pvs1_caution_regions {
 # "One must also be cautious when interpreting truncating
 # variants downstream of the most 3′ truncating variant
 # established as pathogenic in the literature. This is especially
 # true if the predicted stop codon occurs in the last
 # exon or in the last 50 base pairs of the penultimate exon..."
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $f_refflat = $config_genome->{REFSEQ_REFFLAT} || die "no genome config REFSEQ_REFFLAT";

  my $rf = new RefFlatFile("-strip_sharp_annotations" => 1);
  $rf->parse_file(
		  "-refflat" => $f_refflat,
		  "-type" => "refgene",
		 );

  my $rows = $rf->rows();

  my %track;
  foreach my $ref (@{$rows}) {
    my $strand = $ref->{strand} || die;
    my $cds_end;
    if ($strand eq "+") {
      $cds_end = $ref->{cdsEnd} || die;
    } elsif ($strand eq "-") {
      # don't ask  :P
      $cds_end = $ref->{cdsStart} + 1;
    } else {
      die;
    }
    next if $ref->{cdsStart} == $ref->{cdsEnd};
    # non-coding

    my $exons = $ref->{exons};
    next if @{$exons} == 1;
    # only 1 exon, e.g. NM_001004459

    my $idx_cds_end = 0;
    dump_die($ref, "1 exon") if @{$exons} == 1;
    my ($warn_start, $warn_end);

    for ($idx_cds_end = 0; $idx_cds_end < @{$exons}; $idx_cds_end++) {
      my $ex_start = $exons->[$idx_cds_end]->{start};
      my $ex_end = $exons->[$idx_cds_end]->{end};

      #    printf STDERR "exon:%d %d\n", $ex_start, $ex_end;
      if ($cds_end >= $ex_start and $cds_end <= $ex_end) {
	$warn_start = $ex_start;
	$warn_end = $ex_end;
	# start and end of last exon
	last;
      }
    }
    dump_die($ref, "can't find exon with coding end $cds_end") unless $warn_start and $warn_end;

    my $bases = 0;
    if ($strand eq "+") {
      # transcript on +, exons in array order
      my $pex = $exons->[$idx_cds_end - 1];
      my $exon_start = $pex->{start};
      for (my $pos = $pex->{end}; $pos >= $exon_start; $pos--) {
	if (++$bases > ACMG_PVS1_PENULTIMATE_EXON_WARN_BASES) {
	  #	  die "stop";
	  last;
	}
	$warn_start = $pos;
	#	printf STDERR "warn %d\n", $pos;
      }
    } else {
      # transcript on -, exons in reverse array order
      #    dump_die($ref,"debug",1);
      #    printf STDERR "warn last exon: $warn_start $warn_end $idx_cds_end\n";

      my $i_prev = $idx_cds_end + 1;
      if ($i_prev >= @{$exons}) {
	# e.g. NM_003806: stop is in first exon
	next;
	# if only 1 exon is coding, seems reasonable to warn of truncation,
	# so don't record start/stop
      } else {
	my $pex = $exons->[$idx_cds_end + 1] || die "can't move to previous exon";

	my $exon_end = $pex->{start};
	my $exon_start = $pex->{end} || dump_die($pex, "no end");
	# reversed due to - strand

	for (my $pos = $exon_end; $pos <= $exon_start; $pos++) {
	  if (++$bases > ACMG_PVS1_PENULTIMATE_EXON_WARN_BASES) {
	    #	  die "stop";
	    last;
	  }
	  $warn_end = $pos;
	  #	printf STDERR "warn %s\n", join " ", $bases, $warn_start, $warn_end;
	}
      }
    }

    die unless defined $warn_start and defined $warn_end;
    #  printf STDERR "final %s %s:%d-%d\n", $ref->{name}, $ref->{chrom}, $warn_start, $warn_end;

    my $tid = $ref->{name} || die;
    my $chr = cook_chromosome_name($ref->{chrom} || die);

    push @{$track{$tid}}, [ $chr, $warn_start, $warn_end ];
  }

  return \%track;

}

sub get_chr {
  my ($row) = @_;
  return cook_chromosome_name($row->{$FIELD_CHR} || die);
}

sub get_vm_clinvar_reviewable_pathogenic_aa {
  # clinvar pathogenic variants in reviewable genes, AA matches only
  my $gsm_reviewable = new_gsm_lite();
  my $gl_reviewable = get_germline_reviewable_genes("-original" => 1);
  foreach my $gene (@{$gl_reviewable}) {
    $gsm_reviewable->add_gene("-gene" => $gene);
  }

  my $f_clinvar = $FLAGS{"tabix-clinvar"} || die "-tabix-clinvar";

  my $df = new DelimitedFileHP(
			    "-file" => $f_clinvar,
			   );

  my $vm = new VariantMatcher();

  my @f_clone = qw(
		    mRNA_acc
		    CLNREVSTAT
		    CLNSIG
		    ClinVar_Variation_ID
		    PubMed
		    Class
		);

  ram_debug("clinvar_P_AA start");

  my $format_ok = 1;
  foreach my $h ("GeneName",
		 "AAChange",
		 @f_clone) {
    $format_ok = 0 unless $df->has_header($h);
  }

  if ($format_ok) {
    while ($df->next_row()) {
      my $gene = $df->get_value("GeneName");

      next unless $gsm_reviewable->find($gene);
      # filter to only germline-reviewable genes

      my $clnsigs = $df->get_value("CLNSIG");
      my $wanted;
      foreach my $clnsig (split /,/, $clnsigs) {
	$wanted = 1 if $clnsig == CLINVAR_CLNSIG_P;
	# only variants marked pathogenic
	# TO DO: also require expert curation??
      }
      next unless $wanted;

      my $aa = $df->get_value("AAChange") || next;
      my %r = map {$_, $df->get_value($_)} @f_clone;
      printf STDERR "loaded ClinVar %s %s\n", $gene, $aa if $ENV{VERBOSE_CLINVAR};

      #    printf STDERR "clinvar P load %s %s\n", $gene, $aa;
      $vm->add_aa(
		  "-aa" => $aa,
		  "-gene" => $gene,
		  # get rid of this an just use transcript?
		  "-row" => \%r
		 );

    }
  } else {
    printf STDERR "ERROR: ClinVar file doesn't contain gene annotations; AA lookup disabled.\n";
  }

  ram_debug("clinvar_P_AA end");

  return $vm;

}

sub assign_acmg_pm1 {
  my %options = @_;
  my $reasons = $options{"-reasons"} || die;
  push @{$reasons}, generate_acmg_reason(ACMG_2015, "PM1", "COSMIC_hotspot");
}

sub count_medalable_nt {
  my $seq_type = $FLAGS{"count-medalable-nt"};
  my %gene2source;
  if ($seq_type eq "somatic") {
    my $gold_genes = read_simple_file($FLAGS{"gold"} || die "-gold");
    my $silver_genes = read_simple_file($FLAGS{"silver"} || die "-silver");
    my $gl_reportable = read_simple_file(($FLAGS{"gl-reportable-genes-ger"} ||
					  die "-gl-reportable-genes-ger"));
    die "convert to use main reportable list + GSM, gl-reportable-genes-ger obsolete";

    my $f_manual = $FLAGS{"genes-manual"} || die;
    my $df = new DelimitedFile("-file" => $f_manual,
			       "-headers" => 1,
			      );
    my @genes_manual;
    while (my $row = $df->get_hash()) {
      push @genes_manual, $row->{Gene} || die;
    }

    foreach my $g (@{$gold_genes}) {
      $gene2source{$g} = "gold";
    }
    foreach my $g (@{$silver_genes}) {
      $gene2source{$g} = "silver";
    }
    foreach my $g (@genes_manual) {
      $gene2source{$g} = "manual";
    }
    foreach my $g (@{$gl_reportable}) {
      $gene2source{$g} = "germline_reportable";
    }

  } elsif ($seq_type eq "germline") {
    die "germline not implemented";
    # maybe separate options for REVIEWABLE and REPORTABLE
  } else {
    die "specify somatic/germine\n";
  }

  my $rff = new RefFlatFile();
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die;
  my $f_refflat = $config_genome->{REFSEQ_REFFLAT} || die;
  $rff->strip_sharp_annotations(1);
  $rff->parse_file(
		   "-refflat" => $f_refflat,
		   "-type" => "refgene",
		  );
  $rff->cache_base_translations(0);

  my %refflat_genes;
  foreach my $row (@{$rff->rows}) {
    $refflat_genes{$row->{gene}} = 1 if $row->{gene};
  }

  my $gsm = new_gsm();
  $gsm->hgnc_synonym_enable(0);
  foreach my $gene (keys %refflat_genes) {
    $gsm->add_gene(
		   "-gene" => $gene,
		  );
  }

  my $sjpi = get_sjpi();

  my $rpt = new Reporter(
			 "-file" => "bases_report.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   Gene_raw
					   Gene_refflat
					   preferred_isoform
					   used_isoform
					   mapping_count
					   count_coding
					   count_splice_region
					   count_total
					   note
					)
				      ],
			 "-auto_qc" => 1,
			);


  foreach my $gene_raw (sort keys %gene2source) {
    # - find the preferred isoform for this gene
    # - find refFlat entry/entries
    # - count bases
    my $nm_preferred = $sjpi->get_preferred_isoform($gene_raw);
    my $gene_adj = $gene_raw;
    unless ($nm_preferred) {
      $gene_adj = $gsm->resolve("-symbol" => $gene_raw);
      if ($gene_adj) {
	$nm_preferred = $sjpi->get_preferred_isoform($gene_adj) || "";
	# sometimes can't be resolved, e.g. TERC with is NR_
      }
    }

    my %r;
    $r{Gene_raw} = $gene_raw;
    $r{Gene_refflat} = $gene_adj;
    my $note = "";
    my $rows;

    if ($nm_preferred) {
      $r{preferred_isoform} = $nm_preferred;
      my $nm_final = $nm_preferred;

      unless ($rows = $rff->find_by_accession($nm_preferred)) {
	my $set = $sjpi->gene2nms()->{$gene_adj} || die;
	foreach my $nm (@{$set}) {
	  if ($rows = $rff->find_by_accession($nm)) {
	    $note = "no refFlat record for $nm_preferred, using $nm";
	    $nm_final = $nm;
	    last;
	  }
	}
      }
      unless ($rows) {
	# eg. TET3: preferred NM_144993, however this replaced by NM_001287491
	if (my $hits = $rff->find_by_gene($gene_adj)) {
	  die "multi rows" if @{$hits} > 1;
	  # deal with multiple isoforms, etc. if/when this happens
	  my $hit = $hits->[0];
	  my $nm = $hit->{name} || die;
	  $note = "no refFlat record for $nm_preferred, using $nm";
	  $nm_final = $nm;
	  $rows = [ $hit ];
	}
      }
      die "no refflat records for $gene_adj" unless $rows;

      $r{mapping_count} = scalar @{$rows};

    } else {
      $r{preferred_isoform} = "";
      $r{mapping_count} = "";
      $note = "no available SJ preferred";
    }

    if ($rows) {
      my $total_coding = 0;
      my $total_splice_region = 0;
      # across all isoform mappings

      foreach my $row_nm (@{$rows}) {
	my %coding_exons;
	my $summary = $rff->get_feature_summary("-row" => $row_nm);
	foreach my $s (@{$summary}) {
	  my $feature = $s->{feature} || die;
#	  dump_die($s, "debug", 1);
	  if ($feature eq "exon") {
	    $coding_exons{$s->{feature_number}} = 1;
	    my $start = $s->{start} || die;
	    my $end = $s->{end} || die;
	    die if $start > $end;
	    $total_coding += ($end - $start) + 1;
	  }
	}
	my $count_coding = scalar keys %coding_exons;
	dump_die($row_nm, "no coding exons?") if $count_coding < 1;
	my $count_splice_region = ($count_coding - 1) * 2 * 10;
	# - count of coding junctions,
	# - times 2 (donor and acceptor site),
	# - times 10 (10 nt at each edge, i.e. splice + splice_region calls)
	$total_splice_region += $count_splice_region;
      }
      $r{count_coding} = $total_coding;
      $r{count_splice_region} = $total_splice_region;
      $r{count_total} = $total_coding + $total_splice_region;
    } else {
      $r{count_coding} = "";
      $r{count_splice_region} = "";
      $r{count_total} = "";
    }

    $r{used_isoform} = $nm_preferred;
    $r{note} = $note;

    $rpt->end_row(\%r);

#    printf STDERR "%s\n", join " ", $gene_raw, $gene_adj, $nm, $gene_raw eq $gene_adj ? 0 : 1;
  }

  $rpt->finish();
}

sub hack_sjpi {
  my $sjpi = get_sjpi();
  die $sjpi->get_preferred_isoform("KMT2C") || "unknown";
}

sub run_svs_new {
  # 9/2018 update
  my $svc = $FLAGS{"sv"} || die "specify -sv";
  my $lines = read_simple_file($svc);
  # primary SVs

  my $use_single_config = 1;
  # use original gene lists rather than localized ones (eliminate config vars)
  my $use_gsm_index = 1;
  # use previous gene symbols and aliases to find genes
  my $use_ts_onco = 1;
  # use new database for tumor suppressor/oncogene annotations.
  # this may affect some medals, hopefully in a good way.
  # consistency problems may reveal problematic TS/onco annotations
  # which can be repaired.
  my $use_all_somatic_manual_genes = 1;
  # bug fix: also include the "manual" somatic gene list
  my $use_all_medal_matches = 1;
  # bug fix: when reporting genes related to medal, include all events at
  # same medal level (gold/silver) rather than just one

  my ($gene2class, $ts_onco_db);
  if ($use_ts_onco) {
    $ts_onco_db = new TSOncoDB(
			       "-gsm" => new_gsm_lite()
			      );
  } else {
    $gene2class = load_gene_class_files();
    # tumor suppressor/oncogene annotations
    # TO DO: remove gene2class altogether in future version
  }

  #
  #  load reference SVs:
  #

  # manually-curated SVs:
  my $svc_manual = $FLAGS{"sv-manual"} || die "-sv-manual";
  line_delimiter_qc("-file" => $svc_manual, "-delimiter" => "\t");
  printf STDERR "SV manual: %s\n", $svc_manual;
  my $df = new DelimitedFile(
			     "-file" => $svc_manual,
			     "-headers" => 1,
			    );
  while (my $row = $df->get_hash()) {
    my ($gene1, $gene2, $submitter) = @{$row}{qw(pair1 pair2 contact)};
    push @{$lines}, sprintf "%s\n", join "\t", $gene1, $gene2, $submitter;
    # reformat for compatibility w/main set
  }

  my @reference_svs;
  #
  # main reference SVs:
  #
  my $sv_id = 0;
  my %gene2svids;
  my %svid2sv;
  my %saw_load;
  my $skipped_duplicates = 0;
  foreach my $line (@{$lines}) {
    my @f = split /\t/, $line;
    die scalar @f unless @f == 3;
    my ($gene1, $gene2, $stuff) = @f;
    die unless $gene1 or $gene2;

    my $key = join ".", grep {$_ and /\w/} ($gene1, $gene2);
#    die "duplicate $key" if $saw_load{$key};
    if ($saw_load{$key}) {
      $skipped_duplicates++;
      next;
    }
    $saw_load{$key} = 1;
    # skip duplicate entries during loading process to aid ambiguity matching

    my %sv;
    $sv{breakpoints} = [];
    my $sv_id = ++$sv_id;
    $sv{sv_id} = $sv_id;
    my %genes_this;
    foreach ($gene1, $gene2) {
      if (defined($_) and /\w/) {
	push @{$sv{breakpoints}}, $_;
	$genes_this{$_} = 1;
	# each reference SV:
	# - may contain one or two gene breakpoints
	# - each breakpoint refers to exactly one gene
      }
    }
    $sv{breakpoint_count} = scalar @{$sv{breakpoints}};
    push @reference_svs, \%sv;

    # index SVs by ID number and native gene symbol:
    $svid2sv{$sv_id} = \%sv;
    foreach my $g (keys %genes_this) {
      $gene2svids{$g}{$sv_id} = 1;
    }
  }
  printf STDERR "skipped %d SV duplicates during loading\n", $skipped_duplicates if $skipped_duplicates;


  # for mapping user gene symbols to native SV symbols:
  my $gsm_sv_gene = new_gsm_lite();
  foreach (keys %gene2svids) {
    $gsm_sv_gene->add_gene("-gene" => $_);
  }

  #
  #  separate "gold" gene lists:
  #
  my $gsm_gold_genes = new_gsm_lite();

  my @load_flags;
  if ($use_single_config) {
    push @load_flags, "gold";
    if ($use_all_somatic_manual_genes) {
      my $mf = $FLAGS{"genes-manual"} || die "-genes-manual";
      # this is needed to include SMC3 and other manually-added genes
      my $df = new DelimitedFile("-file" => $mf,
				 "-headers" => 1,
				);
      while (my $row = $df->get_hash()) {
	my $g = $row->{Gene} || die;
	$gsm_gold_genes->add_gene("-gene" => $g);
      }
    } else {
      # previous refFlat-localized version contained only SMC3, and
      # not any newer additions
      printf STDERR "SMC compatibility hack\n";
      $gsm_gold_genes->add_gene("-gene" => "SMC3");
      # a gene from the MANUAL gene config
    }
  } else {
    push @load_flags, "gold-genes-mapped-to-fb";
    # OBSOLETE: somatic gold genes mapped to refFlat
  }
  if ($use_single_config) {
    push @load_flags, "gl-reportable-genes";
  } else {
    # OBSOLETE: germline-reportable genes mapped to refFlat
    push @load_flags, "gl-reportable-genes-refflat";
  }

  foreach my $flag (@load_flags) {
    # 3/2016: also include germline-reportable genes
    # TO DO: replace w/standard somatic version(s) to eliminate config entries
    my $list = read_simple_file($FLAGS{$flag} || die $flag);
    foreach (@{$list}) {
      $gsm_gold_genes->add_gene("-gene" => $_);
    }
  }

  #
  #  separate "silver" gene list:
  #
  $df = new DelimitedFile(
			  "-file" => ($FLAGS{"sv-silver-genes-fb"} || die "-sv-silver-genes-fb"),
			  "-headers" => 1,
			 );
  my $gsm_silver_genes = new_gsm_lite();
  while (my $row = $df->get_hash()) {
    my $gene = $row->{gene} || die;
    $gsm_silver_genes->add_gene("-gene" => $gene);
  }

  $gsm_sv_gene->validate_genes("-label" => "SV main");
  $gsm_gold_genes->validate_genes("-label" => "SV gold");
  $gsm_silver_genes->validate_genes("-label" => "SV silver");


  #
  #  classify:
  #
  foreach my $infile (@INPUT_SV_FILES) {
    printf STDERR "processing %s...\n", $infile;
    my $df = new DelimitedFile(
			       "-file" => $infile,
			       "-headers" => 1
			      );

    my @out_rows;
    while (my $row = $df->get_hash()) {
      my $is_truncating;
      my @breakpoints;
      # terrible name: actually "sides" of the event / fusion partner lists
      if (exists $row->{$FUSION_BUILDER_COLUMN}) {
	# FusionBuilder format:
	# - either one or two partner lists
	#   (1 indicates internal rearrangement?)
	# - each side of breakpoint may contain one or more genes
	my @sides = split /_/, $row->{$FUSION_BUILDER_COLUMN};
	my %saw;
	foreach my $thing (@sides) {
	  die unless $thing =~ /\w/;
	  next if $saw{$thing};
	  $saw{$thing} = 1;
	  # special case for internal rearrangements, e.g.
	  # PAX5_PAX5, see below in CREST example

	  my $gsm = new_gsm_lite();
	  foreach my $gene (split /,/, $thing) {
	    $gene = clean_sharp_gene_symbol($gene);
	    $gsm->add_gene("-gene" => $gene);
	  }
	  push @breakpoints, $gsm;
	}
	die "no FusionBuilder Usage column present" unless exists $row->{Usage};
	my $usage = $row->{Usage} || "";
	# may be blank in very rare cases
	# e.g. /rgs01/resgen/prod/tartan/runs/crest-post/fnEMxF4l/intmd/SJRB012408_X2_G1/SJRB012408_X2_G1.counts
	$is_truncating = $usage eq "TRUNCATING";
      } else {
	dump_die($row, "can't find column $FUSION_BUILDER_COLUMN");
      }

      #
      #  match this SV vs. reference set:
      #
      my @sv_hits;
      my @sv_search;
      if ($use_gsm_index) {
	# find reference SVs to check against using gene symbol
	# disambiguation lookup:
#	printf STDERR "GSM index reference SVs\n";
	my %sv_ids;
	foreach my $gsm_sv (@breakpoints) {
	  foreach my $gene (@{$gsm_sv->get_gene_list()}) {
	    if (my $gene_sv = $gsm_sv_gene->find($gene)) {
	      if (my $h = $gene2svids{$gene_sv}) {
		foreach (keys %{$h}) {
		  $sv_ids{$_} = 1;
		}
	      }
	    }
	  }
	}
	foreach (sort {$a <=> $b} keys %sv_ids) {
	  push @sv_search, $svid2sv{$_} || die;
	}
      } else {
	# search everything: much slower, esp. w/gene symbol disambiguation
	@sv_search = @reference_svs;
      }

      foreach my $sv_ref (@sv_search) {
	my @found_i;
	my @fuzzy_i;
#	my @match_genes;
	foreach my $end_gene (@{$sv_ref->{breakpoints}}) {
	  # one gene per reference breakpoint
	  for (my $i = 0; $i < @breakpoints; $i++) {
	    my $gsm = $breakpoints[$i];
	    if (my $hit = $gsm->find($end_gene)) {
	      # reference breakpoint gene found in local breakpoint
#	      printf STDERR "lookup %s: hit %s\n", $end_gene, $hit;
#	      push @{$found_i[$i]}, $end_gene;
	      push @{$found_i[$i]}, $hit;
	      # use the local gene symbol (i.e. refFlat) rather than
	      # the config gene symbol, which may be older
	      # e.g. config MLL3 is refFlat KMT2D
	      $fuzzy_i[$i] = 1 if $gsm->get_gene_count() > SV_MAX_BREAKPOINT_GENES_FOR_GOLD;
	    }
	  }
	}

	my @match_genes = map {@{$_}} grep {$_} @found_i;

	my $local_end_hit_count = scalar grep {$_} @found_i;
	# count of ENDS hit, rather than gene count
	my @reasons;
	my ($is_fuzzy) = grep {$_} @fuzzy_i;
	push @reasons, "large_gene_list" if $is_fuzzy;

	my %hit;
	$hit{sv} = $sv_ref;
	$hit{match_genes} = \@match_genes;
	$hit{reasons} = \@reasons;

	if (@breakpoints == $sv_ref->{breakpoint_count} and
	    # local SV has the same number of breakpoints as the reference SV
	    $local_end_hit_count == $sv_ref->{breakpoint_count}
	    # all genes in the reference SV found in local SV,
	    # and all local breakpoints were covered.

	    # note this case should NOT match the example below,
	    # because all the gene symbols are found in GeneA:

	    # reference SV:
	    # /nfs_exports/genomes/1/projects/ClinicalSeq/FusionGenes/SVCheck.txt contains this SV involving C11orf95 and MAML2:
	    # C11orf95        MAML2   EPD;COSMIC
	    #
	    # In the CREST file /nfs_exports/genomes/1/projects/.allType/ClinicalPilot/.allTumor/SJETV/NextGen/WholeGenome/StructuralVariation/CREST/SJETV092_D_G.predSV there is an entry where GeneA contains a long list of genes including both C11orf95 and MAML2, however GeneB is blank.
	   ) {
	  # - local SV has the same number of breakpoints as the reference SV
	  # - both genes in reference SV found

	  my $handled;
	  if ($sv_ref->{breakpoint_count} == 1 and $is_truncating) {
	    # JZ 5/29/2014:
	    # If any single-gene marked as Recurrent, if it is a
	    # truncation, mark them as silver. For those marked as
	    # TS+Recurr, treat them as non-Recur.
	    my $gene = $sv_ref->{breakpoints}->[0];

	    my ($is_ts, $is_recurrent);
	    if ($use_ts_onco) {
	      # works, but experimental; needs discussion.
	      # maybe exclude tsgene?
	      $is_ts = $ts_onco_db->is_lof("-gene" => $gene);
	      $is_recurrent = $ts_onco_db->is_gof("-gene" => $gene);
	    } else {
	      # old style
	      if (my $class = $gene2class->{$gene}) {
		($is_ts, $is_recurrent) = parse_gene_class($class);
	      }
	    }

	    if ($is_recurrent and not($is_ts)) {
	      # perfect match, but a truncation in an oncogene
	      # which is less interesting than one in a tumor suppressor
	      unshift @reasons, ("perfect_match", "oncogene_truncation");
	      $hit{medal} = CLASS_SILVER;
	      $hit{sort_score} = 6;
	      $handled = 1;
	    }


	  }

	  unless ($handled) {
	    unshift @reasons, "perfect_match";
	    if ($is_fuzzy) {
	      $hit{medal} = CLASS_SILVER;
	      $hit{sort_score} = 5;
	      # pretty good because all genes hit,
	      # but less compelling because there are a lot of genes
	    } else {
	      $hit{medal} = CLASS_GOLD;
	      $hit{sort_score} = 10;
	    }
	  }

	  push @sv_hits, \%hit;
	} elsif (@match_genes) {
	  # a single gene matches, others don't (putative or reference)
	  $hit{sv} = $sv_ref;
	  unshift @reasons, "imperfect_match";
	  $hit{medal} = CLASS_SILVER;
	  $hit{sort_score} = 4;
	  # minor match
	  push @sv_hits, \%hit;
	}
      }

      #
      #  check for matches to separate gold/silver gene lists:
      #
      my @search = map {@{$_->get_gene_list()}} @breakpoints;

      #
      #  gold gene list:
      #
      my @gold_hits;
      foreach my $g (@search) {
	if (my $hit = $gsm_gold_genes->find($g)) {
	  push @gold_hits, $g;
	  # use user symbol (refFlat) which may be newer
	}
      }
      if (@gold_hits) {
	push @sv_hits, {
			"medal" => CLASS_SILVER,
			"sort_score" => 3,
			"reasons" => [ "gold_gene" ],
			"match_genes" => unique_ordered_list(\@gold_hits)
		       };
	# create new record so as not to interact with %hits above,
	# which may have already been added
      }

      #
      #  silver gene list:
      #
      my @silver_hits;
      foreach my $g (@search) {
	my $hit = $gsm_silver_genes->find($g);
	push @silver_hits, $g if $hit;
      }

      if (@silver_hits) {
	push @sv_hits, {
			"medal" => CLASS_SILVER,
			"sort_score" => 2,
			"reasons" => [ "silver_gene" ],
			"match_genes" => unique_ordered_list(\@silver_hits)
		       };
	# create new record so as not to interact with %hits above
      }

      #
      #  set medal based on ranked evidence:
      #
      my $medal;
      my @gene_hits;
      my $reason = "";

      if (@sv_hits) {
	@sv_hits = sort {$b->{sort_score} <=> $a->{sort_score}} @sv_hits;
	my $best = $sv_hits[0];

	if ($use_all_medal_matches) {
	  # report all matching genes at this medal level.  Argubably this
	  # is an improvement/bugfix vs. the original behavior (below).
	  # While this issue doesn't affect medal calls, it does affect
	  # the genes reported in the Evidence field.
	  my $best_medal = $best->{medal} || die;
	  my @sv_matches = grep {$_->{medal} eq $best_medal} @sv_hits;
	  die unless @sv_matches;

	  my $list = unique_ordered_list([ map {@{$_->{match_genes}}} @sv_matches ]);
	  @gene_hits = @{$list};
	} else {
	  printf STDERR "DEBUG: original mode\n";
	  # Just report genes matching the highest-scoring hit
	  # (original behavior).  The trouble here is that sometimes
	  # there may be multiple entries with different genes at the
	  # same medal level, e.g.
	  #
	  # ETV6,FOXP1      TRUNCATING      Gold    perfect_match   FOXP1
	  #
	  # Here only FOXP1 is reported, however ETV6 is also an
	  # equally legitimate match (both genes have a single-sided
	  # SV config entry).  Part of the trouble is these special
	  # comma-delimited single edged SV entries.  As I recall
	  # Michael R didn't like this reporting format as it's
	  # essentially two records reported in one row, but he was
	  # overruled.
	  #
	  # another example, this time for a traditional 2-gene fusion:
	  #
	  # < RELN_NUP98    GENIC   Silver  imperfect_match NUP98
	  # > RELN_NUP98    GENIC   Silver  imperfect_match NUP98,RELN
	  #
	  # NUP98 is not a perfect match to SV config entries (silver),
	  # and RELN is a germline-reportable gene (silver).  Both gene
	  # matches carry equal medaling weight, so why not mention both?
	  @gene_hits = @{$best->{match_genes}};
	}

	if ($ENV{SV_DEBUG}) {
	  printf STDERR "%s %s: dump of %d SV hits:",
	    $row->{$FUSION_BUILDER_COLUMN}, $row->{Usage},
	      scalar(@sv_hits);
	  foreach my $ref (@sv_hits) {
	    my $ref_info;
	    if ($ref->{sv}) {
	      $ref_info = sprintf "SV=%s", join ",", @{$ref->{sv}{breakpoints}};
	    } else {
	      $ref_info = "gold_or_silver_gene";
	    }
	    printf STDERR " %s (hits:%s ref:%s)",
	      $ref->{medal},
		join(",", @{$ref->{match_genes}}),
		  $ref_info;
	  }
	  print STDERR "\n";
	}

	$medal = $best->{medal};
	$row->{sort_score} = $best->{sort_score};
	$reason = join ",", @{$best->{reasons}};
      } else {
	$medal = "Unknown";
	$row->{sort_score} = 0;
      }

      $row->{GSBClass} = $medal;
      $row->{Reason} = $reason;
      $row->{Evidence} = join ",", @gene_hits;
      push @out_rows, $row;
    }

    my $rpt = $df->get_reporter(
				"-file" => get_outfile($infile),
				"-extra" => [
					     qw(
						 GSBClass
						 Reason
						 Evidence
					      )
					    ]
			       );
    $rpt->write_headers();
    # force header line even if empty file
    foreach my $row (sort {$b->{sort_score} <=> $a->{sort_score}} @out_rows) {
      $rpt->end_row($row);
    }
    $rpt->finish();
  }
  printf STDERR "SVs done\n";
}

sub run_cnvs_integrated_new {
  #
  # CNV classification for either somatic or germline.
  #
  my (%options) = @_;
  my $SOMATIC = 1;
  my $GERMLINE = 2;
  my $mode;

  my @cnv_configs;
  my $min_abs_copy_number;

  my $extend_cosmic_match = 1;
  # bug fix: previously, appending extended COSMIC annotation required
  # the main annotation field to be set to "COSMIC".  However there
  # are some cases where the main annotation includes multiple tags,
  # including COSMIC.
  my $use_gsm_annotations = 1;
  # use GSM for annotations rather than relying on hash match
  my $generate_germline_config = 1;
  # generate germline match details on the fly, eliminating a
  # genome config file (-gl-cnv CLINCLS_GL_CNV)
  my $use_original_reviewable_genes = 1;
  # use the native/raw reviewable gene list rather than the localized
  # version.  Hopefully we can eventually eliminate
  # -gl-reviewable-genes-ger / CLINCLS_GERMLINE_REVIEWABLE_GENES_GER

  if ($options{"-somatic"}) {
    printf STDERR "starting SOMATIC CNV classification...\n";
    $mode = $SOMATIC;
    push @cnv_configs, $FLAGS{"cnv"} || die "specify -cnv";
  } elsif ($options{"-germline"}) {
    printf STDERR "starting GERMLINE CNV classification...\n";
    $mode = $GERMLINE;
  } else {
    die "must specify either -somatic or -germline";
  }

  # 3/2016: somatic CNV should also report on germline reportable genes.
  # since germline CNV config only contains reportable genes, we
  # can just use that config for both germline and somatic runs.
  if ($generate_germline_config) {
    # generate on the fly
    push @cnv_configs, generate_germline_cnv_config("-ram" => 1);
  } else {
    # obsolete me
    push @cnv_configs, $FLAGS{"gl-cnv"} || die "-gl-cnv";
  }
  die unless @cnv_configs;

  my %annotations;

  my $gsm_all_genes = new_gsm_lite();
  my $gsm_amp_genes = new_gsm_lite();
  my $gsm_del_genes = new_gsm_lite();
  my $gsm_sv_genes = new_gsm_lite();
  my $gsm_annot = new_gsm_lite();

  #
  #  load target amplifications/deletions:
  #
  foreach my $cnv_check (@cnv_configs) {
    my $lines;
    if (ref $cnv_check) {
      # array of lines rather than a file
      $lines = $cnv_check;
    } else {
      $lines = read_simple_file($cnv_check);
    }
    foreach my $line (@{$lines}) {
      my ($gene, $type, $annot) = split /\t/, $line;
      $gsm_all_genes->add_gene("-gene" => $gene);
      if ($type eq "Amp") {
	$gsm_amp_genes->add_gene("-gene" => $gene);
      } elsif ($type eq "Del") {
	$gsm_del_genes->add_gene("-gene" => $gene);
      } else {
	die;
      }

      if ($annot) {
	die "duplicate annotation for $gene, split by type??" if $annotations{$gene} and $annotations{$gene} ne $annot;
	$annotations{$gene} = $annot;
      }
    }
  }

  if ($mode == $GERMLINE) {
    my $gl_reviewable;
    if ($use_original_reviewable_genes) {
      $gl_reviewable = get_germline_reviewable_genes("-original" => 1);
    } else {
      $gl_reviewable = get_germline_reviewable_genes("-ger" => 1);
    }

    # 3/2016: use mapping directly from Dale's master list which contains
    # both the HUGO and an alternate symbol.  This cuts out the middleman
    # vs. the germline intervals list originally created by Gang,
    # just in case there are issues with older reviewable symbols.
    foreach my $gene (@{$gl_reviewable}) {
      $gsm_all_genes->add_gene("-gene" => $gene);
      # if found in a CNV, these genes will get a silver
      # (we don't have tumor suppressor/oncogene annotations so
      # amplification/deletion related checks won't be performed)
    }
    # TO DO: update to use TSOncoDB???
  } elsif ($mode == $SOMATIC) {
    # load supplemental annotations:
    my $cnv_annot = $FLAGS{"cnv-annotations"} || die "specify -cnv-annotations";
    my $df = new DelimitedFile("-file" => $cnv_annot,
			       "-headers" => 1,
			      );
    # while (my $row = $df->next("-ref" => 1)) {  # headerless
    my %supplemental;
    while (my $row = $df->get_hash()) {
      my $gene = $row->{"Gene Symbol"} || die;
      $supplemental{$gene} = $row->{"Tumour Types(Germline)"};
    }

    foreach my $g (keys %annotations) {
      # supplement COSMIC annotations
      my $annot = $annotations{$g} || "";
      my $match = $extend_cosmic_match ? $annot =~ /COSMIC/ : $annot eq "COSMIC";
      if ($match and my $extra = $supplemental{$g}) {
	$annotations{$g} .= ";" . $extra;
      }
    }
  } else {
    die;
  }

  foreach my $gene (keys %annotations) {
    $gsm_annot->add_gene("-gene" => $gene);
  }

  my $annot = $FLAGS{"gene-exon-region-dir"} || die "specify -gene-exon-region-dir";
  # need genomically-sorted genes in interval for analysis
  printf STDERR "loading %s...", $annot;
  my $ga = new GeneAnnotation(
			      "-style" => "gene_exon_region",
			      "-gene_exon_region_dir" => $annot,
			      "-ignore_non_coding" => 0
			     );
  print STDERR "done\n";

  $gsm_all_genes->validate_genes("-label" => "CNV main");

  if ($mode == $SOMATIC) {
    # somatic only: not applicable to germline since we don't perform
    # germline SV analysis
    my $sv_genes = get_sv_genes();
    foreach my $g (keys %{$sv_genes}) {
      $gsm_sv_genes->add_gene("-gene" => $g);
    }
    $gsm_sv_genes->validate_genes("-label" => "CNV SV");
  }

  foreach my $infile (@INPUT_CNV_FILES) {
    printf STDERR "processing %s...\n", $infile;

    my $df = new DelimitedFile(
			       "-file" => $infile,
			       "-headers" => 1,
			      );
    my @out_rows;
    while (my $row = $df->get_hash()) {
      #      my $genes = $row->{Genes};
      #      foreach (sort keys %{$row}) {
      #	printf "%s: %s\n", $_, $row->{$_};
      #      }
      my $log;
      if ($mode == $SOMATIC) {
	$log = $row->{LogRatio};
	die "ERROR: LogRatio field missing in input file" unless defined $log;
	# LogRatio only appears in paired output (somatic)
      } elsif ($mode == $GERMLINE) {
	$log = $row->{"seg.mean"};
	# unpaired (germline) runs, see
	# http://hc-wiki.stjude.org/display/compbio/Clinical+Genomics+CNV+Manual+Analysis+SOP
	# note this is NOT on log scale so shouldn't be used the same
	# way as in somatic classification.  In germline classification
	# we just use it to detect amplification/deletion status,
	# for high-level amplification the separate "seg.mean" check
	# still applies.
	die "ERROR: seg.mean field missing in input file" unless defined $log;
      } else {
	die;
      }

      my $medal;
      $row->{sort_score} = 0;

      my $chrom = $row->{chrom} || die;
      my $start = $row->{"loc.start"} || dump_die($row, "no loc.start");
      my $end = $row->{"loc.end"} || dump_die($row, "no loc.end");

      if ($start > $end) {
	# some records seem to have start and end swapped,
	# /cgs01/clingen/prod/tartan/runs/classify-cnv-somatic/SCKIrJDS/input/SJOS030272_D1_G1/WHOLE_GENOME/conserting-crest/SJOS030272_D1_G1_CONSERTING_Mapability_100.txt.QualityMerge
	# /home/medmonso/work/dale/2016_12_19_somatic_CNV_crash
	printf STDERR "WARNING: start > end for %s.%s.%s, swapping\n",
	  $chrom, $start, $end;
	my $temp = $start;
	$start = $end;
	$end = $temp;
      }

      $ga->find(
		"-reference" => $chrom,
		"-start" => $start,
		"-end" => $end
	       );
      #      my @genes_report = split /,/, $genes;
      my $genes_genomic = $ga->results_genes_genomic_order();
      # genes in interval, sorted genomically
      my $gene_count = scalar @{$genes_genomic};
      #      printf STDERR "genes: %d\n", $gene_count;
      my $is_intronic = $ga->is_intronic();

      # 10/16/2013:
      # add gene annotation columns, replacing now-retired
      # earlier annotation step
      $row->{"Number ofGenes"} = $gene_count;
      # (sic)
      $row->{Genes} = $gene_count > $CNV_MAX_GENES_ANNOTATE ? sprintf("> %d genes", $CNV_MAX_GENES_ANNOTATE) : join(",", @{$genes_genomic});

      my @medal_genes;
      my @reasons;
      my $is_focal = $gene_count <= CNV_MAX_GENES_FOCAL;

      my $is_deletion;
      my $gsm_search;
      if ($log < 0) {
	$is_deletion = 1;
	$gsm_search = $gsm_del_genes;
      } else {
	$gsm_search = $gsm_amp_genes;
      }
      die "ERROR, missing seg.mean" unless defined $row->{"seg.mean"};
      my $copy_number = $row->{"seg.mean"} * 2;

      if (CNV_HIGH_LEVEL_ENABLE) {
	# if a CNV is highly amplified or deleted
	if ($copy_number >= CNV_MIN_COPY_NUMBER_GAIN_FOR_GOLD or
	    $copy_number <= CNV_HIGH_LEVEL_DELETION_COPIES) {
	  # CNV is a high-level amplification or deletion
	  # search ALL contained genes vs. the watch list

#	  printf "high-level event: %s, genes=%s\n", $copy_number, join ",", @{$genes_genomic};

	  foreach my $gene (@{$genes_genomic}) {
	    if ($gsm_search->find($gene)) {
	      push @medal_genes, $gene;
	    }
	  }
	  if (@medal_genes) {
	    push @reasons, general_label($is_focal, $is_deletion);
	    # downstream users might want/expect this basic info,
	    # so leave it, however it seems a little misleading
	    # because these details are unrelated to the logic behind the
	    # event.
	    push @reasons, sprintf "high_level_%s", $is_deletion ? "deletion" : "amplification";
	    $medal = CLASS_GOLD;
	    $row->{Genes} = join(",", @{$genes_genomic});
	    # if medaled, annotate genes regardless of CNV size
	    # (remove usual cap)
	  }
	}
      }

      my $mode_germline_any = ($mode == $GERMLINE and $CNV_GERMLINE_MEDAL_ANY_SIZE);

      if ($medal) {
	# already handled
      } elsif ($mode_germline_any ? 0 : $gene_count > $CNV_MAX_GENES_MEDAL) {
	$medal = CLASS_UNKNOWN;
      } else {
	#
	#  focal event:
	#
	if ($is_focal) {
	  foreach my $g (@{$genes_genomic}) {
	    push @medal_genes, $g if $gsm_search->find($g);
	  }
	  if (@medal_genes) {
	    push @reasons, general_label($is_focal, $is_deletion);
	    if ($is_deletion) {
	      $medal = CLASS_GOLD;
	    } else {
	      # amplification
	      if ($mode == $SOMATIC) {
		if ($log >= CNV_MIN_LOG2_FOR_HLA) {
		  $medal = CLASS_GOLD;
		  push @reasons, "high_level_amplification";
		} else {
		  $medal = CLASS_SILVER;
		  #		push @reasons, "low_level_amplification";
		  # remove for simplicity: only use "high_level_amplification"
		}
	      } elsif ($mode == $GERMLINE) {
		$medal = CLASS_GOLD;
		# JZ 5/23/2014:
		# I just had a brief discussion with Zhaojie. I think
		# that his concern is valid. We do need to implement CNV
		# for germline somewhat differently.  One big difference
		# is that for any CNVs that affect the coding exons of
		# the 31 genes currently used for reporting, we should
		# mark them as Gold. We can also be a bit more
		# sophisticated to make the tumor suppressor genes
		# “Gold” with deletion and oncogene “Gold” with
		# amplification.  We actually never implemented Germline
		# CNV or SV medal classification. It will be good if we
		# can throw something quick and simple together.
		#
		# JZ 5/27/2014:
		# For germline CNV, we do not need to distinguish
		# high-amplification or low amplification. Just treat
		# every germline amplification as if they were somatic
		# high-amplification.
		#
		# MNE notes:
		# - the genes in germline CNV config file currently
		#   perfectly match the the reportable gene list.
		#   So, this logic will assign gold if the amp/del status
		#   matches the config file (i.e. oncogene/tumor suppressor),
		#   and silver later if only the gene symbol matches.
	      } else {
		die;
	      }
	    }
	  }
	}

	if (not($medal) and
	    ($gene_count == 1 or $gene_count > CNV_MAX_GENES_FOCAL)) {
	  # - broad events
	  # - single-gene focal events that did not already receive a medal
	  my @search_genes;
	  if ($gene_count == 1) {
	    @search_genes = @{$genes_genomic};
	  } elsif ($mode_germline_any) {
	    @search_genes = @{$genes_genomic};
	  } else {
	    @search_genes = ($genes_genomic->[0],
			     @{$genes_genomic}[$#$genes_genomic]);
	  }

	  # first check: look for perfect matches,
	  # (will only ever happen for the non-focal search type)
	  foreach my $g (@search_genes) {
	    push @medal_genes, $g if $gsm_search->find($g);
	  }
	  if (@medal_genes) {
	    die if $is_focal;
	    # "that's unpossible!"

	    if ($is_deletion) {
	      $medal = CLASS_GOLD;
	    } else {
	      # amplification
	      if ($mode == $SOMATIC) {
		if ($log >= CNV_MIN_LOG2_FOR_HLA) {
		  push @reasons, "high_level_amplification";
		  $medal = CLASS_GOLD;
		  # JZ 9/19/2013:
		  # For amplification, I am thinking about marking
		  # amplification >4x as Gold. This will represent those
		  # that have log2ratio >=2. I looked over the list of
		  # amplification and felt that most of them are well
		  # characterized oncogenes expected to have very
		  # high-level of amplification. We can mark anything
		  # with log2ration <=2 as silver.
		  #
		  # 10/9/2013: in light of the below should this be silver??
		} else {
		  # JZ 10/9/2013:
		  # I thought that we can assign bronze to the low-level
		  # amplification (>5 genes)?
		  $medal = CLASS_BRONZE;
		}
	      } elsif ($mode == $GERMLINE) {
		$medal = CLASS_GOLD;
		# see comments in focal code above
	      } else {
		die;
	      }
	    }
	  } else {
	    # second check: look for matches to all genes regardless
	    # of amplification/deletion status
	    foreach my $g (@search_genes) {
	      push @medal_genes, $g if $gsm_all_genes->find($g);
	    }
	    if (@medal_genes) {
	      $medal = CLASS_SILVER;
	      push @reasons, "global_gene_match";
	      # match to global gene list regardless of amp/del status
	    }
	  }
	  if ($medal) {
	    if ($medal and @search_genes > 1) {
	      my $is_breakpoint;
	      foreach my $g (@medal_genes) {
		$is_breakpoint = 1 if $genes_genomic->[0] eq $g or
		  $genes_genomic->[$#$genes_genomic] eq $g;
	      }
	      unshift @reasons, "breakpoint_gene" if $is_breakpoint;
	    }
	    unshift @reasons, general_label($is_focal, $is_deletion);
	  }
	}			# broad

	if (not($medal)) {
	  #
	  #  silver for any gene in SV list:
	  #
	  my @search;
	  if ($is_focal or $gene_count == 1) {
	    @search = @{$genes_genomic};
	  } elsif ($mode_germline_any) {
	    @search = @{$genes_genomic};
	  } else {
	    @search = ($genes_genomic->[0],
		       @{$genes_genomic}[$#$genes_genomic]);
	  }

	  my $hit;
	  foreach my $g (@search) {
	    if ($gsm_sv_genes->find($g)) {
	      push @medal_genes, $g;
	      $hit = 1;
	    }
	  }

	  if ($hit) {
	    $medal = CLASS_SILVER;
	    push @reasons, "sv_gene_match";
	    my $is_breakpoint;
	    foreach my $g (@medal_genes) {
	      $is_breakpoint = 1 if $genes_genomic->[0] eq $g or
		$genes_genomic->[$#$genes_genomic] eq $g;
	      # TO DO: GSM fix??
	    }
	    unshift @reasons, "breakpoint_gene" if $is_breakpoint;
	    unshift @reasons, general_label($is_focal, $is_deletion);
	  }
	}
      }

      $row->{copy_number_gain} = $copy_number;
      if ($copy_number > 0 and $medal and $medal ne CLASS_UNKNOWN) {
	#
	#  supplemental/override logic for amplifications, per JZ 10/8/2013.
	#  only applies if a medal has been assigned (i.e. some kind of match)
	#
	if ($copy_number >= CNV_MIN_COPY_NUMBER_GAIN_FOR_GOLD) {
	  # a) For copy number gain >10, report amplification as gold regardless of num. of genes included in the interval.
	  # [presumably >= 10]
	  $medal = CLASS_GOLD;
	  #	  push @reasons, "high_level_amplification";
	  add_reason(\@reasons, "high_level_amplification");
	  # hopefully not a duplicate
	} elsif ($is_focal) {
	  $medal = CLASS_GOLD;
	  # b)For copy number gain <10, report gold for focal (<=5 genes)
	} else {
	  $medal = CLASS_SILVER unless $medal eq CLASS_GOLD;
	  # - silver for large segment (>5 genes)
	  # - don't override a better medal if assigned earlier
	}
      }

      push @reasons, sprintf "genes=%s", join "|", @medal_genes if @medal_genes;

      $medal = CLASS_UNKNOWN unless $medal;
      $row->{intronic_only} = $is_intronic;

      if ($mode == $GERMLINE and $is_intronic and $medal ne CLASS_UNKNOWN) {
	# would have received a medal, but intronic
	$medal = CLASS_UNKNOWN;
	unshift @reasons, "intronic";
      }

      if ($mode == $GERMLINE and
	  $medal ne CLASS_UNKNOWN and
	  defined $CNV_GERMLINE_MIN_COPY_DELTA) {
	unless (abs($copy_number) >= $CNV_GERMLINE_MIN_COPY_DELTA) {
	  push @reasons, "copy_number_below_threshold";
	  $medal = CLASS_UNKNOWN;
	}
      }

      $row->{GSBClass} = $medal;
      $row->{Reason} = join ",", @reasons;
      $row->{sort_score} = $MEDAL_SORT_WEIGHT{$medal} || die "no weight for $medal";

      my @annotations;
      foreach my $gene (@medal_genes) {
	if ($use_gsm_annotations) {
	  if (my $sym = $gsm_annot->find($gene)) {
	    printf STDERR "huzzah: found $sym by way of $gene\n"
	      if $sym ne $gene;
	    push @annotations, $annotations{$sym} || die "no annotation for $sym";
	  }
	} else {
	  if (my $thing = $annotations{$gene}) {
	    # won't be present if medal is based on a SV gene only
	    push @annotations, $thing;
	  }
	}
      }
      $row->{Evidence} = join "|", @annotations;

      push @out_rows, $row;
    }

    my @sorted = sort {$b->{sort_score} <=> $a->{sort_score}} @out_rows;
    my $rpt = $df->get_reporter(
				"-file" => get_outfile($infile),
				"-extra" => [
					     (
					      "Number ofGenes",
					      "Genes",
					      # 10/16/2013:
					      # we are now doing the
					      # gene annotation
					      "copy_number_gain",
					      "intronic_only",
					      "GSBClass",
					      "Reason",
					      "Evidence",
					     )
					    ]
			       );
    $rpt->write_headers();
    # force header line even if empty
    foreach my $row (@sorted) {
      $rpt->end_row($row);
    }
    $rpt->finish();
  }
  printf STDERR "CNVs done\n";
}

sub get_germline_reviewable_genes {
  # germline reviewable genes
  my (%options) = @_;
  my $flag;
  if ($options{"-original"}) {
    # the original/raw list, not guaranteed to be localized
    $flag = "gl-gold-cancer-ranges";
  } elsif ($options{"-ger"}) {
    # localized to GENE_EXON_REGION; hopefully OBSOLETE THIS
    $flag = "gl-gold-cancer-ranges-ger";
  } else {
    confess "specify -original or -ger";
  }
  my $f_reviewable = $FLAGS{$flag} || die "$flag";

  my $fh = new FileHandle();
  $fh->open($f_reviewable) || die;
  my @genes;
  while (<$fh>) {
    chomp;
    my @f = split /\t/, $_;
    die unless @f == 2;
    push @genes, $f[0];
  }
  return \@genes;
}

sub get_vm_mc_db {
  my %options = @_;
  my $f_db = $FLAGS{"sqlite"} || die "-sqlite";
  my $dbi = DBI->connect("dbi:SQLite:dbname=$f_db","","");
  my $svdb = new SimpleVariantDB(
				 "-dbi" => $dbi,
				);
  my $rows = $svdb->get_rows_with_constants(%options);
  die "no data for " . $options{"-name"} unless @{$rows};
  my $vm = new VariantMatcher();
  my $idx_aa = 0;
  my $idx_genomic = 0;
  my $row_count++;
  foreach my $row (@{$rows}) {
    $row_count++;
    my ($chr, $pos, $ra, $va, $gene, $aa) =
      @{$row}{F_VDB_CHR,
		F_VDB_POS,
		  F_VDB_RA,
		    F_VDB_VA,
		      F_VDB_GENE,
			F_VDB_AA
		      };

    if ($gene and $aa) {
#      printf STDERR "add $gene $aa\n";
      $vm->add_aa("-gene" => $gene, "-aa" => $aa, "-row" => $row);
      $idx_aa++;
    }

    if ($chr and $pos and ($ra or $va)) {
      my $v = new Variant();
      $v->import_generic(
			 "-reference-name" => $chr,
			 "-base-number" => $pos,
			 "-reference-allele" => $ra,
			 "-variant-allele" => $va
			);
      $vm->add_variant($v, "-row" => $row);
      $idx_genomic++;
    }
  }
  printf STDERR "db for %s: %d rows, indexed %d genomic, %d AA\n", $options{"-name"}, $row_count, $idx_genomic, $idx_aa;

  return $vm;
}

sub get_vm_user_variants {
  # DEVELOPMENT
  # - to add?:
  #   - gene/aa lookup (AA can of worms)
  #   - user database name
  #   - GSB call for different match categories
  #   - variant-level GSB call?
  my %options = @_;
  my $param = $options{"-param"} || die "-param";
  die if $param =~ /^\-/;
  my $f_db = $FLAGS{$param};
  my $vm = new VariantMatcher();
  if ($f_db) {
    # optional
    my $df = new DelimitedFile("-file" => $f_db,
			       "-headers" => 1,
			      );
    my $idx_genomic = 0;
    my $row_count = 0;
    my $idx_aa = 0;
    while (my $row = $df->get_hash()) {
      my $v = new Variant();
      $v->import_bambino_row("-row" => $row);
      # genomic details in raw Bambino format, e.g. vcf2tab.pl output
      # if we ultimately go with VCF
      $row_count++;

      $vm->add_variant($v, "-row" => $row);
      $idx_genomic++;
    }
    printf STDERR "db for %s: %d rows, indexed %d genomic, %d AA\n", "user_variants", $row_count, $idx_genomic, $idx_aa;
  }

  return $vm;
}

sub hack_user_variants {
  foreach my $f_gl (@INPUT_GL_FILES) {
    my $df = new DelimitedFile("-file" => $f_gl,
			       "-headers" => 1,
			     );
    my $vm_user_variants = get_vm_user_variants();
    while (my $row = $df->get_hash()) {
      my $is_silent_intron_utr = 0;
      my @reasons;
      my $is_novel_these_db;
      my $medal;
      $medal = assign_aa_match(
				 "-vm" => $vm_user_variants,
				 "-row" => $row,
				 "-is-silent" => $is_silent_intron_utr,
				 "-silent-allowed" => 0,
				 "-protect" => 0,
#				 "-transcript-handshake" => $FIELD_REFSEQ,
				 "-current-medal" => $medal,
				 "-reasons" => \@reasons,
				 "-label" => "user_variant",
				 "-novel-flag" => \$is_novel_these_db,
				 "-medal-perfect" => CLASS_GOLD,
				 "-medal-codon" => CLASS_SILVER,
				 "-medal-site" => CLASS_SILVER,
				 "-genomic-lookup" => 1,
#				 "-pubmed-list" => \@pmid,
#				 "-pubmed-field" => F_VDB_PMID()
				);
      die @reasons;
    }

  }
}
