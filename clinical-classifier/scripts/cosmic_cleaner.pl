#!/bin/env perl
# clean up issues with COSMIC database
# see /nfs_exports/genomes/1/projects/ClinicalSeq/README
# MNE 8/2013
#
#  TO DO:
#  - better parsing of AA change annotations and handling of "?" types
#  - fix inconsistent gene annotations FOR THE SAME VARIANT
#    which are not caught by auto-annotation:
#    e.g. multiple genes for 7:100552435-100552435_c.1186T>A: MUC3A,MUC3B
#  - dynamic NHLBI filtering rather than blacklist file

use strict;
use warnings;

use Getopt::Long;
use Carp qw(confess);
use POSIX qw(ceil);
use File::Basename;

use DelimitedFile;
use GeneAnnotation;
use FileUtils qw(read_simple_file universal_open find_binary);
use Reporter;
use GenomeUtils qw(reverse_complement cook_chromosome_name);
use GeneListCollapser;
use VariantMatcher;
use MiscUtils qw(dump_die);
use List::Util qw(sum);
use DBTools qw(export_query_to_flatfile get_dbi_gedi);
use TdtConfig;
use NucleotideSubstitutionParser;
use TemporaryFileWrangler;
use TabixFile;
use TabixPrep;
use TARTANUtils qw(tartan_genome_put_helper);
use DelimitedFileHP;

use constant INDEL_NUCLEOTIDE_NT_FUZZY_MATCH => 3;

my $EXCLUDE_HYPERMUTABLE = 1;
my $REANNOTATE_GENE_SYMBOLS = 1;
my $TRIM_ENSEMBL_ID_IF_ACCESSION_MISMATCHES = 1;
# Gene name: KRBOX1_ENST00000418176
# Accession Number: ENST00000426937
# - transcript suffix doesn't match, however KRBOX1 is a valid symbol
#my $COSMIC_OUTFILE = "cosmic_fixed.tab";
my $COSMIC_OUTFILE;

my $TOLERATE_NONSENSE = 0;

my $RECURRENCE_SAMPLE_COUNT = 10;
# how many samples a variant should appear in before being considered recurrent

my $GENE_SUMMARY_MAX_AA_COUNT = 10;
# number of AAs to report in gene summary report

#my $IGNORE_SILENT = 0;
#my $IGNORE_SILENT = 1;
my $IGNORE_SILENT = 0;
# 4/2015: back to keeping them in, sometimes these are misclassified
# and actually splice, e.g. TP53 T125T

my $HYPER_MIN_CODING_TO_SEED = 100;
# minimum coding mutations to trigger hypermutable analysis
my $HYPER_TOP_FRACTION_TO_INCLUDE = 0.10;
# examine this fraction of tumors with most mutations in a study
my $HYPER_MAX_TUMOR_FRACTION_HYPER = 0.10;
# warn if fraction of tumors called hypermutable vs. total exceeds
# this threshold (may need to revisit top X fraction)

my $F_TABIX_CHROM = "Chr";
my $F_TABIX_START = "Start";
my $F_TABIX_END = "End";

#my $F_C_HISTOLOGY_SUBTYPE = "Histology subtype";
# v72
my $F_C_HISTOLOGY_SUBTYPE = "Histology subtype 1";
# v76

my %FLAGS;
GetOptions(\%FLAGS,
	   "-cosmic=s",
	   # COSMIC FTP file (compressed OK)

#	   "-gene-exon-region-dir=s",
	   # gene annotation files, e.g.
	   # /nfs_exports/apps/gnu-apps/NextGen/nextgensupport/hg19/GENE_EXON_REGION/

	   "-hypermutable-samples=s",
	   # list of hypermutator sample IDs, e.g.
	   # /nfs_exports/genomes/1/projects/ClinicalSeq/hypermutable_sample_0819_2013.lst

	   "-cancer-gene-list=s",
	   # cancer gene list, e.g.
	   # /nfs_exports/genomes/1/projects/ClinicalSeq/cancer_gene.lst

	   "-bad-gene-list=s",
	   # bad gene list, e.g.
	   # /nfs_exports/genomes/1/projects/ClinicalSeq/bad_gene.lst

	   "-gene2nm=s",
	   # mapping of gene symbol to NM_, e.g.
	   # /nfs_exports/apps/gnu-apps/NextGen/nextgensupport2/NCBI/RefSeq/gene_transcript_matrix.withNM.mod

	   "-max-lines=i",
	   # debug

	   "-refgene-fasta=s",
	   # refGene FASTA files (used to get length), e.g. human.rna.fna.gz
	   # ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.rna.fna.gz

	   "-outfile=s" => \$COSMIC_OUTFILE,

	   "-reannotate=i" => \$REANNOTATE_GENE_SYMBOLS,
	   "-no-summary",
	   "-no-recurrence",

	   "-tolerant=i" => \$TOLERATE_NONSENSE,
	   # specify to prevent crash complaints about multiple gene symbols
	   # should be revisited

	   "-no-nhlbi",

	   "-ignore-silent" => \$IGNORE_SILENT,

	   "-nhlbi-blacklist-mk2=s",

	   "-bad-literature=s",
	   "-bad-variants=s",

	   "-somatic-site-rescue",
	   "-final=s",

	   "-find-hypermutable-samples",
	   # obsolete, xiaotu now doing

	   "-disease-hack",
	   "-dump-gedi-diagnosis",

	   "-hyper-prep=s",
	   # prep for Xiaotu's hypermutator analysis
	   # file is list of hypermutable TUMOR IDs from
	   # -find-hypermutable-samples

	   "-genome=s",
	   "-hypermutable-tumor-config=s",

	   "-reject-summary=s",

	   "-dump-gedi-validated",
	   "-gedi-validated=s",

	   "-hack",

	   "-clean2tabix=s",
	   "-tartan-index=s",
	   # output dir

	   "-f-histology-subtype=s" => \$F_C_HISTOLOGY_SUBTYPE,

	   "-append-occurrences=s",
	  );

if ($FLAGS{"reject-summary"}) {
  reject_summary($FLAGS{"reject-summary"});
  exit(0);
} elsif ($FLAGS{"dump-gedi-validated"}) {
  dump_gedi_validated();
  exit(0);
} elsif (my $if = $FLAGS{clean2tabix}) {
  clean2tabix($if);
  exit(0);
} elsif (my $if2 = $FLAGS{"append-occurrences"}) {
  append_occurrences($if2);
  exit(0);
} elsif ($FLAGS{hack}) {
  hack_test();
  exit(0);
} elsif (my $out_dir = $FLAGS{"tartan-index"}) {
  tartan_genome_put_helper(
			   "-genome" => $FLAGS{genome},
			   "-subdir" => "CLINICAL/COSMIC/COSMIC_tabix",
			   "-out" => $out_dir,
			  );
  exit(0);
}

my $genome = $FLAGS{genome} || die "-genome";
my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";

unless ($COSMIC_OUTFILE) {
  if ($IGNORE_SILENT) {
    $COSMIC_OUTFILE = "cosmic_fixed_no_silent.tab";
  } else {
    $COSMIC_OUTFILE = "cosmic_fixed_including_silent.tab";
  }
}

my %TUMOR_ORIGINS = map {$_, 1} (
				 "primary",
				 "metastasis",
				 "NS",
				 "secondary",
				 "recurrent",
				 "surgery fresh/frozen",
				 "hyperplasia adjacent to primary tumour",
				 "adenoma adjacent to primary tumour"
				);
# known states from report file

my %SOMATIC_STATUS = (
		      "Confirmed somatic variant" => 1,
		      "Reported in another cancer sample as somatic" => 1,
		      "Confirmed germline variant" => 0,
		      "Not specified" => 0,
		      "Reported in another sample as germline" => 0,
		      "Variant of unknown origin" => 0,
		     );

if ($FLAGS{"somatic-site-rescue"}) {
  somatic_site_rescue();
  exit(0);
} elsif ($FLAGS{"find-hypermutable-samples"}) {
  find_hypermutable_samples();
  exit(0);
} elsif ($FLAGS{"dump-gedi-diagnosis"}) {
  dump_gedi_diagnosis();
  exit(0);
} elsif ($FLAGS{"disease-hack"}) {
  disease_hack();
  exit(0);
} elsif ($FLAGS{"hyper-prep"}) {
  hypermutable_prep();
  exit(0);
}

printf STDERR "excluding silent mutations?: %s\n", $IGNORE_SILENT ? "yes" : "no";

my $vm_pcgp = get_vm_pcgp();

my $bad_literature = load_bad_literature();
my $bad_variants = load_bad_variants();

my $nm_lengths = load_nm_lengths();
my $gene2nm = load_gene2nm();

my $cgl = read_simple_file($FLAGS{"cancer-gene-list"} || die "-cancer-gene-list");
my %cancer_genes = map {$_, 1} @{$cgl};

my $badl = read_simple_file($FLAGS{"bad-gene-list"} || die "-bad-gene-list");
my %bad_genes = map {$_, 1} @{$badl};

my $vm_nhlbi_blacklist;
if ($FLAGS{"no-nhlbi"}) {
  printf STDERR "DEBUG, NHLBI blacklist disabled\n";
  $vm_nhlbi_blacklist = get_new_vm();
} else {
  $vm_nhlbi_blacklist = parse_nhlbi_blacklist();
}

#
#  start:
#
my %hypermutable_samples;
my %hypermutable_tid;
if ($FLAGS{"hypermutable-samples"}) {
  my $hs = read_simple_file($FLAGS{"hypermutable-samples"} || die "-hypermutable-samples");
  %hypermutable_samples = map {$_, 1} @{$hs};
} elsif (my $hf = $FLAGS{"hypermutable-tumor-config"}) {
  my $df = new DelimitedFile("-file" => $hf,
			     "-headers" => 1,
			     );
  while (my $row = $df->get_hash()) {
    my $tid = $row->{tumor_id} || die;
    my $call = $row->{hypermut_xma};
    $hypermutable_tid{$tid} = $call;
  }

  my $hyper_count = 0;
  foreach (values %hypermutable_tid) {
    $hyper_count++ if $_;
  }
  printf STDERR "hypermutable percent: %.1f\n", $hyper_count * 100 / scalar keys %hypermutable_tid;
} else {
  die "-hypermutable-tumor-config";
}

my $f_cosmic = $FLAGS{cosmic} || die "-cosmic";
die "where is $f_cosmic?" unless -s $f_cosmic;

my $df = new DelimitedFile(
			   "-file" => $f_cosmic,
			   "-headers" => 1,
			  );

my @all_rows;

my $line_count = 0;
my $max_lines = $FLAGS{"max-lines"};
my $null_positions = 0;

my $rpt_fixed = $df->get_reporter(
				  "-file" => $COSMIC_OUTFILE,
				  "-extra" => [
					       qw(
						   sj_diagnosis
						   sj_subtype
						   sj_subgroup
						)
					      ]
				 );

my $rpt_null = $df->get_reporter(
				 "-file" => "null_positions.tab",
				);

my %reject;
# under construction

my $F_POS;
my $F_STRAND;

my $ga;
if ($REANNOTATE_GENE_SYMBOLS) {
  if ($REANNOTATE_GENE_SYMBOLS == 2) {
    my $f = $config_genome->{REFSEQ_REFFLAT};
    $ga = new GeneAnnotation(
			     "-style" => "refgene_flatfile",
			     "-refgene_flatfile" => $f,
			     "-single_gene_mode" => 1,
			    );
  } else {
    my $ger_dir = $config_genome->{GENE_EXON_REGION_DIR} || die;
    $ga = new GeneAnnotation(
			     "-style" => "gene_exon_region",
			     "-gene_exon_region_dir" => $ger_dir,
			     "-single_gene_mode" => 1,
			    );
  }
}


my %reannot_note;
printf STDERR "start initial parse...\n";

my $tfw = new TemporaryFileWrangler();
if (0) {
  printf STDERR "DEBUG: not unlinking tempfiles!\n";
  $tfw->auto_unlink(0);
}

my %tabix_spill;

while (my $row = $df->get_hash()) {
#  my $tumor_origin = $row->{"Tumour origin"} || die;
#  die "unknown origin type $tumor_origin" unless $TUMOR_ORIGINS{$tumor_origin};
#  next unless $tumor_origin eq "primary";
  # leave this alone for now

#  printf STDERR "%d...\n", ++$line_count;

  unless ($F_POS) {
    $F_POS = get_position_column($row);
    $F_STRAND = get_strand_column($row);
  }

  my $pmid = $row->{Pubmed_PMID};
  if ($pmid and my $info = $bad_literature->{$pmid}) {
    if ($info->{all_bad}) {
      $info->{all_bad} = 2;
      $reject{bad_publication}++;
      printf STDERR "reject bad_publication\n";
      next;
    } else {
      my $gene = $row->{"Gene name"};
      if ($info->{bad_genes}{$gene}) {
	$reject{bad_publication_gene}++;
	$info->{bad_genes}{$gene} = 2;
	# actually found in data
	printf STDERR "reject bad_publication_and_gene\n";
	next;
      }
    }
  }

  if (my $gene = $row->{"Gene name"} and my $cds = $row->{"Mutation CDS"}) {
    if ($bad_variants->{$gene}{$cds}) {
      printf STDERR "reject bad_variant\n";
      $reject{bad_variant}++;
      $bad_variants->{$gene}{$cds} = 2;
      next;
    }
  }

  if ($row->{$F_POS} !~ /\w/) {
    $reject{null_position}++;
    printf STDERR "reject null_pos\n";
    $null_positions++;
    $rpt_null->end_row($row);
    next;
  }

  my $is_silent = 0;
  my $aa = $row->{"Mutation AA"};
  if ($aa =~ /^p\.([A-Z])\d+([A-Z])$/) {
    $is_silent = 1 if $1 eq $2;
  }
  my $desc = $row->{"Mutation Description"} || die;
  $is_silent =1 if $desc =~ /silent/i;
  # sometimes AA annotation is incomplete but description mentions
  # silent, e.g.:
  # ANO3    ENST00000256737 2946    14004   TCGA-13-1487-01 1474985 1398684 ovary   NS      carcinoma       serous_carcinoma        y       79500   c.2946G>A       p.*982* Substitution - coding silent    het     11:26638567-26638567    +

  # filter NHLBI variants:
  my $cds = $row->{"Mutation CDS"};
  if ($cds =~ /\d+([ACGT])>([ACGT])$/) {
    # simple SNV
    my ($ref_base, $var_base) = ($1, $2);
    my $strand = $row->{$F_STRAND};
    my $is_minus = $strand eq "-";
    if ($is_minus) {
      $ref_base = reverse_complement($ref_base);
      $var_base = reverse_complement($var_base);
    }

    my ($chrom, $pos) = get_chrom_and_pos($row, $F_POS);

    if ($vm_nhlbi_blacklist->find_snv(
				      "-reference" => $chrom,
				      "-base-number" => $pos,
				      "-reference-base" => $ref_base,
				      "-variant-base" => $var_base
				     )) {
      printf STDERR "reject NHLBI %s, strand=%s, aa=%s, silent=%d\n",
	join(".", $chrom, $pos, $ref_base, $var_base),
	  $strand, $aa, $is_silent;
      $reject{nhlbi}++;
      next;
    }
  }

  if ($IGNORE_SILENT and $is_silent) {
    printf STDERR "reject silent\n";
    $reject{silent}++;
    next;
    # put this filter after NHLBI check so we can count those events
  }

  my $somatic_status = $row->{"Mutation somatic status"};
  my $wanted = $SOMATIC_STATUS{$somatic_status};
  die "unknown somatic status" unless defined $wanted;
  my $rescue = $somatic_status eq "Confirmed somatic variant";

  unless ($wanted) {
    # require somatic variants only

    my $nsp = new NucleotideSubstitutionParser();
    my $cds_raw = $row->{"Mutation CDS"} || "";
    # can be blank in v72
    my $strand = $row->{$F_STRAND} || die;
    $nsp->auto_strand_fix($strand);

    if ($nsp->parse($cds_raw)) {
      my ($chrom, $pos) = get_chrom_and_pos($row, $F_POS);
      my $ref = $nsp->reference_sequence;
      my $var = $nsp->variant_sequence;
      if (length($ref) == length($var)) {
	# sub
	if (my $hits = $vm_pcgp->find_snv(
				      "-reference" => $chrom,
				      "-base-number" => $pos,
				      "-reference-base" => $ref,
				      "-variant-base" => $var
					 )) {
	  printf STDERR "salvage PCGP sub %s\n", join ".", $chrom, $pos, $ref, $var;
	  $rescue = 1;
	}
      } elsif ($ref eq "-" or $var eq "-") {
	# insertion or deletion
	if ($ref =~ /\S/ and $var =~ /\S/) {
	  # full allele details not always specified, e.g. c.924_925ins13
	  my $v = new Variant();
#	  printf STDERR "debug %s %s\n", $cds_raw, join ".", $chrom, $pos, $ref, $var;
	  $v->import_generic(
			     "-reference-name" => $chrom,
			     "-base-number" => $pos,
			     "-reference-allele" => $ref,
			     "-variant-allele" => $var,
			    );

	  if (my $hits = $vm_pcgp->find_indel(
					      "-variant" => $v,
					      "-match-basic-type" => 1,
					      "-match-size" => 1,
					      "-fuzz-bases" => INDEL_NUCLEOTIDE_NT_FUZZY_MATCH,
					     )) {
	    $rescue = 1;
	    printf STDERR "salvage PCGP indel %s\n", join ".", $chrom, $pos, $ref, $var;
	  }
	}
      } else {
#	dump_die($row, "ERROR: unhandled indel $cds_raw: $chrom $pos $ref $var", 1);
	printf STDERR "ERROR: unhandled indel $cds_raw: $chrom $pos $ref $var\n";
	# typically complex, code does not handle these well.
	# very few of these in PCGP anyway.
      }
    }

    unless ($rescue) {
      printf STDERR "reject non-somatic status %s\n", $somatic_status;
      $reject{non_somatic_status}++;
      next;
    }
  }

  if ($EXCLUDE_HYPERMUTABLE and not($rescue)) {
    #
    # filter hypermutable samples unless:
    # - validated PCGP variant
    # - marked as confirmed somatic
    #
    if (%hypermutable_samples) {
      if ($hypermutable_samples{$row->{ID_sample}}) {
	$reject{hypermutable}++;
	printf STDERR "reject hypermutable sample %s\n", $row->{ID_sample};
	next;
      }
    } elsif (%hypermutable_tid) {
      my $tid = $row->{ID_tumour};
      my $is_hyper = $hypermutable_tid{$tid};
#      dump_die($row, "no hyper call for tumor_id $tid") unless defined $is_hyper;
      # sanity check: needs reanalysis w/new builds
      # CAN'T: not all tumor ids reported ("Unknown" description, etc)
      if ($is_hyper) {
	$reject{hypermutable}++;
	printf STDERR "reject hypermutable tumor %s\n", $tid;
	next;
	# exclude variants from hypermutators
      }
    } else {
      die;
    }
  }

  add_disease_map($row);

  if ($REANNOTATE_GENE_SYMBOLS and $row->{$F_POS}) {
    my $gene_raw = $row->{"Gene name"} || die;
#    if (not($ga->is_valid_gene("-gene" => $gene_raw))) {
    # revised: just reannotate everything

    my $pos_string = $row->{$F_POS};

    my ($chr, $range) = split /:/, $pos_string;
    my ($start, $end) = split /\-/, $range;

    my $patched;

    if ($ga->find(
		  "-reference" => $chr,
		  "-start" => $start,
		  "-end" => $end
		 )
       ) {
      my $genes = $ga->results_genes();
      foreach (@{$genes}) {
	s/_loc\w+$//;
      }
      $patched = join ",", @{$genes};
      #      printf STDERR "patch %s => %s\n", $gene_raw, $patched;
    }

    unless ($patched) {
      my @f = split /_/, $gene_raw;
      if (@f == 2) {
	if ($f[1] eq "HUMAN") {
	  $patched = $f[0];
	} elsif ($f[1] =~ /^\d+$/ and $f[0] =~ /^[A-Z]+\d+$/) {
	  # possibly accession.version reformatted to accession_version?
	  # e.g. AB019437_1, 
	  $patched = join ".", @f;
	} elsif ($f[1] =~ /^(ENST\d+)/ and
		 $row->{"Accession Number"} =~ /^$1/) {
	  # appended transcript ID
	  # might not be a perfect match, e.g.
	  # AC073343.1_ENST00000544825 where accession is
	  # ENST00000544825_v68
	  $patched = $f[0];
	} elsif ($TRIM_ENSEMBL_ID_IF_ACCESSION_MISMATCHES and
		 $f[1] =~ /^ENST\d+/) {
	  # Gene name: KRBOX1_ENST00000418176
	  # Accession Number: ENST00000426937
	  # - transcript suffix doesn't match, however KRBOX1 is a valid symbol
	  $patched = $f[0];
	  #	  printf STDERR "hey now: %s %s %s\n", $gene_raw, $row->{"Accession Number"}, $patched;
	} else {
	  # 	  die join ",",$gene_raw, $row->{"Accession Number"};
	}
      }
    }

    unless ($patched) {
      printf STDERR "failed to find sym for %s at %s.%d\n",
	$gene_raw, $chr, $start;
      $patched = $gene_raw;
    }

    if ($row->{"Gene name"} ne $patched) {
      my $old = $row->{"Gene name"};
      my $key = join "_", $old, $patched;
      unless ($reannot_note{$key}) {
	printf STDERR "reannotate %s => %s\n", $old, $patched;
	$reannot_note{$key} = 1;
      }
      $row->{"Gene name"} = $patched;
    }
  }  # reannotate

  $rpt_fixed->end_row($row);

  last if $max_lines and ++$line_count > $max_lines;
}
$rpt_fixed->finish();
$rpt_null->finish();

printf STDERR "rejected record summary:\n";
foreach (sort keys %reject) {
  printf STDERR "  %s: %d\n", $_, $reject{$_};
}

#
#  confirm bad publication filtering:
#
foreach my $pmid (keys %{$bad_literature}) {
  my $info = $bad_literature->{$pmid};
  if ($info->{all_bad}) {
    printf STDERR "bad publication filter result for PMID %d: %s\n",
      $pmid, $info->{all_bad} == 2 ? "ok" : "missed!";
  } else {
    foreach my $bad_gene (keys %{$info->{bad_genes}}) {
      printf STDERR "bad publication filter result for PMID %d gene %s: %s\n", $pmid, $bad_gene, $info->{bad_genes}{$bad_gene} == 2 ? "ok" : "missed!";
      # 2 = actually saw in data
    }
  }
}

#
#  confirm bad variant filtering:
#
foreach my $gene (keys %{$bad_variants}) {
  foreach my $variant (keys %{$bad_variants->{$gene}}) {
    printf STDERR "bad variant filter result for %s %s: %s\n",
      $gene, $variant,
	$bad_variants->{$gene}{$variant} == 2 ? "ok" : "missed!";
    # 2 = actually saw in data
  }
}

unless ($FLAGS{"no-summary"}) {
  # report_summary_by_gene(\@all_rows);
  die "needs update, all rows no longer in RAM";
}

unless ($FLAGS{"no-recurrence"}) {
#  report_recurrence(\@all_rows);
  die "needs update, all rows no longer in RAM";
}



sub get_unique_variant_id {
  my ($row) = @_;
  my $vkey;
  if (0) {
    $vkey = $row->{"Mutation ID"} || die;
  } else {
    my @bits;
    foreach my $f ($F_POS, "Mutation CDS") {
      my $v = $row->{$f};
      $v = "NA" unless defined $v;
#      dump_die($row, "undef/empty field $f") unless $v and length($v) and $v =~ /\w/;
      # Mutation CDS may be blank in COSMIC v72
      push @bits, $v;
    }
    $vkey = join "_", @bits;
    # position alone might NOT be unique,
    # e.g. 17:7579472 has c.215G>C and c.215G>A
  }
  return $vkey;
}

sub report_recurrence {
  my ($all_rows) = @_;
  my %id2row;
  foreach my $row (@{$all_rows}) {
    my $vid = get_unique_variant_id($row);
    push @{$id2row{$vid}}, $row;
  }

  my $rpt = new Reporter(
			 "-file" => "recurrence.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       "Gene name",
				       "Mutation genome position",
				       "Mutation CDS",
				       "Mutation AA",

				       "variant_key",
				       "sample_count",
				       "aa_ambiguous",
				       "pubmed_count",
				      ]
			);

  foreach my $vkey (keys %id2row) {
    my $rows = $id2row{$vkey};
    my %aa = map {$_->{"Mutation AA"}, 1} @{$rows};
    my %sid = map {$_->{"ID_sample"}, 1} @{$rows};
    my %pmid = map {$_->{"Pubmed_PMID"}, 1} @{$rows};
    die unless %sid;
    die unless %pmid;

    printf STDERR "WARNING: ambiguous AA %s for %s\n", join(",", keys %aa), $vkey if scalar keys %aa > 1;

    my %r;
    %r = %{$id2row{$vkey}->[0]};

    $r{"Mutation genome position"} = $r{$F_POS};
    $r{"Mutation AA"} = join ",", sort keys %aa;
    $r{aa_ambiguous} = scalar keys %aa > 1 ? 1 : 0;
    $r{variant_key} = $vkey;
    $r{sample_count} = scalar keys %sid;
    $r{pubmed_count} = scalar keys %pmid;

    $rpt->end_row(\%r);
  }
  $rpt->finish();

  #  echo "$geneName       $AAchange       $pubmed_count   $occurence      $chrloc" >>${OUTPUT_FILE}
}

sub report_summary_by_gene {
  my ($all_rows) = @_;
  my %id2row;
  foreach my $row (@{$all_rows}) {
    my $vid = get_unique_variant_id($row);
    # FIX ME:
    # for this report,
    # should variants be counted only by AA change annotation??
    # - the problem with doing it this way is that AA annotation
    #   formats are less consistent than CDS annotations
    push @{$id2row{$vid}}, $row;
  }

  my %gene_info;

  foreach my $vkey (keys %id2row) {
    my $rows = $id2row{$vkey};
    my %genes = map {$_->{"Gene name"}, 1} @{$rows};

    my $glc = new GeneListCollapser();
    $glc->verbose(1);
    $glc->collapse("-hash" => \%genes);

    if (scalar keys %genes > 1) {
      my $msg = sprintf "ERROR: multiple genes for %s: %s", $vkey, join ",", sort keys %genes;
      $glc->collapse("-hash" => \%genes, "-single" => 1);
      if ($TOLERATE_NONSENSE) {
	printf STDERR "$msg\n";
      } else {
	die join ",", $msg;
      }
    }
    my ($gene) = keys %genes;
    my %pmid = map {$_->{"Pubmed_PMID"}, 1} @{$rows};
    my %aa = map {$_->{"Mutation AA"}, 1} @{$rows};
    my %samples = map {$_->{"ID_sample"}, 1} @{$rows};

    if (scalar keys %samples >= $RECURRENCE_SAMPLE_COUNT) {
      #
      # variant was observed recurrently
      #
      increment_tracker(\%gene_info, $gene, "recurrent_samples", \%samples);
      increment_tracker(\%gene_info, $gene, "recurrent_pmid", \%pmid);
      increment_tracker(\%gene_info, $gene, "recurrent_variants", $vkey);
      increment_tracker(\%gene_info, $gene, "recurrent_aa", \%aa);
      increment_aa_tracker(\%gene_info, $gene, "recurrent_aa_samples", \%aa, \%samples);
    }
    increment_tracker(\%gene_info, $gene, "all_variants", $vkey);
    increment_tracker(\%gene_info, $gene, "all_samples", \%samples);
    increment_tracker(\%gene_info, $gene, "all_pmid", \%pmid);
    increment_tracker(\%gene_info, $gene, "all_aa", \%aa);
    increment_aa_tracker(\%gene_info, $gene, "all_aa_samples", \%aa, \%samples);
  }

  if (0) {
    # debug: report a NON-recurrent variant for sanity checking
    foreach my $sid (keys %{$gene_info{TP53}{all_samples}}) {
      if ($gene_info{TP53}{recurrent_samples}{$sid}) {
	print STDERR "recurrent TP53 sample $sid\n";
      } else {
	print STDERR "non-recurrent TP53 sample $sid\n";
      }
    }
    die;
  }

  #
  #  set final counts for each category:
  #
  foreach my $gene (keys %gene_info) {
    foreach my $label (keys %{$gene_info{$gene}}) {
      my $cl = "count_" . $label;
      $gene_info{$gene}{$cl} = scalar keys %{$gene_info{$gene}{$label}};
    }
  }

  my $rpt = new Reporter(
			 "-file" => "gene_summary.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       "gene",
				       "cancer_gene",
				       "refseq",
				       "refseq_length",
				       "bad_gene",
				       "variant_count_recurrent",
				       "variant_count_all",
				       "sample_count_recurrent",
				       "sample_count_all",
				       "pubmed_count_recurrent",
				       "pubmed_count_all",
				       "aa_count_all",
				       "aa_count_recurrent",
				       "aa_count_stop_recurrent",
				       "aa_count_stop_all",
				       "aa_count_frameshift_recurrent",
				       "aa_count_frameshift_all",
				       "aa_recurrent",
				       "aa_all",
				      ]
			);

  foreach my $gene (sort {
    ($gene_info{$b}{count_recurrent_samples} || 0) <=>
      ($gene_info{$a}{count_recurrent_samples} || 0)
    } keys %gene_info) {
    # sort genes by highest count of recurrent samples
    my %r;
    $r{gene} = $gene;
    if (my $nm = $gene2nm->{$gene}) {
      $r{refseq} = $nm;
      $r{refseq_length} = $nm_lengths->{$nm} || "n/a";
    } else {
      $r{refseq} = "n/a";
      $r{refseq_length} = "n/a";
    }
    $r{cancer_gene} = $cancer_genes{$gene} ? 1 : 0;
    $r{bad_gene} = $bad_genes{$gene} ? 1 : 0;
    $r{variant_count_recurrent} = $gene_info{$gene}{count_recurrent_variants};
    $r{variant_count_all} = $gene_info{$gene}{count_all_variants};
    $r{sample_count_recurrent} = $gene_info{$gene}{count_recurrent_samples};
    $r{sample_count_all} = $gene_info{$gene}{count_all_samples};
    $r{pubmed_count_recurrent} = $gene_info{$gene}{count_recurrent_pmid};
    $r{pubmed_count_all} = $gene_info{$gene}{count_all_pmid};
    $r{aa_count_recurrent} = $gene_info{$gene}{count_recurrent_aa};
    $r{aa_count_all} = $gene_info{$gene}{count_all_aa};
    $r{aa_count_stop_recurrent} = count_stops($gene_info{$gene}{recurrent_aa});
    $r{aa_count_stop_all} = count_stops($gene_info{$gene}{all_aa});
    $r{aa_count_frameshift_recurrent} = count_frameshifts($gene_info{$gene}{recurrent_aa});
    $r{aa_count_frameshift_all} = count_frameshifts($gene_info{$gene}{all_aa});


#    $r{aa_recurrent} = join ",", sort keys %{$gene_info{$gene}{recurrent_aa}};
#    $r{aa_all} = join ",", sort keys %{$gene_info{$gene}{all_aa}};

    $r{aa_recurrent} = get_sorted_aa_list(\%gene_info, $gene, "recurrent_aa_samples");
    $r{aa_all} = get_sorted_aa_list(\%gene_info, $gene, "all_aa_samples");

    $rpt->end_row(\%r);
  }

  $rpt->finish();

  # gene
  # sample_count_recurrent
  # sample_count_all
  # pmid_recurrent
  # pmid_all
  # gene_class
  # cancer_gene
  # broken_gene

}

sub increment_tracker {
  my ($tracker, $gene, $label, $thing) = @_;
  if (ref $thing eq "HASH") {
    foreach my $key (keys %{$thing}) {
      $tracker->{$gene}{$label}{$key}++;
    }
  } elsif (ref $thing){
    die;
  } else {
    $tracker->{$gene}{$label}{$thing}++;
  }
}

sub increment_aa_tracker {
  my ($tracker, $gene, $label, $aas, $samples) = @_;
  foreach my $aa (keys %{$aas}) {
    foreach my $sample (keys %{$samples}) {
      $tracker->{$gene}{$label}{$aa}{$sample}++;
    }
  }
}

sub count_stops {
  # e.g. p.E1097*
  my ($hash) = @_;
  my $stops = 0;
  foreach my $aa (keys %{$hash}) {
    $stops++ if $aa =~ /p\.[A-Z]\d+\*$/;
  }
  return $stops;
}

sub count_frameshifts {
  # e.g. p.C917fs*8
  my ($hash) = @_;
  my $fs = 0;
  foreach my $aa (keys %{$hash}) {
    $fs++ if $aa =~ /p\.[A-Z]\d+fs\*/;
  }
  return $fs;
}

sub get_sorted_aa_list {
  my ($tracker, $gene, $label) = @_;

  my $root = $tracker->{$gene}{$label};
  my $result;

  if ($root) {
    my %aa2count;
    foreach my $aa (sort keys %{$root}) {
      $aa2count{$aa} = scalar keys %{$root->{$aa}};
    }
    my @aa_by_count = sort {$aa2count{$b} <=> $aa2count{$a}} keys %aa2count;
    my @things;
    foreach my $aa (@aa_by_count) {
      #    printf STDERR "%s: %d\n", $aa, $aa2count{$aa};
      push @things, sprintf "%s=%d", $aa, $aa2count{$aa};
      last if @things == $GENE_SUMMARY_MAX_AA_COUNT;
    }
    push @things, "..." if @aa_by_count > $GENE_SUMMARY_MAX_AA_COUNT;
    $result = join ", ", @things;
  }
  return $result;
}

sub load_gene2nm {
#  my $f = $FLAGS{gene2nm} || die "-gene2nm";
  my $f = $config_genome->{GENE_TRANSCRIPT_MATRIX} || die;
  open(IN, $f) || die;
  my %gene2nm;
  while (<IN>) {
    chomp;
    my ($gene, $dunno_maybe_gi, $nm) = split /\t/, $_;
    $gene2nm{$gene} = $nm unless $gene2nm{$gene};
  }
  return \%gene2nm;
}

sub load_nm_lengths {
  my $nm_fa = $FLAGS{"refgene-fasta"} || die "refgene-fasta";
  die $nm_fa unless -s $nm_fa;
  my $fh = universal_open($nm_fa) || die;
  my $nm;
  my %lengths;
  print STDERR "loading refGene lengths...\n";
  while (<$fh>) {
    chomp;
    if (/^>/) {
      #    if (/(NM_\d+)/) {
      if (/([NX][MR]_\d+)/) {
	# use unversioned
	$nm = $1;
	$lengths{$nm} = 0;
      } else {
	die "can't identify accession in $_";
      }
    } else {
      $lengths{$nm} += length($_);
    }
  }
  return \%lengths;
}

sub parse_nhlbi_blacklist {
  #
  #  database of sites we DON'T want to report
  #  (sites found at freqency of > .001+)
  #  (...taking it again from the top)
  #
  # copied from medal_ceremony.pl, ugh; make into a module??
  print STDERR "parsing NHLBI blacklist (mkII)...";
  my $blacklist = $FLAGS{"nhlbi-blacklist-mk2"} || die "-nhlbi-blacklist-mk2";
  my $vm = get_new_vm();

  my $df = new DelimitedFile("-file" => $blacklist,
			     "-headers" => 1,
			     );
  my $count = 0;
  while (my $row = $df->get_hash()) {
    print STDERR "." if ++$count % 25000 == 0;

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
  print STDERR "\n";

  return $vm;
}

sub get_new_vm {
#  my ($force_rsc) = @_;
#  my $rsc = $force_rsc || $REFERENCE_SANITY_CHECK;
#  my @options;
#  push @options, ("-fasta_dir" => $FLAGS{"fasta-dir"} || die "-fasta-dir") if $rsc;
  return new VariantMatcher();
}

sub load_bad_literature {
  my $infile = $FLAGS{"bad-literature"} || die "-bad-literature";

  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );

  my %bad_pmid;
  while (my $row = $df->get_hash()) {
    my $pmid = $row->{pmid} || die;
    my $gene = $row->{gene} || die;
    my $all_bad = $row->{all_bad};
    die unless defined $all_bad;

    $bad_pmid{$pmid}{all_bad} = $all_bad;
    $bad_pmid{$pmid}{bad_genes}{$gene} = 1;
  }
  return \%bad_pmid;
}

sub load_bad_variants {
  my $infile = $FLAGS{"bad-variants"} || die "-bad-variants";
  open(IN, $infile) || die;
  my %bad;
  while (<IN>) {
    chomp;
    my @f = split /\s+/, $_;
    die unless @f == 3;
    # wack format: tab + space
    my ($gene, $change, $whatevs) = @f;
    $bad{$gene}{$change} = 1;
  }
  return \%bad;
}

sub somatic_site_rescue {
  my $f_cosmic = $FLAGS{cosmic} || die "-cosmic";
  my $f_final = $FLAGS{final} || die "-final";

  #
  # final cleaned results file:
  #
  my $df = new DelimitedFile(
			     "-file" => $f_final,
			     "-headers" => 1,
			     );
  my %final;

  my $f_pos;

  while (my $row = $df->get_hash()) {
    $f_pos = get_position_column($row) unless $f_pos;

    my $transcript = $row->{"Accession Number"};
    # RAW transcript unaffected by reannotation
    my $cds = $row->{"Mutation CDS"};
    my $pos = $row->{$f_pos};
    if ($transcript and $cds and $pos) {
      $final{$transcript}{$cds}{$pos} = 1;
    }
  }

  #
  #  find sites rejected for somatic status but at same site
  #  as sites in final file
  #
  $df = new DelimitedFile(
			  "-file" => ($FLAGS{cosmic} || die "-cosmic"),
			  "-headers" => 1,
			 );
  my %salvage;
  while (my $row = $df->get_hash()) {
    my $somatic_status = $row->{"Mutation somatic status"};
    my $wanted = $SOMATIC_STATUS{$somatic_status};
    die "unknown somatic status" unless defined $wanted;
    unless ($wanted) {
      # require somatic variants only
      my $transcript = $row->{"Accession Number"};
      my $cds = $row->{"Mutation CDS"};
      my $pos = $row->{$f_pos};
      if ($transcript and $cds and $pos) {
	if ($final{$transcript}{$cds}{$pos}) {
	  my $key = join "_", $transcript, $cds, $pos;
	  $salvage{$key}{transcript} = $transcript;
	  $salvage{$key}{cds} = $cds;
	  $salvage{$key}{pos} = $pos;
	  $salvage{$key}{count}++;
	  $salvage{$key}{status}{$somatic_status}++;
	}
      }
    }
  }

  my $rpt = new Reporter(
			 "-file" => "rescue_somatic_filter.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   transcript
					   cds
					   pos
					   count
					   status
					)
				      ],
			 "-auto_qc" => 1,
			);

  foreach my $key (sort keys %salvage) {
    my $info = $salvage{$key};
    my %r = %{$info};

    my @stat;
    foreach my $status (sort keys %{$info->{status}}) {
      push @stat, join "=", $status, $info->{status}{$status};
    }
    $r{status} = join ",", @stat;
    $rpt->end_row(\%r);
  }
  $rpt->finish();

}

sub find_column {
  my ($row, $try) = @_;
  my $result;
  foreach my $f (@{$try}) {
    if (exists $row->{$f}) {
      $result = $f;
      last;
    }
  }
  return $result;
}

sub get_position_column {
  my ($row) = @_;
  my @try = (
	     "Mutation GRCh37 genome position",
	     # v66
	     "Mutation genome position",
	     # v72 backported to GRCh37
	    );
  return find_column($row, \@try);
}

sub get_strand_column {
  my ($row) = @_;
  my @try = (
	     "Mutation GRCh37 strand",
	     # v66
	     "Mutation strand",
	     # v72 backported to GRCh37
	    );
  return find_column($row, \@try);
}

sub find_hypermutable_samples {
  # reimplementation of JZ's process, described as:
  #
  # "To mitigate the impact of hymutable samples on global mutation
  # profile, hypermutable tumors were determined by the following process
  # and their mutations were filtered. Initial hypermutable samples were
  # determined by requiring # of coding mutations >=100 and we then check
  # whether the tumor was ranked amongst top 10% highest mutated tumors of
  # their published study. If >10% of the tumors were classified by this
  # threshold, the threshold will be re-adjusted to reflect higher
  # background mutation rate in this process. A total of 448 tumors were
  # identified as hypermutable by this process."
  #
  # QUESTIONS:
  # 1. intital count of >= 100: within a study or global?
  # 2. is tracking by SAMPLE ID or TUMOR ID?  TUMOR, even though sample
  #    ID is final value reported?
  # 3. do does tracking include silent?  yes, as in coding region?
  #    count all rows?
  # 4. how to track study, PMID?  What if blank?
  # 5. any duplicate handling?
  # 6. TID 1524187: high mutation count but only 8 tumors in the study.
  #    only 1 study with this tumor.  Discarded since not in 10% of
  #    PMID.  Rules don't work well if study reports small absolute
  #    number of samples.
  # 7. "if >10% of the tumors..." what is this?  final check at end?
  my %jz_hyper;
  if (my $f = $FLAGS{"hypermutable-samples"}) {
    my $hs = read_simple_file($f);
    die unless @{$hs};
    %jz_hyper = map {$_, 1} @{$hs};
  }

  my $f_cosmic = $FLAGS{cosmic} || die "-cosmic";

  # 1. count all coding mutations, bucketing by PMID and tumor ID/sample ID?
  # 2. foreach PMID:
  # 3. rank all tumors within PMID by mutation count
  # 4. examine all tumors above seed threshold (100)
  # 5. is tumor in top 10%?  if so, hypermutable
  #    output is SAMPLE ID not TUMOR ID (???)
  # compare count of hyper tumor IDs vs. all tumor IDs: if > 10% called
  # hyper, threshold needs to be adjusted.

  $df = new DelimitedFile(
			  "-file" => $f_cosmic,
			  "-headers" => 1,
			 );

  my %descs_wanted = (
		      "Substitution - Nonsense" => 1,
		      "Substitution - Missense" => 1,
		      "Substitution - coding silent" => 1,
		      "Deletion - Frameshift" => 1,

		      "Complex" => 1,
		      # ?
		      "Nonstop extension" => 1,

		      "Insertion - Frameshift" => 1,
		      "Deletion - In frame" => 1,
		      "Insertion - In frame" => 1,
		      "Complex - compound substitution" => 1,
		      "Complex - deletion inframe" => 1,
		      "Complex - insertion inframe" => 1,
		      "Complex - frameshift" => 1,
		      "Whole gene deletion" => 1,

		      "Unknown" => 0,
		      "No detectable mRNA/protein" => 0,
		     );

  my %pmid;
  my %tumor2mutation;
  my %tumor2pmid;

  #
  #  parse database, accumulating variants by publication/tumor/sample:
  #
  while (my $row = $df->get_hash()) {
    my $desc = $row->{"Mutation Description"} || "";
    my $wanted = $descs_wanted{$desc};
    die "unknown desc $desc" unless defined $wanted;

    if ($wanted) {
      my $id_sample = gimme($row, "ID_sample");
      my $id_tumor = gimme($row, "ID_tumour");
      my $id_mutation = gimme($row, "Mutation ID");
      my $pmid = gimme($row, "Pubmed_PMID");
      if ($pmid !~ /\w/ and $row->{ID_STUDY}) {
	# proxy if no PMID
	$pmid = sprintf "STUDY_%s", $row->{ID_STUDY};
      }
      unless ($pmid =~ /\w/) {
	# proxy if no PMID or STUDY ID
	my $hist_primary = gimme($row, "Primary histology");
	my $hist_subtype = gimme($row, $F_C_HISTOLOGY_SUBTYPE);
	$pmid = join ".", $hist_primary, $hist_subtype;
      }

      if ($pmid =~ /\w/) {
	printf STDERR "note: duplicate mutation %s in tumor ID %s\n", $id_mutation, $id_tumor if $tumor2mutation{$id_tumor}{$id_mutation};
	$tumor2mutation{$id_tumor}{$id_mutation}++;
	# hash rather than simple count:
	# a mutation may be reported more than once per tumor
	# (i.e. in different samples)
	$pmid{$pmid}{$id_tumor}{$id_sample}{$id_mutation}++;
	$tumor2pmid{$id_tumor}{$pmid} = 1;
      } else {
	dump_die($row, "no PMID or proxy");
      }
    }
  }

  my $verbose = 1;
  my %hyper_tid;
  my %hyper_samples;

  foreach my $tid (keys %tumor2mutation) {
    my $pmids = $tumor2pmid{$tid} || die "no PMIDs for TID $tid";
    my $tumor_count = scalar keys %{$tumor2mutation{$tid}};
    if ($tumor_count >= $HYPER_MIN_CODING_TO_SEED) {
      # tumor is a potential hypermutator
      printf STDERR "candidate hypermutator, tumor_id:%d, count:%d, study_count:%d (%s)\n",
	$tid,
	  $tumor_count,
	    scalar(keys %{$pmids}),
	      join(",", sort keys %{$pmids});

      my $is_hyper;
      foreach my $pmid (sort keys %{$pmids}) {
	#
	#  within each publication, rank tumors by mutation count:
	#
	my %counts;
	foreach my $tid (keys %{$pmid{$pmid}}) {
	  my $samples = $pmid{$pmid}{$tid};
	  my %mid;
	  foreach my $sid (keys %{$samples}) {
	    foreach my $mid (keys %{$samples->{$sid}}) {
	      $mid{$mid} = 1;
	    }
	  }
	  my $count = scalar keys %mid;
	  $counts{$tid} = $count;
	}

	my @sorted = sort {$counts{$b} <=> $counts{$a}} keys %counts;
	my $max_idx = int((scalar @sorted - 1) * $HYPER_TOP_FRACTION_TO_INCLUDE);
	# -1 to convert to index space,
	# e.g. for array of 100, last index to use is 9 (0-9)

	#
	#  top X percent to search for hypermutators
	#
	if ($verbose) {
	  printf STDERR "PMID:%s tumors:%d max_idx:%d\n",
	    $pmid, scalar @sorted, $max_idx;
	  my $idx = 0;
	  foreach my $t (@sorted) {
	    printf STDERR "  idx:%d tid:%d count:%d samples:%s%s\n",
	      $idx++,
		$t,
		  $counts{$t},
		    join(",", sort keys %{$pmid{$pmid}{$t}}),
		      ($t eq $tid ? " ***" : "");
	  }
	}

	my $max_value = $counts{$sorted[$max_idx]};

	for (my $i = 0; $counts{$sorted[$i]} >= $max_value; $i++) {
	  # rather than stopping at $max_idx, use the value
	  # at that position in case there are multiple entries
	  # with that value at the border
	  if ($sorted[$i] == $tid) {
	    $is_hyper = 1;
	    my @samples = keys %{$pmid{$pmid}{$tid}};
	    printf STDERR "hypermutable: tumor=%d samples=%s\n", $tid, join ",", @samples if $verbose;
	    $hyper_tid{$tid} = 1;
	    printf STDERR "NOTE: hypermutable tumor %d has %d samples\n", $tid, scalar @samples if @samples > 1;
	    foreach (@samples) {
	      $hyper_samples{$_} = 1;
	    }
	  }
	}
	print STDERR "\n" if $verbose;
      }
    } else {
      my %samples;
      my @pmids;
      foreach my $pmid (sort keys %{$pmids}) {
	push @pmids, $pmid;
	my $samples = $pmid{$pmid}{$tid} || die;
	foreach (keys %{$samples}) {
	  $samples{$_} = 1;
	}
      }
      printf STDERR "not considering tumor_id:%s count:%d pmids:%s samples:%s\n", $tid, $tumor_count, join(",", @pmids), join(",", sort keys %samples) if $verbose;
    }
  }

  my %all_tid;

  foreach my $pmid (keys %pmid) {
    foreach my $tid (keys %{$pmid{$pmid}}) {
      $all_tid{$tid} = 1;
    }
  }

  my $total_tumor_count = scalar keys %all_tid;
  my $hyper_tumor_count = scalar keys %hyper_tid;
  my $ratio = $hyper_tumor_count / $total_tumor_count;
  printf STDERR "total tumors:%d hypermutable:%d ratio:%.2f\n",
    $total_tumor_count, $hyper_tumor_count, $ratio;

  if (%jz_hyper) {
    # compare results vs. JZ's
    my %all = (%jz_hyper, %hyper_samples);
    foreach my $sid (sort {$a <=> $b} keys %all) {
      my $label;
      if ($hyper_samples{$sid} and $jz_hyper{$sid}) {
	$label = "both";
      } elsif ($hyper_samples{$sid} and !$jz_hyper{$sid}) {
	$label = "mne_exclusive";
      } elsif (!$hyper_samples{$sid} and $jz_hyper{$sid}) {
	$label = "jz_exclusive";
      } else {
	die;
      }

      printf STDERR "compare %s: %s\n", $sid, $label;
    }
  }

  open(OUT, ">hypermutable_samples.txt") || die;
  foreach (sort {$a <=> $b} keys %hyper_samples) {
    printf OUT "%s\n", $_;
  }

  open(OUT, ">hypermutable_tumors.txt") || die;
  foreach (sort {$a <=> $b} keys %hyper_tid) {
    printf OUT "%s\n", $_;
  }
  close OUT;

  die "ERROR: total hypermutable fraction exceeds max threshold!" if $ratio > $HYPER_MAX_TUMOR_FRACTION_HYPER;
  # adjust for higher background rate?



}

sub gimme {
  my ($row, $key) = @_;
  my $thing = $row->{$key};
  dump_die($row, "no value for $key") unless defined $thing;
  return $thing;
}

sub disease_hack {
  # development for COSMIC->SJ disease mapping

  my $infile = "histology.txt";
  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			     );
  my $outfile = $infile . ".disease.tab";
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   qw(
					       sj_diagnosis
					       sj_subtype
					       sj_subgroup
					    )
					  ]
			     );

  while (my $row = $df->get_hash()) {
    add_disease_map($row);
    $rpt->end_row($row);
  }
  $rpt->finish();
}

sub add_disease_map {
  # attempt to translate COSMIC disease annotations to SJ
  my ($row) = @_;
  my $c_primary = gimme($row, "Primary histology");
  my $c_subtype = gimme($row, $F_C_HISTOLOGY_SUBTYPE);

  my $sj_diagnosis = "";
  my $sj_subtype = "";
  my $sj_subgroup = "";
  # see gedi diagnosis/subtype/subgroup tables

  my %paired_map;
  $paired_map{"primitive_neuroectodermal_tumour-medulloblastoma"}{"SHH_subtype"} = [ qw(MB MB SHH) ];
  $paired_map{"primitive_neuroectodermal_tumour-medulloblastoma"}{"WNT_subtype"} = [ qw(MB MB WNT) ];

  $paired_map{"haematopoietic_neoplasm"}{"acute_myeloid_leukaemia"} = [ qw(AML AML) ];
  $paired_map{"haematopoietic_neoplasm"}{"acute_myeloid_leukaemia_therapy_related"} = [ qw(AML TAML) ];
  $paired_map{"haematopoietic_neoplasm"}{"acute_myeloid_leukaemia_myelodysplastic_syndrome_therapy_related_NOS"} = [ qw(AML TAML) ];

  my %primary_map = (
		     "osteosarcoma" => "OS",
		     "adrenal_cortical_carcinoma" => "ACT",

		     "Ewings_sarcoma-peripheral_primitive_neuroectodermal_tumour" => "EWS",
		     # not sure about this one: wikipedia:
		     # "peripheral primitive neuroectodermal tumours are generally not associated with bones, while Ewing sarcomas are most commonly related to bone."

		     "rhabdomyosarcoma" => "RHB",
		     "retinoblastoma" => "RB",
		     "neuroblastoma" => "NBL",
		     "malignant_melanoma" => "MEL",
		     # NOTE: "malignant_melanoma_of_soft_parts-clear_cell_sarcoma" should NOT be considered a melanoma:
		     #
		     # Alberto Pappo:
		     # "Malignant melanoma of soft parts is not a melanoma it is better known as clear cell sarcoma of soft parts.
		     # It has an EWS/CREB fusion and it is a soft tissue sarcoma I would definitely exclude"

		     "primitive_neuroectodermal_tumour-medulloblastoma" => "MB",
		    );

  my %subtype_map  = (
#		      "acute_lymphoblastic_leukaemia" => "ALL",
# no such SJ
		      "acute_lymphoblastic_T_cell_leukaemia" => "TALL",
		      "choroid_plexus_carcinoma" => "CPC",
		      "acute_lymphoblastic_B_cell_leukaemia" => "BALL",
		     );

  my %subtype_regexp = (
			"ependymoma" => "EPD",
			# fuzzier search required, e.g.
			#  - ependymoma
			#  - ependymoma_Grade_II
			#  - ependymoma_Grade_III-IV
		       );

  my %complex = (
		 "glioma" => "GG",
		);

  unless ($sj_diagnosis) {
    if (my $stub = $complex{$c_primary}) {
      my %grade2type = (
		       "I" => "L",
		       "II" => "L",
		       "III" => "H",
		       "IV" => "H",
		       "III-IV" => "H",
		      );
      if ($c_subtype =~ /grade_(.*)$/i) {
	my $grade = $1;
	my $type = $grade2type{$grade} || die "can't get type for grade $grade";
	my $code = sprintf '%s%s', $type, $stub;
	$sj_diagnosis = $sj_subtype = $code;
      }
    }
  }

  unless ($sj_diagnosis) {
    if (my $hit = $paired_map{$c_primary}{$c_subtype}) {
      ($sj_diagnosis, $sj_subtype, $sj_subgroup) = @{$hit};
      foreach ($sj_diagnosis, $sj_subtype, $sj_subgroup) {
	$_ = "" unless defined $_;
	# e.g. no subgroup beyond AML -> TAML
      }
    }
  }

  unless ($sj_diagnosis) {
    if (my $code = $primary_map{$c_primary} || $subtype_map{$c_subtype}) {
      # simple mapping where primary or subtype code is enough to
      # identify SJ diagnosis and subtype codes
      $sj_diagnosis = $sj_subtype = $code;
    }
  }

  unless ($sj_diagnosis) {
    foreach my $regexp (sort keys %subtype_regexp) {
      if ($c_subtype =~ /$regexp/i) {
	$sj_diagnosis = $sj_subtype = $subtype_regexp{$regexp};
	last;
      }
    }
  }


  $row->{sj_diagnosis} = $sj_diagnosis;
  $row->{sj_subtype} = $sj_subtype;
  $row->{sj_subgroup} = $sj_subgroup;
}

sub dump_gedi_diagnosis {
  my $dbi = get_dbi_gedi("-type" => "research");
  my $sql = 'select * from diagnosis,subtype,subgroup where diagnosis.diagnosis_id=subtype.diagnosis_id and subtype.subtype_id = subgroup.subtype_id;';
  export_query_to_flatfile(
			   "-dbi" => $dbi,
			   "-sql" => $sql,
			   "-outfile" => "gedi_diagnosis.tab",
			  );
  $dbi->disconnect();
}

sub dump_gedi_validated {
  my $dbi = get_dbi_gedi("-type" => "research");
  my $sql = 'select * from snpindel_find_t1';
  export_query_to_flatfile(
			   "-dbi" => $dbi,
			   "-sql" => $sql,
			   "-outfile" => "gedi_validated.tab",
			  );
  $dbi->disconnect();

}

sub hypermutable_prep {
# I will need to obtain a data matrix with following columns:
#  1. SampleID
#  2. Number_of_Tier1_Coding_SNV
#  3. Disease_Type
#  4. WGS_or_WES
#  5. prior_def_whether_sample_is_hypermutator
# A few optional columns would help also
#  6. Number_of_tier2/3_SNV (if WGS)
#  7. Age info (pediatric or adult)


  my $ht = $FLAGS{"hyper-prep"} || die;

  my $hyper_tid = read_simple_file($ht, "-hash1" => 1);

  my $f_cosmic = $FLAGS{cosmic} || die "-cosmic";
  my $ml = $FLAGS{"max-lines"};

  $df = new DelimitedFile(
			  "-file" => $f_cosmic,
			  "-headers" => 1,
			 );

  my %descs_wanted = (
		      "Substitution - Nonsense" => 1,
		      "Substitution - Missense" => 1,
		      "Substitution - coding silent" => 1,
		      "Deletion - Frameshift" => 1,

		      "Complex" => 1,
		      # ?
		      "Nonstop extension" => 1,

		      "Insertion - Frameshift" => 1,
		      "Deletion - In frame" => 1,
		      "Insertion - In frame" => 1,
		      "Complex - compound substitution" => 1,
		      "Complex - deletion inframe" => 1,
		      "Complex - insertion inframe" => 1,
		      "Complex - frameshift" => 1,
		      "Whole gene deletion" => 1,

		      "Unknown" => 0,
		      "No detectable mRNA/protein" => 0,
		     );

  my %tumor2mutation;
  my %tumor2disease;
  my %tumor2screen;
  my %tumor2age;

  #
  #  digest database:
  #
  my $row_count = 0;
  while (my $row = $df->get_hash()) {
    $row_count++;
    last if $ml and $row_count > $ml;
    my $desc = $row->{"Mutation Description"} || "";
    my $wanted = $descs_wanted{$desc};
    die "unknown desc $desc" unless defined $wanted;

    if ($wanted) {
      my $id_tumor = gimme($row, "ID_tumour");
      my $id_mutation = gimme($row, "Mutation ID");
      my $histology = gimme($row, "Primary histology");
      my $subtype = gimme($row, $F_C_HISTOLOGY_SUBTYPE);

      my $disease = join "/", $histology, $subtype;
      # hack
      $tumor2disease{$id_tumor}{$disease} = 1;
      $tumor2mutation{$id_tumor}{$id_mutation}++;
      # record generated tier here?
      #
      # hash rather than simple count:
      # a mutation may be reported more than once per tumor
      # (i.e. in different samples)

      my $screen = gimme($row, "Genome-wide screen");
      if ($screen eq "y") {
	$screen = 1;
      } elsif ($screen eq "n") {
	$screen = 0;
      } else {
	die;
      }
      $tumor2screen{$id_tumor}{$screen} = 1;

      my $age = gimme($row, "Age");
      $tumor2age{$id_tumor}{$age} = 1;
    }
  }

  # sanity checks:
  foreach my $tid (keys %tumor2disease) {
    die "multiple diseases for $tid" if scalar keys %{$tumor2disease{$tid}} > 1;
  }
  foreach my $tid (keys %tumor2screen) {
    printf STDERR "NOTE: multiple screens for TID %s\n", $tid
      if scalar keys %{$tumor2screen{$tid}} > 1;
    # very rarely more than one, e.g. TID 1656777 has
    # records for each type (maybe e.g. combination of WGS +
    # validation/sanger followup?)
  }
  foreach my $tid (keys %tumor2age) {
    die "multiple age for $tid" if scalar keys %{$tumor2age{$tid}} > 1;
  }

  my $rpt = new Reporter(
			 "-file" => "hyper_prep.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   tumor_id
					   mutation_count
					   primary_histology
					   histology_subtype
					   genome_wide_screen
					   hypermutator_call
					   age
					   sj_diagnosis
					   sj_subtype
					   sj_subgroup
					)
				      ],
			 "-auto_qc" => 1,
			);

  foreach my $tid (sort {$a <=> $b} keys %tumor2mutation) {
    my %r;
    $r{tumor_id} = $tid;
    $r{mutation_count} = scalar keys %{$tumor2mutation{$tid}};

    my @diseases = keys %{$tumor2disease{$tid}};
    die unless @diseases == 1;
    my ($d1, $d2) = split /\//, $diseases[0];

    $r{primary_histology} = $d1;
    $r{histology_subtype} = $d2;
    $r{"Primary histology"} = $d1;
    $r{$F_C_HISTOLOGY_SUBTYPE} = $d2;
    add_disease_map(\%r);

    my $screens = $tumor2screen{$tid};
    $r{genome_wide_screen} = $screens->{1} ? 1 : 0;
    # to handle > 1 screen type case

    $r{hypermutator_call} = $hyper_tid->{$tid} ? 1 : 0;

    my $ages = $tumor2age{$tid};
    my @ages = keys %{$ages};
    die unless @ages;
    $r{age} = join ",", @ages;

    $rpt->end_row(\%r);
  }
  $rpt->finish();
}

sub reject_summary {
  my ($file) = @_;
  open(IN, $file) || die;
  my %saw;
  foreach my $r (qw(
		     non_somatic_status
		     bad_variant
		     NHLBI
		     hypermutable_sample
		     hypermutable_tumor
		     bad_publication
		     bad_publication_and_gene
		     null_pos
		  )) {
    $saw{$r} = 1;
  }

  my %counts;
  while (<IN>) {
    if (/^reject (.*)/) {
      my $reason = $1;
      $reason = "non_somatic_status" if $reason =~ /^non-somatic status/;
      $reason = "NHLBI" if $reason =~ /^NHLBI/;
      $reason = "hypermutable_sample" if $reason =~ /hypermutable sample/;
      $reason = "hypermutable_tumor" if $reason =~ /hypermutable tumor/;
      die $reason unless $saw{$reason};
      $counts{$reason}++;
    }
  }
  foreach (sort keys %counts) {
    printf "%s: %d\n", $_, $counts{$_};
  }
}

sub get_vm_pcgp {
  my $f_val = $FLAGS{"gedi-validated"} || die "-gedi-validated";
  # TO DO: live db query?
  my $df = new DelimitedFile("-file" => $f_val,
			     "-headers" => 1,
			     );
  my $vm = new VariantMatcher();
  while (my $row = $df->get_hash()) {
    $vm->add_gedi_row("-row" => $row);
  }
  return $vm;
}

sub get_chrom_and_pos {
  my ($row, $f_pos) = @_;
  my $pos_raw = $row->{$f_pos};
  $pos_raw =~ s/\-\d+$//;
  my ($chrom, $pos) = split /:/, $pos_raw;
  die unless $pos =~ /^\d+$/;
  return ($chrom, $pos);
}

sub get_cooked_site {
  my ($row, $f_pos) = @_;
  my $pos_raw = $row->{$f_pos};
  my @f = split /:/, $pos_raw;
  die unless @f == 2;
  my ($chr_raw, $range) = @f;
  @f = split /\-/, $range;
  die unless @f == 2;
  my ($start, $end) = @f;
  return(cook_chromosome_name($chr_raw), $start, $end);
}

sub hack_test {
  my $file = "cosmic_fixed_including_silent.tab.gz";
  my $tf = new TabixFile("-file" => $file,
			 "-index" => 1,
			 "-f_chr" => $F_TABIX_CHROM,
			 "-f_start" => $F_TABIX_START,
			 "-f_end" => $F_TABIX_END,
			);
  die "done";
}

sub clean2tabix {
  my ($infile) = @_;

  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );

#  my @headers = @{$df->headers_raw};
#  push @headers, $F_TABIX_CHROM, $F_TABIX_START, $F_TABIX_END;

  my @headers = ($F_TABIX_CHROM, $F_TABIX_START, $F_TABIX_END);
  push @headers, @{$df->headers_raw};
  # "sort" command apparently confused by content and column
  # indexes don't work properly if position columns appended to the
  # output.  Prepend to work around.

  my $outfile = sprintf '%s.tabix.gz', basename($infile);

  my $tp = new TabixPrep(
			 "-outfile" => $outfile,
			 "-headers" => \@headers,
			 "-header_chr" => $F_TABIX_CHROM,
			 "-header_start" => $F_TABIX_START,
			 "-header_end" => $F_TABIX_END,
			);
  if (0) {
    print STDERR "DEBUG: no scratch/tempclean\n";
    $tp->delete_tempfiles(0);
    $tp->use_scratch(0);
  }


  while (my $row = $df->get_hash()) {
    unless ($F_POS) {
      $F_POS = get_position_column($row);
    }

    my ($chrom, $start, $end) = get_cooked_site($row, $F_POS);
    $row->{$F_TABIX_CHROM} = $chrom;
    $row->{$F_TABIX_START} = $start;
    $row->{$F_TABIX_END} = $end;

    $tp->add_row("-row" => $row);
  }
  $tp->finish();

}

sub append_occurrences {
  # add a new column to cleaned file of occurrence count for each variant.
  my ($infile) = @_;
  my $of = $infile . ".occurrences.gz";

  $df = new DelimitedFileHP(
			     "-file" => $infile,
			    );
  $df->prepare_query("-fields" => [qw(Chr WU_HG19_Pos ReferenceAllele MutantAllele)]);

  # first pass:
  my %counts;
  my $count = 0;
  while ($df->next_row()) {
    my $key = join ".", @{$df->get_query()};
    $counts{$key}++;
    printf STDERR "%d...\n", $count if ++$count % 100000 == 0;
  }

  $df = new DelimitedFileHP(
			     "-file" => $infile,
			     "-headers_extra" => [
						  "occurrence_count"
						 ]
			    );
  $df->write_init("-file" => $of, "-compress" => 1);
  $df->prepare_query("-fields" => [qw(Chr WU_HG19_Pos ReferenceAllele MutantAllele)]);

  $count = 0;
  while ($df->next_row()) {
    my $key = join ".", @{$df->get_query()};
    my %extra;
    $extra{occurrence_count} = $counts{$key} || die "no count for $key";
    $df->write_row("-extra" => \%extra);
    printf STDERR "%d...\n", $count if ++$count % 100000 == 0;
  }
  $df->write_finish();

}


