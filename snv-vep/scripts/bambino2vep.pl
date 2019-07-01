#!/bin/env perl
# use VEP to produce SJ-post-ish output
# MNE 6/2017+
#
# TO DO:
# - frameshift SJ AA format
# - splice_region AA annotations for non-coding
# - generate AAchange/exon annotation for coding_sequence_variant2.tab.hgvs.tab
# - manually check proteinIns, e.g. protein_altering2.tab.hgvs.tab
#   bambino-style insertions correct?  I think so, but many exceptions?
# - clarify: EXON and INTRON: continue uppercase?
# - break very large files into smaller chunks?
# - 1.2494330.G.A: hits for both NM_003820.2 and NM_003820.3, make this a
#   settable param

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Cwd qw(abs_path getcwd);
use Getopt::Long;
use File::Basename;
use File::Copy;
use File::Find;
use File::Path;
use List::Util qw(min);
use Data::Dumper;
use FileUtils qw(universal_open write_simple_file);
use Bio::SeqIO;

use TdtConfig;
use VEP;
use MiscUtils qw(dump_die shell_cmd);
use SJPreferredIsoform;
use MiscUtils qw(dump_die build_argv_list get_hash_option log_message);
use DelimitedFile;
use Reporter;
use Variant;
use VariantParsePolicy;
use AAParser;
use RefFlatFile qw(
		    CICERO_CODING
		    CICERO_UTR_3
		    CICERO_INTERGENIC
		 );
use ConfigUtils qw(config_or_manual);
use GeneAnnotation;
use GeneSymbolMapper qw(new_gsm_lite);
use TemporaryFileWrangler;
use EUtilities;

dump_die(\%INC) if $ENV{DEBUG_INC};

my $REFFLAT_CACHE_LIMIT = 25;

my $SPLICE_REGION_MAX_BP_FROM_EXON = 10;
my $SPLICE_MAX_BP_FROM_EXON = 3;
# narrow search when adding exon annotations to variants Annovar
# itself calls as splice
my $SPLICE_PROMOTION_MAX_BP_DONOR = 2;
my $SPLICE_PROMOTION_MAX_BP_ACCEPTOR = 3;
# SJ splice rules: within 2 nt of donor site and 3 nt of acceptor site.
# Annovar defaults are 2 nt from either, so special handling
# needed here.

my $UNKNOWN_AA_CHAR = "?";
# code to report when amino acid code is not provided (i.e. VEP "-").

my $CODON_ANNOTATION_UPSTREAM = 0;
# if 0, use the nearest codon, which is also used for the exon number.
# if 1, search UPSTREAM for the previous codon, which may be far away.
#       exon number still uses the closest exon.  This might be the
#       policy followed by SJ post but I don't know for sure.

my $SANITY_CHECK_AA_WITH_REFFLAT = 1;
my $SANITY_CHECK_AA_MIN_DISTANCE_TO_REPORT = 3;
#my $SANITY_CHECK_AA_MIN_DISTANCE_TO_REPORT = 1;
# small legitimate changes are possible due to e.g. "dup" annotation,
# so only report discrepancies of at least this many codon #s away

my $VERBOSE;

use constant VEP_NULL => "-";

my %CONSEQUENCE_NULL = map {$_, 1} (
				    "intergenic_variant",
				    "downstream_gene_variant",
				    "upstream_gene_variant",

				    "incomplete_terminal_codon_variant",
				    # i.e. damaged/imcomplete transcript
				    # annotation
				   );

my $tmp_rank = 1;
my %CONSEQUENCE_RANK = (
			"frameshift_variant" => $tmp_rank++,
			# if both frameshift and stopgain, call frameshift
			# rather than nonsense

			"stop_gained" => $tmp_rank++,

			"splice_acceptor_variant" => $tmp_rank++,
			"splice_donor_variant" => $tmp_rank++,
			"inframe_insertion" => $tmp_rank++,
			"inframe_deletion" => $tmp_rank++,

			"start_lost" => $tmp_rank++,
			"stop_lost" => $tmp_rank++,

			"missense_variant" => $tmp_rank++,
			"synonymous_variant" => $tmp_rank++,
			# prioritize synonymous over splice_region
			# because we have our own custom splice logic,
			# some of these may be promoted to splice

			"splice_region_variant" => $tmp_rank++,
			# may beb promoted to splice

			"coding_sequence_variant" => $tmp_rank++,
			# noncoding only??
			"protein_altering_variant" => $tmp_rank++,
			# insertion that alters protein rather
			# than a clean proteinIns

			"non_coding_transcript_exon_variant" => $tmp_rank++,
			"non_coding_transcript_variant" => $tmp_rank++,

			"stop_retained_variant" => $tmp_rank++,
			# silent stop
			"start_retained_variant" => $tmp_rank++,
			# silent start?

			"intron_variant" => $tmp_rank++,
			# want higher than UTR because raw intronic
			# variants may be promoted to splice or splice_region

			"5_prime_UTR_variant" => $tmp_rank++,
			"3_prime_UTR_variant" => $tmp_rank++,
			# want lower than non-coding so that non-coding
			# gets preference in annotation
		       );

foreach my $csq (sort keys %CONSEQUENCE_NULL) {
  $CONSEQUENCE_RANK{$csq} = $tmp_rank++;
}

my %CONSEQUENCE_TO_SJ = (
			 "5_prime_UTR_variant" => "UTR_5",
			 "3_prime_UTR_variant" => "UTR_3",
			 "coding_sequence_variant" => "exon",
			 "frameshift_variant" => "frameshift",
			 "inframe_deletion" => "proteinDel",
			 "inframe_insertion" => "proteinIns",
			 "intron_variant" => "intron",
			 "missense_variant" => "missense",
			 "splice_acceptor_variant" => "splice",
			 "splice_donor_variant" => "splice",
			 "splice_region_variant" => "splice_region",
			 "start_lost" => "missense",
			 # SJ: chr1.21616906.A.C
			 "stop_gained" => "nonsense",

			 "stop_lost" => "nonsense",
			 # maybe frameshift?

			 "stop_retained_variant" => "silent",
			 "start_retained_variant" => "silent",
			 "synonymous_variant" => "silent",
			 "non_coding_transcript_exon_variant" => "EXON",
			 "non_coding_transcript_variant" => "INTRON",
			 # TO DO: needs more testing
			);

my %SJ_CLASS_RANK;
$tmp_rank = 1;
foreach my $c (qw(
		   frameshift
		   nonsense
		   splice
		   proteinIns
		   proteinDel
		   splice_region
		   missense
		   silent
		   intron
		   exon
		   UTR_3
		   UTR_5
		   EXON
		   INTRON
		)) {
  $SJ_CLASS_RANK{$c} = $tmp_rank++;
}
$SJ_CLASS_RANK{""} = $tmp_rank++;
# for any variant that doesn't yield a class, e.g. intergenic

my $F_PMID = "vep_pubmed";

my %FLAGS;

my $LBC;
my $LBC_LOWER_PRIORITY_FOR_READTHROUGHS = 1;
my $LIMIT_NULL_CLASS;

my @clopts = (
	      #
	      #  input variants:
	      #
	      "-file=s",

	      "-genome=s",
	      # if set, retrieve -fasta and -sjpi from specified genome

	      "-cache=s",
	      # VEP cache directory
	      "-fasta=s",
	      # genome FASTA (single .fai-indexed file)
	      "-sjpi=s",
	      # SJ preferred isoforms file
	      "-refflat=s",

	      #
	      #  options:
	      #
	      "-filter-sjpi=i",

	      "-save-intermediate",
	      "-prep-only=s",

	      "-diff=s",
	      "-find-multi=s",
	      # find variants with > 1 output row for the same gene
	      # (possibly invalid preferred isoforms)

	      "-cache-setup=s",
	      "-no-unpack",
	      "-no-tabix-conversion",

	      "-hack",
	      "-prune-tabix-cache=s",

	      "-codon-annotation-upstream=i" => \$CODON_ANNOTATION_UPSTREAM,
	      "-verbose=i" => \$VERBOSE,
	      "-out=s",

	      "-debug-ram",

	      VariantParsePolicy::get_command_line_parameters(),

	      "-vep-buffer-size=i",
	      "-vep-fork=i",
	      "-vep-only",

	      "-vep-prebuilt=s",
	      "-dump-modules",

	      "-pubmed-only",
	      "-f-pubmed=s" => \$F_PMID,

	      "-limit-by-consequence" => \$LBC,
	      "-limit-null-class" => \$LIMIT_NULL_CLASS,

	      "-extract-refseq-ids",
	      "-check-concordance=s",
	      "-patch-sjpi",
	      "-refseqs=s",

	      "-vep-extra-params=s",
	      "-refflat-cache-limit=i" => \$REFFLAT_CACHE_LIMIT,
	     );

GetOptions(
	   \%FLAGS,
	   @clopts
	  );

if ($FLAGS{"cache-setup"}) {
  cache_setup();
  exit(0);
} elsif ($FLAGS{diff}) {
  result_diff();
  exit(0);
} elsif ($FLAGS{"find-multi"}) {
  find_multi_transcript();
  exit(0);
} elsif ($FLAGS{"prune-tabix-cache"}) {
  prune_tabix_cache();
  exit(0);
} elsif ($FLAGS{hack}) {
  hack();
  exit(0);
} elsif ($FLAGS{"dump-modules"}) {
  print Dumper \%INC;
  exit(0);
} elsif ($FLAGS{"check-concordance"}) {
  check_concordance();
  exit(0);
} elsif ($FLAGS{"patch-sjpi"}) {
  patch_sjpi();
  exit(0);
}

my @h_extra;
printf STDERR "codon annotation policy: %s\n",
  $CODON_ANNOTATION_UPSTREAM ? "always upstream" : "nearest";

$FLAGS{"filter-sjpi"} = 1 unless defined $FLAGS{"filter-sjpi"};
my $SJPI_MODE = $FLAGS{"filter-sjpi"};

printf STDERR "filtering to SJ-preferred isoforms only?: %s\n", $SJPI_MODE ? "yes" : "no";

my $VPP = new VariantParsePolicy("-flags" => \%FLAGS);
#my $v = $VPP->get_variant_from_row($row);

my $genome = $FLAGS{genome} || die "-genome";

my $vrd = config_or_manual(
			      "-config-type" => "genome",
			      "-config-name" => $genome,
			      "-parameter" => "VEP_RESOURCES_DIR",
			      "-manual" => $FLAGS{cache}
			     );
my $cache_dir = $vrd || die "specify -cache CACHE_DIR or -genome GENOME\n";

if ($FLAGS{"extract-refseq-ids"}) {
  my ($main) = glob($cache_dir . "/*refseq/[0-9]*");
  my $outfile = sprintf 'refseqs_%s.txt', basename($main);

  my %refseqs;

  my $callback = sub {
    if (/^\d+\-\d+\.gz$/) {
      my $fn = $_;
      local $/;
      my $fh = universal_open($_);
      my $blob = <$fh>;
      my %hits = map {$_, 1} $blob =~ /(NM_\d+\.\d+)/g;
      foreach (keys %hits) {
	$refseqs{$_} = 1;
      }
      printf STDERR "%s: total:%d here:%s\n", $fn, scalar(keys %refseqs), join ",", sort keys %hits;
    }
  };

  find($callback, $main);
  write_simple_file([ sort keys %refseqs ], $outfile);

  exit(0);
}


#
#  reference FASTA setup:
#
my $FASTA;
unless ($FASTA = $FLAGS{fasta}) {
  if ($vrd) {
    ($FASTA) = glob($vrd . "/*fa");
    die "no fasta file found in $vrd" unless $FASTA and -s $FASTA;
  } else {
    die "specify -fasta FILE or -genome GENOME\n";
  }
#} elsif ($config_genome) {
#  $fasta = $config_genome->{FASTA};
# *** DISABLE: code wants to write an index file, will NOT work with
# TARTAN output!
}
#die "specify -fasta FILE or -genome GENOME\n" unless $fasta;
die "specify -fasta FILE (DO NOT specify a tartan version as VEP wants to write an index file)\n" unless $FASTA;

my $sjpi = get_sjpi();

my $vep = new VEP(
		  "-cache_dir" => $cache_dir,
		  "-fasta" => $FASTA,
		  "-prep_only" => $FLAGS{"prep-only"},
		  "-streaming_mode" => 1,
		 );
$vep->buffer_size($FLAGS{"vep-buffer-size"}) if $FLAGS{"vep-buffer-size"};
$vep->fork_count($FLAGS{"vep-fork"}) if $FLAGS{"vep-fork"};
my $vep_only = $FLAGS{"vep-only"};
$vep->vep_only(1) if $vep_only;
if (my $of = $FLAGS{"vep-prebuilt"}) {
  $vep->prebuilt($of);
}
$vep->vep_extra_params($FLAGS{"vep-extra-params"}) if $FLAGS{"vep-extra-params"};

my $infile = $FLAGS{file} || die "-file";

my $f_refflat = config_or_manual(
				 "-config-type" => "genome",
				 "-config-name" => $genome,
				 "-parameter" => "REFSEQ_REFFLAT",
				 "-manual" => $FLAGS{refflat},
				);
my $AAP = new AAParser();

my $outfile = $FLAGS{out} || basename($infile) . ".hgvs.tab";

my $df = new DelimitedFile("-file" => $infile,
			   "-headers" => 1,
			  );
# while (my $row = $df->next("-ref" => 1)) {  # headerless
ram_debug("VEP prep start");
while (my $row = $df->get_hash()) {
  my $v = $VPP->get_variant_from_row($row);
  $vep->add_variant($v);
}
ram_debug("VEP start");
$vep->run_vep();
ram_debug("VEP end");

if ($vep_only) {
  copy($vep->vep_out, $outfile . ".vep_out.txt") || die;
  exit(0);
} elsif ($FLAGS{"save-intermediate"}) {
  copy($vep->vep_in, $outfile . ".vep_in.txt") || die;
  copy($vep->vep_out, $outfile . ".vep_out.txt") || die;
  printf STDERR "command line: %s\n", $vep->vep_command();
} elsif (my $of = $FLAGS{"prep-only"}) {
  copy($vep->vep_in, $of) || die;
  exit(0);
}

printf STDERR "refflat: %s\n", $f_refflat;
ram_debug("refflat start");
my $rff = new RefFlatFile();
$rff->parse_file(
		 "-refflat" => $f_refflat,
		 "-type" => "refgene",
		);
$rff->cache_base_translations(1);
$rff->cache_limit($REFFLAT_CACHE_LIMIT);
ram_debug("refflat end");

foreach my $r (@{$rff->rows}) {
  $r->{name2} = $r->{bin};
  # hack for GeneAnnotation.pm
}

ram_debug("GA start");
my $GA = new GeneAnnotation(
			    "-style" => "refgene_flatfile",
			    "-refgene_flatfile" => $rff->rows(),
			    "-ignore_non_coding" => 0,
			   );
ram_debug("GA end");

$df = new DelimitedFile("-file" => $infile,
			"-headers" => 1,
		       );

my %h_ignore = map {$_, 1} (
			    "Uploaded_variation",
			    "Location",
			    "Allele",
			   );

my $PUBMED_ONLY = $FLAGS{"pubmed-only"};
my $rpt;

my %missing_transcript_warn;

ram_debug("merge start");
my $headers_initialized;
while (my $row = $df->get_hash()) {
  my $v = $VPP->get_variant_from_row($row);
  my $parsed = get_result_rows($vep, $v);
  # get results, filtered for SJPI
  if ($PUBMED_ONLY and @{$parsed} > 1) {
    # collapse results to just one row to guarantee we don't produce
    # multiple output rows per input row.
    my %pmid;
    foreach my $r (@{$parsed}) {
      if (my $pmid = $r->{Extra_hash}->{PUBMED}) {
	$pmid{$pmid} = 1;
      }
    }
    $parsed->[0]->{Extra_hash}{PUBMED} = join ",", sort keys %pmid;
    $parsed = [ $parsed->[0] ];
  }

  unless ($headers_initialized) {
    # header setup
    if ($PUBMED_ONLY) {
      if (grep {$_ eq $F_PMID} @{$df->headers_raw}) {
	printf STDERR "field %s detected in input file, keeping existing column\n", $F_PMID;
	#    die "Exists";
      } else {
	push @h_extra, $F_PMID;
      }
    } else {
      die "no VEP headers yet??" unless @{$vep->vep_headers};
      @h_extra = grep {!$h_ignore{$_}} @{$vep->vep_headers};
      # add all VEP headers, with a few exceptions

      push @h_extra, "vep_HGVSp";
      push @h_extra, "vep_HGVSc";
      push @h_extra, "vep_sj_gene";
      #  push @h_extra, "vep_sj_gi";
      # NCBI dropped gi numbers on 6/16/2016
      push @h_extra, "vep_sj_class";
      push @h_extra, "vep_sj_aachange";
      push @h_extra, "vep_sj_cdna";
      #push @h_extra, "vep_sj_exception";
      push @h_extra, $F_PMID;
      push @h_extra, "vep_sj_note";
      push @h_extra, "vep_sj_filter_isoform";
      push @h_extra, "vep_sj_filter_isoform_unversioned";
      push @h_extra, "vep_sj_filter_isoform_preferred";
      push @h_extra, "vep_result_count";
    }

    $headers_initialized = 1;
  }


#  die scalar @{$parsed} if @{$parsed} > 1;
  # several reasons this might happen:
  #  - multiple genes overlap.  Here multiple output rows are okay,
  #    downstream tools (e.g. annovar2medals) can drop undesirable rows
  #  - preferred isoform filtering failed, e.g. if SJ-preferred isoform
  #    has been withdrawn by NCBI.  How to deal with this?  Drop all
  #    but one record, or leave these in too?

  my %raw = %{$row};
  my @out;

  my $F_VEP_ROW = "__vep_parsed_row";
  my $F_VEP_CONSEQUENCE_RANK = "__vep_csq_rank";

  if (@{$parsed} == 0) {
    # VEP didn't return a result for this variant.  This seems to
    # happen sometimes esp. with non-canonical chromosomes.
    # Maybe pre-filter these out as in Annovar+?
    my %r = %raw;
    die "no VEP headers yet?" unless @h_extra;
    foreach (@h_extra) {
      $r{$_} = "";
    }
    $r{vep_sj_note} = "no_VEP_results";
    push @out, \%r;
  }

  foreach my $vr (@{$parsed}) {
    dump_die($vr, "new VR", 1) if $VERBOSE;

    my %r = (%raw, %{$vr});
    $r{$F_VEP_ROW} = $vr;

#    my @exceptions;
    my @notes;
    push @notes, @{$vr->{note_queue}} if $vr->{note_queue};

    my $extra = $vr->{Extra_hash};

    #
    # HGVS:
    #
    my $hgvs_c = $extra->{HGVSc} || "";
    my $hgvs_p = $extra->{HGVSp} || "";

    $hgvs_p =~ s/%3D$/=/;
    # unescape "=" character, e.g. p.Ser286%3D
    # http://lists.ensembl.org/pipermail/dev/2015-April/011024.html

    $r{vep_HGVSp} = $hgvs_p;
    $r{vep_HGVSc} = $hgvs_c;

    #
    #  coding sequence and isoform:
    #
    my $cds_change = "";
    my $filter_isoform = "";
    my $filter_isoform_unv = "";
    my $is_preferred = "";
    my $is_noncoding;

    if ($hgvs_c) {
      my @f = split /:/, $hgvs_c;
      ($filter_isoform, $cds_change) = @f;
    } else {
      $filter_isoform = $vr->{Feature};
    }
    if ($filter_isoform) {
      $filter_isoform_unv = $filter_isoform;
      $filter_isoform_unv =~ s/\.\d+$//;
      my $status = $sjpi->is_preferred_nm($filter_isoform_unv);
      $is_preferred = $status ? 1 : 0 if defined $status;
      $is_noncoding = $filter_isoform_unv =~ /^NR_/;
    }
    $r{vep_sj_cdna} = $cds_change;
    $r{vep_sj_filter_isoform} = $filter_isoform;
    $r{vep_sj_filter_isoform_unversioned} = $filter_isoform_unv;
    $r{vep_sj_filter_isoform_preferred} = $is_preferred;
    $r{$F_PMID} = $extra->{PUBMED} || "";

    #
    #  gene symbol:
    #
    my $sj_gene = "";
    if ($filter_isoform_unv) {
      if (my $hits = $rff->find_by_accession($filter_isoform_unv)) {
	$sj_gene = $hits->[0]->{bin} || dump_die($hits->[0], "no bin");
      } else {
	my $warn = $missing_transcript_warn{$filter_isoform_unv} ? 0 : 1;
	$warn = 0 unless $filter_isoform_unv =~ /^NM_/;
	# only main refGene records are actionable

	if ($warn) {
	  # non-coding accessions won't be in refflat
	  printf STDERR "ERROR: can't find refflat entry for %s\n", $filter_isoform_unv;
	  $missing_transcript_warn{$filter_isoform_unv} = 1;
	}
	my $g = $sjpi->get_nm_gene($filter_isoform_unv);
	$sj_gene = $g if $g;
	# rescue from preferred isoforms file
      }
    }
    $r{vep_sj_gene} = $sj_gene;

    #
    #  AAChange:
    #
    my $aa_change = "";
    if ($hgvs_p) {
      $aa_change = parse_hgvs_p($hgvs_p);
      if ($aa_change =~ /^[A-Z][a-z][a-z]\d+delinsTer$/ or
	  # pi_stop.tab: p.Tyr225delinsTer
	  $aa_change =~ /^Ter(\d+)([A-Z][a-z][a-z])extTer\?$/
	  # stop_lost.tab: p.Ter303LeuextTer?
	 ) {
	($aa_change) = derive_aa("-row" => $row, "-vr" => $vr);
      }
    } else {
      ($aa_change) = derive_aa("-row" => $row, "-vr" => $vr);
      push @notes, "no_hgvsp_derived_aa" if $aa_change;
      # TO DO: generate "delins" type HGVSp??
    }

    my $rfr;
    # refFlat record related to this position
    if ($filter_isoform_unv) {
      # accession for this row
      my $raw_hits = find_refflat_records($row);
      my @hits = grep {$_->{name} eq $filter_isoform_unv} @{$raw_hits};
      if (@hits) {
	unless (@hits == 1) {
	  # e.g. X.49208469.T.G for newer refFlat,
	  # multiple mappings for NM_001098407 / GAGE2D
	  printf STDERR "WARNING: multiple refFlat hits for %s!:\n", $filter_isoform_unv;
	  foreach my $h (@hits) {
	    dump_die($h, "debug", 1);
	  }
	}
	$rfr = $hits[0];
      }
    }

    #
    #  class:
    #
    my $class = "";
    my $best;

    if (my $csq = $vr->{Consequence}) {
      my @csq = split /,/, $csq;

      foreach my $csq (@csq) {
	dump_die($vr, "consequence rank not defined for $csq") unless exists $CONSEQUENCE_RANK{$csq};
      }

      my @by_rank = sort {$CONSEQUENCE_RANK{$a} <=>
			    $CONSEQUENCE_RANK{$b}} @csq;

      $best = $by_rank[0];

      if (exists $CONSEQUENCE_TO_SJ{$best}) {
	$class = $CONSEQUENCE_TO_SJ{$best};
	if ($aa_change) {
	  # if AA change provided, shorten protein codes
	  $aa_change = tighten_aa($aa_change);
	} else {
	  ($aa_change) = derive_aa("-row" => $row,
				   "-vr" => $vr,
				   "-class" => $class,
				   "-rfr" => $rfr,
				  );
	}
	if ($class eq "frameshift" and ($aa_change || "") =~ /\*$/) {
	  $aa_change =~ s/\*$/fs/;
	  # preserve SJ policy of formatting as frameshift even if AA
	  # change codes directly to a stop
	}
      } elsif ($best eq "protein_altering_variant") {
	# protein is altered then extended/shortened;
	# need to examine details to class
	($aa_change, $class) = derive_aa("-row" => $row, "-vr" => $vr);
      } elsif ($CONSEQUENCE_NULL{$best}) {
	$aa_change = $class = "";
	# no value in SJ post output for these cases
      } else {
	dump_die($row, "unhandled consequence for $best");
      }

      #
      #  frameshift-related logic:
      #
      if ($class eq "frameshift") {
	if ($aa_change =~ /^\*\d+\*$/) {
	  # frameshift in stop that retains stop
	  $class = "silent";
	} elsif ($aa_change =~ /\*$/) {
	  # some frameshifts code directly to stop, e.g.
	  # RP1     chr8    55537979        -       T       frameshift      N513fs  NM_006269
#	  $class = "nonsense";
	  # 3/2019: leave as frameshift
	} else {
	  ($aa_change) = derive_aa("-row" => $row, "-vr" => $vr,
				   "-class" => $class,
				   "-aa" => $aa_change,
				   "-notes" => \@notes
				  );
	}
      } elsif ($class eq "nonsense") {
	if ($aa_change =~ /^([A-Z]+\d+).*[A-Z]fs\*\d+$/) {
	  # has both a frameshift and a nonsense consequence,
	  # chr1    150315791       -       GTATACTAATATCTCTGCCTGACAGTTGACA
	  # N430Sfs*10: inserted sequence causes a stop, but not at the
	  # beginning
	  $class = "frameshift";
	  $aa_change = sprintf "%sfs", $1;
	} elsif ($aa_change =~ /^(\*\d+[A-Z])ext\*\d+$/) {
	  # *775Rext*14 => *775R
	  $aa_change = $1;
	}
      }
    } else {
      #	dump_die({ %{$row}, %{$vr} }, "no Consequence in output");
      # sometimes blank, e.g. chr1.915389.T.A apparently in
      # some contexts only, e.g. 100.tab
      printf STDERR "WTF: blank Consequence field\n";
    }

    $r{$F_VEP_CONSEQUENCE_RANK} = $best ? ($CONSEQUENCE_RANK{$best} || die "no rank for $best") : $tmp_rank + 100;

    my $pos = get_pos_from_row($row);
    my $soft_promotion;
    my $splice_exon_pos;
    my $splice_exon_pos_upstream;

    if ($rfr) {
      my $strand = $rfr->{strand} || die;

      if ($class eq "silent" or $class eq "missense") {
	#
	#  promote silent/missense to splice at exon edges:
	#
	my $info = $rff->get_base_translations(
					       "-row" => $rfr,
					       "-generate-codons" => 1,
					       "-fasta" => $FASTA,
					       "-exons-only" => 1,
					       "-hash" => 1
					      );
	if (my $site = $info->{$pos}) {
	  if ($site->{is_split_codon_exon_boundary}) {
	    # full promotion with traditional AA call
	    push @notes, sprintf "promoted_%s_to_splice", $class;
	    $class = "splice";
	    $splice_exon_pos = $pos;
	  } elsif ($site->{is_splice_edge}) {
	    # soft promotion, leave AAChange alone (e.g. T125T)
	    push @notes, sprintf "promoted_%s_to_soft_splice", $class;
	    $class = "splice";
	    $splice_exon_pos = $pos;
	    $soft_promotion = 1;
	  }
	}
      } elsif ($class eq "EXON") {
	#
	# generate non-coding exon AA annotation
	#
	my ($ftype, $fno) = $rff->get_annotation_for_position("-row" => $rfr,
							      "-base" => $pos,
							      "-extended" => 1,
							      #								"-intergenic-ok" => 1,
							      #								"-noncoding-synthesize-coding" => 1
							     );
	$aa_change = sprintf 'E%d_exon', $fno;
      } elsif ($class eq "intron" or
	       $class eq "splice" or
	       $class eq "splice_region") {
	# - attempt to promote various regions to splice/splice_region
	# - find base of nearest exon for later annotation step
	my $max = $class eq "splice" ? $SPLICE_MAX_BP_FROM_EXON : $SPLICE_REGION_MAX_BP_FROM_EXON;
	my $ref_seq = get_reference_from_row($row) || "-";
	# in raw bambino output for insertions, reference is blank

	my @accum;
	push @accum, [ $pos, $pos - $max, -1 ];
	push @accum, [ $pos, $pos + ((length($ref_seq) - 1)) + $max, 1 ];
	# in case of a MNV or deletion, look further to make sure we
	# seek properly from both the start and the end of the event.
	# (A) false negative caused by large deletion:
	#     demo_cosmic_jz_splice.tab
	# (B) false positive caused by earlier naive application of
	#     event length in both directions:
	#     demo_large_deletion_not_splice_region.tab

	foreach my $ref (@accum) {
	  my ($pos, $end, $direction) = @{$ref};
	  printf STDERR "start:$pos end:$end direction:$direction\n" if $VERBOSE;
	  my $start_pos = $pos;
	  while (1) {
	    my $ft2 = $rff->get_annotation_for_position(
							"-row" => $rfr,
							"-base" => $pos,
							"-intergenic-ok" => 1,
							"-noncoding-synthesize-coding" => 1,
							#					"-extended" => 1,
						       );
	    printf STDERR "annotation for %s: %s\n", $pos, $ft2 if $VERBOSE;
	    if ($ft2 eq CICERO_CODING) {
	      #		  or ($is_noncoding and $ft2 eq CICERO_UTR_3)
	      # non-coding transcripts will only have 3' UTR.
	      # Using the UTR_3 annotation is problematic because
	      # splice_region will be promoted to splice!
	      #
	      # found exonic sequencing within target search interval
	      #
	      $splice_exon_pos = $pos;
	      if ($class eq "intron") {
		# promote intron to splice region
		$class = "splice_region";
		push @notes, "promoted_to_splice_region";
	      }
	      my $distance = abs($start_pos - $pos);
	      $distance -= (length($ref_seq) - 1) if $direction == 1;
	      # for MNV/deletion, calculate distance from END
	      # of variant rather than start when searching forward

	      if ($class ne "splice") {
		if ($distance <= $SPLICE_PROMOTION_MAX_BP_ACCEPTOR and
		    ($rfr->{strand} eq "-" ?
		     $direction < 0 : $direction > 0)
		   ) {
		  # acceptor side splice:
		  # if +, seeking ahead to acceptor boundary,
		  # if -, seeking back
		  $class = "splice";
		  # promote splice_region to splice
		  my $note = sprintf "promoted_to_acceptor_splice";
		  $note .= "_multi" if length($ref_seq) > 1;
		  $note .= sprintf "_%s", $rfr->{strand};
		  push @notes, $note;
		}

		if ($distance <= $SPLICE_PROMOTION_MAX_BP_DONOR and
		    ($rfr->{strand} eq "+" ?
		     $direction < 0 : $direction > 0)
		   ) {
		  # donor side splice: not usually needed
		  # if +, seeking back to donor boundary
		  # if -, seeking ahead
		  $class = "splice";
		  my $note = sprintf "promoted_to_donor_splice";
		  $note .= "_multi" if length($ref_seq) > 1;
		  $note .= sprintf "_%s", $rfr->{strand};
		  push @notes, $note;
		}
	      }

	      last;
	    }

	    last if $pos == $end;
	    # processed last entry (+ or - direction)
	    $pos += $direction;
	  }
	  last if $splice_exon_pos;
	}
      }

      if ($CODON_ANNOTATION_UPSTREAM and $splice_exon_pos) {
	# found a nearby exon which will be used for exon numbering,
	# however with this option enabled we want to report the exon #
	# UPSTREAM, which may be far away
	my $dir;
	my $spos;
	if ($strand eq "+") {
	  # to search upstream in a + transcript, move - starting
	  # from variant start position
	  $dir = -1;
	  $spos = $pos;
	} elsif ($strand eq "-") {
	  # to search upstream in a - transcript, move + starting
	  # from END of variant start position (e.g. if multi-base)
	  $dir = 1;
	  my $ref_seq = get_reference_from_row($row) || "-";
	  $spos = $pos + (length($ref_seq) - 1);
	}

	for (my $try = 0; $try < 100000; $try++) {
	  my ($ft2) = $rff->get_annotation_for_position(
							"-row" => $rfr,
							"-base" => $spos,
							"-intergenic-ok" => 1,
							"-noncoding-synthesize-coding" => 1,
							#					"-extended" => 1,
						       );
	  if ($ft2 eq CICERO_CODING) {
	    $splice_exon_pos_upstream = $spos;
	    last;
	  } elsif ($ft2 eq CICERO_INTERGENIC) {
	    # moved outside of gene model: give up
	    last;
	  }
	  $spos += $dir;
	}
      }
    }

    if (($class eq "splice" or $class eq "splice_region") and
	not($soft_promotion)) {
      if ($splice_exon_pos) {
	my $info = $rff->get_base_translations(
					       "-row" => $rfr,
					       "-generate-codons" => 1,
					       "-fasta" => $FASTA,
					       "-exons-only" => 1,
					       # TO DO: this doesn't work
					       # properly for non-coding
					       "-hash" => 1
					      );
	if (my $site = $info->{$splice_exon_pos}) {
	  my $site_up;
	  if ($CODON_ANNOTATION_UPSTREAM and $splice_exon_pos_upstream) {
	    $site_up = $info->{$splice_exon_pos_upstream};
	  }

	  if ($site_up) {
	    # report the nearest exon number and the nearest *upstream* codon
	    $aa_change = sprintf '%s%d_E%d%s',
	      $site_up->{codon_code},
		$site_up->{codon_number},
		  $site->{feature_number},
		    $class;
	  } else {
	    # report the nearest exon number and codon
	    $aa_change = sprintf '%s%d_E%d%s',
	      $site->{codon_code},
		$site->{codon_number},
		  $site->{feature_number},
		    $class;
	  }
	} else {
	  dump_die($row, "no base translation info at $splice_exon_pos", 1)
	    unless $is_noncoding;
	}
      } elsif ($VERBOSE) {
	dump_die($row, "no splice exon pos for $class", 1)
	  unless $is_noncoding;
	# FIX ME: this will happen for 3' UTR in non-coding
      }
    }

    if ($SANITY_CHECK_AA_WITH_REFFLAT and
	$rfr and $vr and $vr->{Protein_position}) {
      # there are some cases particularly in GRCh37 where VEP reports
      # incorrect information due to refSeq import problems.
      # attempt to flag these cases by comparing protein positions
      # with those derived from refFlat.
      # Might legitimately happen if:
      #  - local refFlat file uses different transcript version than VEP
      #  - VEP adjusts protein coordinates to e.g. make a "dup" annotation
      #  - VEP/HGVS standardization policy, e.g. left-align
      my $info = $rff->get_base_translations(
					     "-row" => $rfr,
					     "-generate-codons" => 1,
					     "-fasta" => $FASTA,
					     "-exons-only" => 1,
					     # TO DO: this doesn't work
					     # properly for non-coding
					     "-hash" => 1
					    );
      my ($ppos_start, $ppos_end) = split /\-/, ($vr->{Protein_position} || die);

      foreach ($ppos_start, $ppos_end) {
	$_ = "" if $_ and $_ eq "?";
      }

      if (my $stuff = $info->{$pos}) {
	my $rf_codon = $stuff->{codon_number};
	if ($ppos_start and $rf_codon) {
	  my @diffs;
	  push @diffs, abs($ppos_start - $rf_codon);
	  push @diffs, abs($ppos_end - $rf_codon) if $ppos_end;
	  my $diff = min(@diffs);

	  if ($diff >= $SANITY_CHECK_AA_MIN_DISTANCE_TO_REPORT) {
	    push @notes, sprintf "refflat_codon_number_discrepancy=%d", $rf_codon;
	    push @notes, sprintf "refflat_codon_number_discrepancy_distance=%d", $diff;
	  }
	}
      }
    }

    $r{vep_sj_class} = $class;

    #
    # final AA change reformatting for SJ legacy compatibility,
    # required by gedi db loader
    # e.g. /home/medmonso/work/mark/2019_03_04_VEP_dup/all/
    #

    # we may want to disable this as new data model will be more tolerant
    # (per mark wilkinson 3/7/2019)
    if ($aa_change =~ /dup$/) {
      #
      # duplicated amino acids
      #
      if ($aa_change =~ /^([A-Z])(\d+)dup$/) {
	# duplication of single AA, e.g. K477dup becomes K477>KK
	$aa_change = sprintf '%s%d>%s', $1, $2, $1 x 2;
      } elsif ($aa_change =~ /^([A-Z])(\d+)_([A-Z])(\d+)dup$/) {
	my ($aa1, $pos1, $aa2, $pos2) = ($1, $2, $3, $4);
	if ($pos1 == $pos2 - 1) {
	  # duplication of consecutive 2 AAs, e.g. H195_P196dup
	  # becomes P196>PHP report 2nd AA and position, unchanged 2nd
	  # AA, plus new insertion of the 2-AA sequence
	  $aa_change = sprintf '%s%d>%s%s%s', $aa2, $pos2, $aa2, $aa1, $aa2;
	} else {
	  # event spanning more than 2 codons
	  my ($aa_ref, $aa_alt) = split /\//, $vr->{Amino_acids} || "";
	  my $pp = $vr->{Protein_position};
	  if ($aa_ref and $aa_alt and $pp and
	      $aa_ref eq $aa2 and
	      $pos2 == $pp and
	      length($aa_alt) > length($aa_ref) and
	      substr($aa_alt, 0, 1) eq substr($aa_ref,0,1)
	     ) {
	    # e.g. H532_E537dup becomes E537>EHVDSQE (repeat of HVDSQE).
	    # be VERY conservative here as VEP may sometimes shift positions,
	    # potentially catastrophically corrupting reformatting output
	    $aa_change = sprintf '%s%d>%s', $aa2, $pos2, $aa_alt;
	    push @notes, "caution_long_dup_reformat";
	  } else {
	    push @notes, "unhandled_dup_reformat=" . $aa_change;
	  }
	}
      } else {
	# logic will break down for longer insertions.
	# Need to study VEP output closely: can we recover needed info?
	# do we need to refer to refgene2protein?  What happens if
	# VEP shifts/leftaligns inserted sequence during the call,
	# which can happen?
	push @notes, "unhandled_dup_reformat=" . $aa_change;
	# punt and warn
      }
    }

    if ($aa_change =~ /^([A-Z])(\d+)_[A-Z]\d+ins([A-Z]+)$/) {
      # an insertion where the inserted sequence affects 2 codons,
      # and the final result is still in-frame
      my ($aa1, $pos1, $inserted) = ($1, $2, $3);
      $aa_change = sprintf '%s%d>%s%s', $aa1, $pos1, $aa1, $inserted;
      # e.g. V3365_S3366insNYI becomes V3365>VNYI
    }


    $r{vep_sj_aachange} = $aa_change;

    #
    #  exceptions and notes:
    #
#    $r{vep_sj_exception} = join ";", @exceptions;
    $r{vep_sj_note} = join ";", @notes;

#    $rpt->end_row(\%r);
    push @out, \%r;
  }

  if ($SJPI_MODE and @{$parsed} > 1) {
    my %genes = map {$_->{vep_sj_gene} => 1} @out;
    my @genes = keys %genes;
    my $ok;
    $ok = 1 if @genes == @{$parsed};
    # don't warn if the variant simply interacts with multiple genes,
    # as that's a legitimate outcome
    dump_die($row, "SJ-preferred isoform filtering enabled but multiple VEP output rows", 1) unless $ok;
  }

  if ($LIMIT_NULL_CLASS and @out > 1) {
    # if some rows have a generated Class value and some don't, just
    # keep those with an annotation.  This suppresses e.g. upstream/downstream
    # records when a coding annotation exists, and also preserves rare
    # cases where a variant interacts with multiple genes.
    my @has_class;
    my @no_class;

    foreach my $r (@out) {
      if ($r->{vep_sj_class}) {
	push @has_class, $r;
      } else {
	push @no_class, $r;
      }
    }

    if (@has_class and @no_class) {
      @out = @has_class;
    }
  }

  if ($LBC) {
    foreach my $or (@out) {
      my $class = $or->{vep_sj_class};
      dump_die($or, "vep_sj_class not defined") unless defined $class;
      # blank OK
      dump_die($or, "no SJ_CLASS_RANK for $class") unless exists $SJ_CLASS_RANK{$class};
    }

    my @sorted = sort {$SJ_CLASS_RANK{$a->{vep_sj_class}} <=> $SJ_CLASS_RANK{$b->{vep_sj_class}}} @out;
    # sort by final SJ class rank rather than raw VEP rank,
    # in case of promotion (e.g. missense/silent to splice,
    # intron to splice_region)
    my $chosen_class = $sorted[0]->{vep_sj_class};

    if ($LBC_LOWER_PRIORITY_FOR_READTHROUGHS) {
      # e.g. chr1.74715215.T.C has results for both TNNI3K and FPGT-TNNI3K.
      # if the consequences are the same, prioritize the standard gene

      my @filtered = grep {$_->{vep_sj_class} eq $chosen_class} @sorted;
      if (@filtered > 1) {
	# multiple rows with the selected rank
	my %g_std;
	my %g_readthrough;
	foreach my $r (@sorted) {
	  my $g = $r->{vep_sj_gene} || "unknown";
	  if ($g =~ /\-/) {
	    $g_readthrough{$g} = 1;
	  } else {
	    $g_std{$g} = 1;
	  }
	}

	if (scalar keys %g_std == 1) {
	  my ($g) = keys %g_std;
	  @sorted = grep {$_->{vep_sj_gene} eq $g} @sorted;
	}
      }
    }

    @out = $sorted[0];
  }

  unless ($rpt) {
    # delay instantiation since we have to parse first set to get VEP headers
    $rpt = $df->get_reporter(
				"-file" => $outfile,
				"-extra" => \@h_extra,
			       );
    $rpt->auto_qc(1);
  }

  foreach my $r (@out) {
    $r->{vep_result_count} = scalar @out;
    $rpt->end_row($r);
  }



}
$rpt->finish();

ram_debug("done");

sub get_result_rows {
  my ($vep, $v) = @_;

#  my $rows = $vep->get_results($v);

  my $rows = $vep->get_result();
#  dump_die($rows->[0]);
#  die scalar @{$rows};
  unless ($rows) {
    printf STDERR "ERROR: no VEP results returned for %s\n", $v->get_snv4();
    return [];
  }
  my $keys = $rows->[0]->{Uploaded_variation} || die;
  my ($vnum, $snv4) = $vep->parse_tracking_key($keys);

  if ($snv4 eq $v->get_snv4) {
#    printf STDERR "sync OK for %s\n", $snv4;
  } else {
    # VEP does not always report an output row for every input row,
    # breaking synchronization.  This only seems to happen for
    # non-primary references, e.g. "_random" contigs.
    printf STDERR "sync error, infile row is %s, results are for %s\n", $v->get_snv4, $snv4;
    $vep->held_result($rows);
    # hold this VEP result, it should be usable when we reach it
    # in the input file
    return [];
  }

  foreach my $r (@{$rows}) {
    # Extra: IMPACT=HIGH;STRAND=-1;REFSEQ_MATCH=rseq_mrna_match,rseq_ens_match_wt;HGVSc=NM_005646.3:c.2902A>T;HGVSp=NP_005637.3:p.Glu968Ter
#    dump_die($r, "debug", 1);
    my %info;
    foreach my $f (split /;/, $r->{Extra} || die) {
      my @f = split /=/, $f;
      die unless @f == 2;
      my ($k, $v) = @f;
      die if $info{$k};
      $info{$k} = $v;
    }
    $r->{Extra_hash} = \%info;
  }

  if ($SJPI_MODE) {
    #
    #  filter results to SJ preferred isoforms.  If the preferred isoform
    #  is not available, secondary isoforms will be used in the order
    #  specified in the SJPI file.
    #
    my @pref;
    my %gene2rank;
    my @unknown;

    foreach my $r (@{$rows}) {
      if (my $acc = extract_transcript_id($r, 1)) {
	my $rank = $sjpi->get_nm_rank($acc);
	if ($rank) {
	  my $gene = $sjpi->get_nm_gene($acc) || die;
	  # for accessions not in preferred file
	  push @{$gene2rank{$gene}{$rank}}, $r;
	} else {
	  push @unknown, $r;
	}
      } else {
#	dump_die($r, "no HGVSc", 1);
      }
    }

    foreach my $gene (keys %gene2rank) {
      # some variants touch more than one gene: pick the highest-rank
      # isoform in each
      my ($best_rank) = sort {$a <=> $b} keys %{$gene2rank{$gene}};
      my $best_rows = $gene2rank{$gene}{$best_rank};
      unless (scalar @{$best_rows} == 1) {
	# 1.2494330.G.A: hits for both NM_003820.2 and NM_003820.3 (!)

	my %nm2v;
	foreach my $r (@{$best_rows}) {
	  my ($nm, $version) = split /\./, $r->{Feature} || die;
	  $nm2v{$nm}{$version} = $r;
	}
	die "multiple rows in same rank and multi nm" unless scalar keys %nm2v == 1;
	my ($unv) = keys %nm2v;
	my @avail = sort {$b <=> $a} keys %{$nm2v{$unv}};
	my $v_use = $avail[0];
	$best_rows = [ $nm2v{$unv}{$v_use} || die ];
#	printf STDERR "multiple hits for %s, using latest\n", $unv;
	# fairly common occurrence in VEP 91
      }
      push @pref, @{$best_rows};
    }

    if (not(@pref) and @unknown) {
      # in some cases no preferred isoform can be found, may
      # indicate issues with preferred isoform list:
      # - refSeqs withdrawn from NCBI
      # - refSeq IDs newer than in SJPI list
      my ($pr, @suppressed) = @unknown;
      my @notes;
      push @notes, "no_preferred_isoforms_in_results";
      if (@suppressed) {
	push @notes, sprintf "suppressed=%s",
	  join ",", map {extract_transcript_id($_) || "?"} @suppressed;
      }
      $pr->{note_queue} = \@notes;
      push @pref, $pr;
    }

    $rows = \@pref if @pref;
  }

  return $rows;
}

sub get_row_snv4 {
  my ($row) = @_;
  return $VPP->get_variant_from_row($row)->get_snv4();
}

sub get_row_variant {
  my ($row) = @_;
  return $VPP->get_variant_from_row($row);
}

sub cache_setup {
  # should be run in the output_dir of a TARTAn run
  #
  # TO DO:
  # - run tabix conversion step!
  # - prune tabix cache?  Maybe too risky esp. for local use
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $fasta = $config_genome->{FASTA} || die;
  my $archive = abs_path($FLAGS{"cache-setup"} || die);

  my $start_dir = getcwd() || die;
  my $subdir = $FLAGS{out} || "vep_cache";
  mkpath($subdir);
  chdir($subdir) || die "can't cd to $subdir";

  foreach my $thing ($fasta, $fasta . ".fai") {
    die "where is $thing" unless -s $thing;
    my $link = basename($thing);
    die if $thing eq $link;
    unlink $link;
    symlink($thing, $link) || die;
  }

  my $cmd;

  unless ($FLAGS{"no-unpack"}) {
    # debug only
    printf STDERR "unpacking %s...\n", $archive;
    $cmd = sprintf 'tar zxf %s', $archive;
    shell_cmd("-cmd" => $cmd);
  }

  #
  # convert cache to tabix format, critical for good performance:
  #
  unless ($FLAGS{"no-tabix-conversion"}) {
    # debug only
    printf STDERR "running tabix cache conversion...\n";
    $cmd = sprintf 'convert_cache.pl -species all -version all -dir . -force';
    shell_cmd("-cmd" => $cmd);
  }

  # generate single-variant test file:
  # FIX ME: change site depending on 37 or 38
  my $tfw = new TemporaryFileWrangler();
  my $tmp_in = $tfw->get_tempfile("-append" => ".vep_in.txt");
  my $rpt = new Reporter(
			 "-file" => $tmp_in,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   Chr
					   Pos
					   Chr_Allele
					   Alternative_Allele
					)
				      ],
			 "-auto_qc" => 1,
			);
  $rpt->end_row(
		{
		 "Chr" => 17,
		 "Pos" => 7579903,
		 "Chr_Allele" => "",
		 "Alternative_Allele" => "TCA",
		}
	       );
  $rpt->finish();

  if (0) {
    print STDERR "debug copy\n";
    copy($tmp_in, "variant_input.tab") || die;
  }

  $cmd = sprintf '%s -genome %s -file %s -bambino -cache .  -fasta %s', $0, $genome, $tmp_in, basename($fasta);
  shell_cmd("-cmd" => $cmd);
  # perform test run using this cache and FASTA, which will also index
  # the FASTA file. Indexing must be performed before TARTAn run
  # completes and directory becomes read-only.

  unlink(basename($tmp_in) . ".hgvs.tab");

  #
  #  after tabix conversion, prune raw files to save space.
  #  For homo_sapiens_refseq/88_GRCh37, reduces from ~6.7 gb to ~4.2 gb.
  #
  $FLAGS{"prune-tabix-cache"} = ".";
  $FLAGS{force} = 1;
  prune_tabix_cache();

  chdir($start_dir);
  # session can now be ended and run added to index
  # under VEP_RESOURCES_DIR
}

sub get_pos_from_row {
  my ($row) = @_;
  return $row->{$VPP->f_pos || die "no f_pos"};
}

sub get_reference_from_row {
  my ($row) = @_;
  return $row->{$VPP->f_ra || die "no f_ra"};
}

sub find_refflat_records {
  my ($row) = @_;
  my $chr = $row->{$VPP->f_chr || die "no chr"};
  my $pos = get_pos_from_row($row);
  my $results;
  if ($GA->find(
		"-reference" => $chr,
		"-start" => $pos,
		"-end" => $pos
	       )) {
    $results = $GA->results_rows();
  }
  return $results;
}

sub derive_aa {
  my %options = @_;
  my $row = $options{"-row"} || die "-row";
  my $vr = $options{"-vr"} || die "-vr";
  my $assigned_class = $options{"-class"};
  my $assigned_aa = $options{"-aa"};
  my $notes = $options{"-notes"};

  my $aa_change = "";
  my $class = "";
  my $aas = $vr->{Amino_acids} || dump_die($row, "no Amino_acids");
  if ($aas eq VEP_NULL) {
    if ($assigned_class and $assigned_class =~ /UTR/) {
      if (my $rfr = $options{"-rfr"}) {
	my $pos = get_pos_from_row($row);
	my ($ftype, $fno) = $rff->get_annotation_for_position("-row" => $rfr,
							      "-base" => $pos,
							      "-extended" => 1,
							      #								"-intergenic-ok" => 1,
							      #								"-noncoding-synthesize-coding" => 1
							     );
	$aa_change = sprintf 'E%d_%s', $fno, $assigned_class;
      } else {
	# e.g. utr_crash.tab: site does not actually appear to be in UTR!
	$aa_change = sprintf 'E?_%s', $assigned_class;
      }
    }
  } else {
    my ($aa_ref, $aa_var) = split /\//, $aas;
    my ($ppos_start, $ppos_end) = split /\-/, ($vr->{Protein_position} || die);

    if (($assigned_class || "") eq "frameshift") {
      # frameshift.tab: SJ says Q130fs, VEP says p.Thr131AspfsTer24 but
      # event starts at Q130
      # TO DO: more testing (w/deletion to avoid insertion ambiguity)
      if ($aa_ref eq VEP_NULL) {
	# in some cases VEP doesn't provide the amino acids, e.g.
	# chr1    22927220        -       GGCC    frameshift      S819fs
	# is -/GX, even though AA is p.Ser819GlyfsTer83
	die "-aa" unless $assigned_aa;
	die "-notes" unless $notes;
#	if ($assigned_aa =~ /^([A-Z][a-z][a-z])(\d+)/) {
	if ($assigned_aa =~ /^([A-Z\*])(\d+)/) {
	  # now already reformatted to short code
	  $aa_change = sprintf '%s%dfs', $1, $2;
	  push @{$notes}, "frameshift_annotation_rescue";
	} elsif ($assigned_aa =~ /^(\-)(\d+).*X/) {
	  # NRAP    chr10   115405584       -       GCAC    frameshift      E370fs  NM_006175
	  # vep reports -370>VX
	  $aa_change = sprintf '%s%dfs', $UNKNOWN_AA_CHAR, $2;
	  push @{$notes}, "frameshift_annotation_rescue";
	} else {
	  dump_die($row, "ERROR, can't generate fs for $assigned_aa");
	}
      } elsif (not($ppos_end) or $ppos_start == $ppos_end) {
	$aa_change = sprintf '%s%dfs', $aa_ref, $ppos_start;
      } else {
	$aa_change = sprintf '%s%d_%s%dfs',
	  substr($aa_ref, 0, 1), $ppos_start,
	    substr($aa_ref, -1), $ppos_end;
      }
    } elsif ($aas eq "*") {
      # VPS37B  chr12   123351665       -       TC      exon    E4_exon Unknown
      # = frameshift_variant,stop_retained_variant
      $aa_change = sprintf '*%d*', $ppos_start;
      $class = "silent";
    } elsif (length($aa_ref) == 1 and length($aa_var) == 1) {
      # single AA change
      $aa_change = sprintf '%s%d%s', $aa_ref, $ppos_start, $aa_var;
      # to do: set class based on affected AAs?
    } elsif (length($aa_ref) == 1 and length($aa_var) > 1) {
      # final protein is in frame and longer, so call proteinIns
      $aa_change = sprintf '%s%d>%s', $aa_ref, $ppos_start, $aa_var;
      $class = "proteinIns";
    } elsif (length($aa_ref) > 1 and length($aa_var) == 1) {
      $aa_change = sprintf '%s%d_%s%d%s',
	substr($aa_ref, 0, 1), $ppos_start,
	  substr($aa_ref, -1), $ppos_end,
	    ($aa_var eq VEP_NULL ? "del" : ">" . $aa_var);
      $class = "proteinDel";
      # final protein is in frame and shorter, so call proteinDel
    } else {
      # e.g. 3.195511354.GTGTCACCT.CTGCAC
      # complex protein sub: there is no SJ class for this!
      $aa_change = parse_hgvs_p($vr->{Extra_hash}{HGVSp});
      $class = length($aa_ref) > length($aa_var) ? "proteinDel" : "proteinIns";
      # NOT ACCURATE but what can we do?
    }
  }
  return ($aa_change, $class);
}

sub result_diff {
  my $infile = $FLAGS{diff} || die "-diff";
  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  while (my $row = $df->get_hash()) {
    my $nm_old = $row->{mRNA_acc} || dump_die($row, "no acc");
    my $nm_new = $row->{vep_sj_filter_isoform_unversioned} || "";
    # might not be present, e.g. modifier
    next unless $nm_old eq $nm_new;
    # VEP may report results for a different isoform, or multiple isoforms

    my $class_old = $row->{Class} || die;
    my $class_new = $row->{vep_sj_class} || "";
    # no_class.tab: sometimes there is no Consequence field

    my $aa_old = $row->{AAChange} || die;
    my $aa_new = $row->{vep_sj_aachange} || "";

    my $snv4 = join ".", @{$row}{qw(Chr WU_HG19_Pos ReferenceAllele MutantAllele)};

    my @notes;
    my @diffs;

    my $new_splice_logic;
    unless (uc($class_old) eq uc($class_new)) {
      if (($class_old eq "missense" or $class_old eq "silent") and
	  $class_new = "splice") {
	$new_splice_logic = 1;
      } else {
	push @diffs, "class";
      }
    }

    unless ($aa_old eq $aa_new) {
      if ($new_splice_logic) {
	# don't flag as mismatch in this case
	push @notes, "new_splice_logic";
      } else {
	push @diffs, "AA";
      }
    }

    if (my $note = $row->{vep_sj_note}) {
      push @notes, $1 if $note =~ /(codon_number_discrepancy_distance=\d+)/;
    }

    printf "%s\n", join "\t",
      $snv4, $class_old, $class_new, $aa_old, $aa_new,
	$row->{vep_sj_filter_isoform},
	@notes,
	  @diffs ? @diffs : "identical";

  }
}

sub tighten_aa {
  # convert 3-character protein codes to 1-character where possible,
  # and condense delins formatting
  my ($aa_change) = @_;

  my $raw = $aa_change;
  $aa_change =~ s/([A-Z][a-z][a-z])/shorten_aa_code($1)/eg;
  $aa_change =~ s/delins/>/;
  # Y418_R420delinsW becomes Y418_R420>W
  if ($VERBOSE) {
    printf STDERR "reformatted %s => %s\n", $raw, $aa_change
      unless $raw eq $aa_change;
  }

  # if ($aa_change =~ /^([A-Z][a-z][a-z])(\d+)(dup|del)$/) {
  #   # Ser777del => S777del
  #   my ($aa_long, $aa_pos, $event) = ($1, $2, $3);
  #   my $ct = new Bio::Tools::CodonTable();
  #   my @codons = $ct->revtranslate($aa_long);
  #   my $cooked = $ct->translate($codons[0]);
  #   $aa_change = sprintf '%s%d%s', $cooked, $aa_pos, $event;
  # } elsif ($aa_change =~ /^([A-Z][a-z][a-z])(\d+)_([A-Z][a-z][a-z])(\d+)(dup|del)$/) {
  #   # Trp8_Ala12dup
  #   my ($aa_long1, $aa_pos1, $aa_long2, $aa_pos2, $event) = ($1, $2, $3, $4, $5);
  #   my $ct = new Bio::Tools::CodonTable();
  #   my @codons = $ct->revtranslate($aa_long1);
  #   my $cooked1 = $ct->translate($codons[0]);
  #   @codons = $ct->revtranslate($aa_long2);
  #   my $cooked2 = $ct->translate($codons[0]);
  #   $aa_change = sprintf '%s%d_%s%d%s',
  #     $cooked1, $aa_pos1,
  # 	$cooked2, $aa_pos2,
  # 	  $event;
  # }

  return $aa_change;
}

sub shorten_aa_code {
  my ($aa_long) = @_;
  my $ct = new Bio::Tools::CodonTable();
  if (my @codons = $ct->revtranslate($aa_long)) {
    $aa_long = $ct->translate($codons[0]);
  }
  return $aa_long;
}

sub hack {
  my @things;
  push @things, "Glu39dup";
  push @things, "Xyz39del";

  foreach my $str (@things) {
    my $raw = $str;
    $str =~ s/([A-Z][a-z][a-z])/shorten_aa_code($1)/eg;
    printf "%s => %s %d\n", $raw, $str, $raw eq $str ? 0 : 1;
  }
}

sub prune_tabix_cache {
  my $root = $FLAGS{"prune-tabix-cache"};

  unless ($FLAGS{force}) {
    $| = 1;
    printf "about to prune %s of non-tabix files: confirm? [y/n]: ", $root;
    my $resp = <STDIN>;
    chomp $resp;
    unless (uc($resp || "") eq "Y") {
      print "quitting.\n";
      exit(1);
    }
  }

  my %dirs;
  my $callback = sub {
    $dirs{$File::Find::dir} = 1;
  };

  find($callback, $root);

  foreach my $dir (keys %dirs) {
    if (-f $dir . "/all_vars.gz") {
      my @all = glob($dir . "/*.gz");
      my @other = grep {/_var.gz$/} @all;
      if (@other) {
	printf STDERR "unlink %s\n", join ",", @other;
	unlink @other;
      }
    }
  }


}

sub find_multi_transcript {
  my $infile = $FLAGS{"find-multi"};

  $VPP = new VariantParsePolicy("-flags" => \%FLAGS);

  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  my %track;
  while (my $row = $df->get_hash()) {
    my $snv4 = get_row_snv4($row);
    my $gene = $row->{vep_sj_gene} || next;
    push @{$track{$gene}{$snv4}}, $row;
  }

  foreach my $gene (sort keys %track) {
    foreach my $snv4 (keys %{$track{$gene}}) {
      printf "%s: %s %s\n", $gene, $snv4, join(",", map {$_->{vep_sj_filter_isoform}} @{$track{$gene}{$snv4}}) if scalar @{$track{$gene}{$snv4}} > 1;
    }
  }
}

sub extract_transcript_id {
  my ($row, $return_unversioned) = @_;
  my $acc = "";
  if (my $c = $row->{Extra_hash}{HGVSc}) {
    my @f = split /:/, $c;
    $acc = $f[0];
  } else {
    $acc = $row->{Feature};
  }
  $acc =~ s/\.\d+$// if $return_unversioned;
  return $acc;
}

sub parse_hgvs_p {
  my ($hgvs_p) = @_;
  my $aa_change = "";
  if ($hgvs_p) {
    my @f = split /:/, $hgvs_p;
    if (@f == 2) {
      $aa_change = $f[1];
    } elsif ($hgvs_p =~ /(p\.\S+)$/) {
      # NC_000022.10:IGLV5-52:u_t_1.1:p.His117Gln
      # NC_000022.10:IGLV2-8:u_t_1.1:p.Leu6=
      $aa_change = $1;
    } else {
      die "can't parse HGVSp $hgvs_p";
    }
    $aa_change =~ s/^p\.//;
    if (my $cooked = $AAP->parse_substitution($aa_change)) {
      $aa_change = $cooked;
    }
  }
  return $aa_change;
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

    log_message(sprintf "RAM at %s: RSS:%d VSZ:%d", $label, @info{qw(RSS VSZ)});
  }
}

sub check_concordance {
  # compare original SJ-post style output w/VEP+ output
  # - TO DO:
  #   - detect SJPI problems
  #   - detect gene symbol swaps

  my $vpp = new VariantParsePolicy("-flags" => \%FLAGS);

  my $vep_refseqs = get_vep_refseqs();

  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $ger_dir = $config_genome->{GENE_EXON_REGION_DIR} || die;
  my @ger = glob($ger_dir . "/*");
  my %ger_genes;
  foreach my $f (@ger) {
#    print STDERR "$f\n";
    open(IN, $f) || die;
    while (<IN>) {
      my @f = split /\t/, $_;
      $ger_genes{$f[0]}=1;
    }
  }
  my $gsm_ger = new_gsm_lite();
  foreach (sort keys %ger_genes) {
#    printf STDERR "add GER %s\n", $_;
    $gsm_ger->add_gene("-gene" => $_);
  }

  my $infile = $FLAGS{"check-concordance"} || die;
  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );


  my %compare;
  $compare{GeneName} = "vep_sj_gene";
  $compare{Class} = "vep_sj_class";
  $compare{mRNA_acc} = "vep_sj_filter_isoform_unversioned";
  $compare{AAChange} = "vep_sj_aachange";

  my @extra;
  push @extra, "mismatch_field_count";
  foreach (sort keys %compare) {
    push @extra, "mismatch_" . $_;
  }
  push @extra, "concordance_note";

  my $rpt = $df->get_reporter(
			      "-file" => $infile . ".concordance.tab",
			      "-extra" => \@extra,
			     );

  my @f_compare = sort keys %compare;
  my $identical = 0;
  my $not_identical = 0;
  my %gene_rescue;
  my $gene_rescue = 0;
  my %mismatch_pattern;
  my $ignore_invalid_refseq = 0;
  my %invalid_refseq;
  while (my $row = $df->get_hash()) {
    my @notes;

    my $nm = $row->{mRNA_acc} || die;
    unless ($vep_refseqs->{$nm}) {
      $ignore_invalid_refseq++;
      $invalid_refseq{$nm} = 1;
      next;
    }

    my %mismatches;
    foreach my $f (@f_compare) {
#      $mismatches{$f} = 1 unless $row->{$f} eq $row->{$compare{$f}};
      $mismatches{$f} = 1 unless lc($row->{$f}) eq lc($row->{$compare{$f}});
      # lowercase comparison for slightly-altered data,
      # e.g. proteinDel vs. proteindel
    }

    if ($mismatches{GeneName}) {
      # if the gene symbol is different, is this just a symbol change?
      my $old = $row->{GeneName};
      my $new = $row->{vep_sj_gene};

      if ($new and my $mapped = $gsm_ger->find($new)) {
	# gene symbol is different between SJ post/GENE_EXON_REGION and VEP,
	# but it's equivalent
	$gene_rescue{hugo_symbol_change}++;
	delete $mismatches{GeneName};
	# consider resolved
      } else {
	# ?
	if ($row->{mRNA_acc} eq $row->{vep_sj_filter_isoform_unversioned}) {
	  $gene_rescue{same_NM_but_symbol_change}++;
	  delete $mismatches{GeneName};
	  # consider resolved
	} else {
	  printf STDERR "gene rescue failed for %s: orig=%s new=%s\n",
	    $vpp->get_variant_from_row($row)->get_snv4(),
	      $old, $new;
	}
      }
    }

    if ($mismatches{Class} and $mismatches{AAChange}) {
      my $class_old = $row->{Class};
      my $class_new = $row->{vep_sj_class};
      if ($class_old eq "frameshift" and $class_new eq "nonsense") {
	$mismatches{Class} = $mismatches{AAChange} = 2;
	# likely explained
	push @notes, "likely_frameshift_direct_stopgain_policy_change";
      }
    }

    $row->{mismatch_field_count} = scalar keys %mismatches;
    if (%mismatches) {
      $not_identical++;

      my $key = join "_", sort keys %mismatches;
      $mismatch_pattern{$key}++;
    } else {
      $identical++;
    }

    foreach my $f (keys %compare) {
      $row->{"mismatch_" . $f} = $mismatches{$f} || 0;
    }

    $row->{concordance_note} = join ",", @notes;
    $rpt->end_row($row);

  }
  $rpt->finish();

  printf STDERR "ignored because refseq not in VEP cache: %d\n", $ignore_invalid_refseq;
  write_simple_file([sort keys %invalid_refseq], "invalid_refseq.txt");

  printf STDERR "identical:%d mismatches:%d\n", $identical, $not_identical;
  printf STDERR "mismatch patterns:\n";
  foreach my $k (sort {$mismatch_pattern{$b} <=> $mismatch_pattern{$a}} keys %mismatch_pattern) {
    printf STDERR "  %s: %d\n", $k, $mismatch_pattern{$k};
  }

  printf STDERR "gene symbol rescues:\n";
  foreach (sort keys %gene_rescue) {
    printf STDERR "  %s: %d\n", $_, $gene_rescue{$_};
  }

}

sub get_sjpi {
  my $genome = $FLAGS{genome} || die "-genome";
  my $f_sjpi = config_or_manual(
				"-config-type" => "genome",
				"-config-name" => $genome,
				"-parameter" => "GENE_TRANSCRIPT_MATRIX",
				"-manual" => $FLAGS{sjpi}
			       );

  die "specify -sjpi FILE or -genome GENOME\n" unless $f_sjpi;
  return new SJPreferredIsoform("-file" => $f_sjpi);
}

sub patch_sjpi {
  # patch SJ preferred isoforms file for compatibility w/refseqs in VEP cache
  my $sjpi = get_sjpi();

  my $gene2nms = $sjpi->get_gene2nms();

  my $genes = $sjpi->get_genes_ordered();
  # genes as ordered in original file
  my $vep_refseqs = get_vep_refseqs();

  my $count_in_vep = 0;
  my $count_no_vep = 0;
  my $count_no_vep_preferred = 0;
  my %problem_genes;
  my @missing_nm;

  foreach my $gene (@{$genes}) {
    my $entry = 0;
    foreach my $nm (@{$gene2nms->{$gene}}) {
      $entry++;
      if ($vep_refseqs->{$nm}) {
	# OK
	$count_in_vep++;
      } else {
	$count_no_vep++;
	$count_no_vep_preferred++ if $entry == 1;
	printf STDERR "missing $gene $nm $entry\n";
	$problem_genes{$gene} = 1;
	push @missing_nm, $nm;
      }
    }
  }

  printf STDERR "genes:%d problem_genes:%d in_vep:%d not:%d preferred_not:%d\n",
    scalar(@{$genes}),
      scalar(keys %problem_genes),
	$count_in_vep, $count_no_vep, $count_no_vep_preferred;

  my $eu = new EUtilities();
  my $result_files = $eu->fetch_genbank("-ids" => \@missing_nm);
  # fetch old records from NCBI, and cache

  my %missing_to_tv;
  #
  # for missing IDs, find transcript variant # (if any):
  #
  foreach my $rf (@{$result_files}) {
    my $stream = Bio::SeqIO->new("-file" => $rf,
				 -format => 'GenBank');
    while (my $record = $stream->next_seq()) {
      # Bio::SeqIO::genbank
      my $accession = $record->accession_number();
      my $version = $record->version();
      my $desc = $record->desc() || die "no desc for $accession";
      my $transcript_variant = $desc =~ /transcript variant (\d+)/ ? $1 : 0;

      $missing_to_tv{$accession} = $transcript_variant;
    }
  }

  # see if we have an ID for this transcipt variant #, or failing that
  # the gene symbol:


  die scalar keys %missing_to_tv;

}

sub get_vep_refseqs {
  my $f_vep_refseqs = $FLAGS{refseqs} || die "-refseqs";
  open(RS, $f_vep_refseqs) || die;
  my %vep_refseq;
  while (<RS>) {
    chomp;
    s/\.\d+$//;
#    print STDERR "save $_\n";
    $vep_refseq{$_}=1;
  }
  return \%vep_refseq;
}

