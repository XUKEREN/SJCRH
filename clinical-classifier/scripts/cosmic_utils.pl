#!/bin/env perl
# COSMIC db utilities
# MNE 3/2017
#
# TO DO:
# - eliminate JZ blacklisted publications?: may not be needed if starting
#   with cleaned version since these should be gone already?
# - verify all gold/silver genes have entries?: not needed here,
#   genomic matching ultimately required

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use FileUtils qw(write_simple_file read_simple_file find_binary);
use MiscUtils qw(dump_die build_argv_list get_hash_option);
use DelimitedFile;
use Reporter;
use TdtConfig;
use DelimitedFileHP;
use ReviewableReportable;

my $MAX_PUBMED_RECORDS = 5;

my $MIN_COUNT_FOR_RECURRENT = 5;


my @SUSPICIOUS_ABSTRACT_TERMS = (
				 "single-cell",
				 "intratumoral",
				 "intratumour",
				 "intratumor",
				 "intra-tumor",
				 "intraindividual",
				 "intraclonal",
#				 "ITH",
				 "metastases",
				 # "We performed whole-exome sequencing of primary melanomas and multiple matched metastases from eight patients"

				);
# maybe "heterogeneity"? cell lines?

my %FLAGS;
my @clopts = (
	      "-genome=s",
	      "-build-pubmed-summary",
	      "-cosmic-tabix=s",

	      "-hack-compare",

	      "-max-pubmed=i" => \$MAX_PUBMED_RECORDS,

	      "-find-recurrent-silent",
	      "-cosmic=s",

	      "-max=i",

	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

if ($FLAGS{"build-pubmed-summary"}) {
  build_pubmed_summary();
  merge_pubmed_summary();
} elsif ($FLAGS{"hack-compare"}) {
  hack_compare();
} elsif ($FLAGS{"find-recurrent-silent"}) {
  find_recurrent_silent();
}

sub build_pubmed_summary {
  # original version, bucketing by snv4/aa/pmid.
  # more accurate but may over-separate evidence with same AA.
  #
  # analyze (cleaned) COSMIC to identify recurrent variants,
  # used by medal ceremony to report important recurrent variants
  # in gold/silver genes, as well as generally in the Evidence field
  # for all COSMIC matches.
  #
  # medal ceremony needs to filter by minimum validated sample count
  # (3+) which seems to be the level in JZ's gold/silver lists.

  my $tf = $FLAGS{"cosmic-tabix"};
  unless ($tf) {
    my $genome = $FLAGS{genome} || die "-genome";
    my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
    my $tabix_dir = $config_genome->{CLINCLS_COSMIC_TABIX_DIR} || die;
    ($tf) = glob($tabix_dir . "/*.gz");
  }

  open(IN, sprintf "tabix -l %s|", $tf) || die;
  my @chroms;
  while (<IN>) {
    chomp;
    push @chroms, $_;
  }
  die unless @chroms;
  @chroms = sort @chroms;

  my $rpt = new Reporter(
			 "-file" => "cosmic_pubmed_summary.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   Gene
					   AAchange
					   Chr
					   WU_HG19_Pos
					   ReferenceAllele
					   MutantAllele
					   TotalSample
					   TotalVerifiedSample
					   PubmedInfo
					)
				      ],
			 "-auto_qc" => 1,
			);


  foreach my $chrom (@chroms) {
    my $cmd = sprintf 'tabix -h %s %s|', $tf, $chrom;
    open(IN, $cmd) || die;

    my $df = new DelimitedFile("-fh" => \*IN,
			       "-headers" => 1,
			      );

    my %track;
    while (my $row = $df->get_hash()) {
      die unless exists $row->{Pubmed_PMID};

#      dump_die($row);
#      my $sid = $row->{ID_sample} || die;
# http://cancer.sanger.ac.uk/cosmic/download
# "A sample id is used to identify a sample within the COSMIC database. There can be multiple ids, if the same sample has been entered into the database multiple times from different papers."
# so, maybe best to just use the sample name, so that e.g. TCGA-XXXX
# will not be duplicated if submitted by multiple users.
      my $sid = $row->{"Sample name"} || die;
      my $pmid = $row->{Pubmed_PMID} || 0;
      my $chr = $row->{"#Chr"} || die;
      my $pos = $row->{"WU_HG19_Pos"} || die;
      my $ra = $row->{ReferenceAllele} || die;
      my $va = $row->{MutantAllele} || die;
      my $aa = $row->{"Mutation AA"} || "";
      # AA might not be present
      my $snv4 = join ".", $chr, $pos, $ra, $va;
      my $status = $row->{"Mutation somatic status"} || die;
      $track{$snv4}{$aa}{$pmid}{all}{$sid} = 1;
      # bucket by AA as well for easier comparison with JZ's version.
      # SNV4 will be used for linkage to main data though.
      $track{$snv4}{$aa}{$pmid}{validated}{$sid} = 1 if $status eq "Confirmed somatic variant";
      $track{$snv4}{$aa}{$pmid}{demo_row} = $row;
    }

    printf "%s: unique_sites:%d\n", $chrom, scalar keys %track;

    foreach my $snv4 (sort keys %track) {
      foreach my $aa (sort keys %{$track{$snv4}}) {
	my @pmi;
	my %all_samples;
	my %all_verified_samples;
	my %r;
	@r{qw(Chr WU_HG19_Pos ReferenceAllele MutantAllele)} = split /\./, $snv4;
	# PubmedInfo: 18224685;56;56,16900509;35;35,23103869;35;35,18772397;33;33,11290569;32;32

	foreach my $pmid (sort keys %{$track{$snv4}{$aa}}) {
	  $r{AAchange} = $aa;
	  my @total_samples = keys %{$track{$snv4}{$aa}{$pmid}{all}};
	  my @total_verified_samples = keys %{$track{$snv4}{$aa}{$pmid}{validated}};
	  my $total = scalar @total_samples;
	  my $valid = scalar @total_verified_samples;
	  foreach (@total_samples) {
	    $all_samples{$_} = 1;
	  }
	  foreach (@total_verified_samples) {
	    $all_verified_samples{$_} = 1;
	  }

	  my %pmi;
	  $pmi{string} = join ";", $pmid, $total, $valid;
	  $pmi{pmid} = $pmid;
	  $pmi{total} = $total;
	  $pmi{valid} = $valid;
	  push @pmi, \%pmi unless $pmid == 0;

	  my $row = $track{$snv4}{$aa}{$pmid}{demo_row} || die;

	  my $gene = $row->{"Gene name"};
	  $gene =~ s/_loc\w+$//;
	  # trailing refflat "sharp" uniquifier

	  $r{Gene} = $gene;
	}

	my @sorted = sort {$b->{valid} <=> $a->{valid} || $b->{total} <=> $a->{total} || $b->{pmid} <=> $a->{pmid}} @pmi;
	my $idx_end = $#sorted;
	$idx_end = $MAX_PUBMED_RECORDS - 1 if $idx_end >= $MAX_PUBMED_RECORDS;

	$r{PubmedInfo} = join ",", map {$_->{string}} @sorted[0 .. $idx_end];
	$r{TotalSample} = scalar keys %all_samples;
	$r{TotalVerifiedSample} = scalar keys %all_verified_samples;
	$rpt->end_row(\%r);
      }
    }
  }
  $rpt->finish;

  # TO DO: verify that all cosmic gold/silver genes observed

}

sub hack_compare {
  my $df = new DelimitedFile(
			     "-file" => "/research/rgs01/resgen/dev/tartan/runs/tartan_import/E5Gl2YTO/output/reference/Homo_sapiens/GRCh37-lite/CLINICAL/COSMIC/COSMIC_SNV_indel_pubmed_info.txt",
			     "-headers" => 1,
			    );
  my %old;
  while (my $row = $df->get_hash()) {
    my $key = join "_", @{$row}{qw(Gene AAchange)};
    $old{$key} = $row->{PubmedInfo};
  }

  $df = new DelimitedFile(
			  "-file" => "cosmic_pubmed_summary.tab",
			  "-headers" => 1,
			 );
  my %new;
  while (my $row = $df->get_hash()) {
    my $key = join "_", @{$row}{qw(Gene AAchange)};
    $new{$key} = $row->{PubmedInfo};
  }

  foreach my $key (sort keys %old) {
    my $tag;
    if ($new{$key}) {
      $tag = "found";
    } else {
      my $info = $old{$key};
      my @entries = split /,/, $info;
      $tag = sprintf 'missing_%d', scalar @entries;
    }
    printf STDERR "%s: %s\n", $key, $tag;
  }
}

sub merge_pubmed_summary {
  my $infile = "cosmic_pubmed_summary.tab";
  my $outfile = $infile . ".merged.tab";

  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
  			      "-auto_qc" => 1,
			     );

  my %track;
  while (my $row = $df->get_hash()) {
    # to more closely emulate JZ's process, run 2nd pass flattening
    # sites with same gene, AA, and position.
    # - collapse JAK3 M511I, leave BRAF V600E separate.
    # - merge and re-sort COSMIC records: e.g. JAK3 should be 22237106;2;2
    my $key = join "_", @{$row}{qw(Gene AAchange WU_HG19_Pos)};
    push @{$track{$key}}, $row;
  }

  foreach my $key (sort keys %track) {
    my $set = $track{$key};
    my $ts = 0;
    my $tv = 0;
    my %pmid;
    foreach my $row (@{$set}) {
      $ts += $row->{TotalSample};
      $tv += $row->{TotalVerifiedSample};

      my @pmi = split /,/, $row->{PubmedInfo};
      foreach my $pmi (@pmi) {
	my ($pmid, $total, $verified) = split /;/, $pmi;
	$pmid{$pmid}{total} += $total;
	$pmid{$pmid}{valid} += $verified;
      }
    }
    delete $pmid{0};

    my @sorted_pmid = sort {$pmid{$b}{valid} <=> $pmid{$a}{valid} ||
			 $pmid{$b}{total} <=> $pmid{$a}{total} ||
			   $b <=> $a} keys %pmid;
    my @sorted = map {join ";", $_, $pmid{$_}{total}, $pmid{$_}{valid}} @sorted_pmid;
    my $idx_end = $#sorted;
    $idx_end = $MAX_PUBMED_RECORDS - 1 if $idx_end >= $MAX_PUBMED_RECORDS;

    my $pmi_new = join ",", @sorted[0 .. $idx_end];

    foreach my $row (@{$set}) {
      $row->{TotalSample} = $ts;
      $row->{TotalVerifiedSample} = $tv;
      $row->{PubmedInfo} = $pmi_new;
      $rpt->end_row($row);
    }
  }


  $rpt->finish();
}

sub find_recurrent_silent {
  my $f_cosmic = $FLAGS{cosmic} || die "-cosmic";

  find_binary("pmid2abstract.pl", "-die" => 1);

#  my $df = new DelimitedFile("-file" => $f_cosmic,
#			     "-headers" => 1,
#			     );

  my $df = new DelimitedFileHP(
			     "-file" => $f_cosmic,
			    );
  my %silent;
  my %pmid2silent;
  my %study2silent;
  my %study2total;
  my %study2unique;
  my %study2pmid;
  my %study2sample;
  my %study2tumor;
  my %silent2pmid;
  my $rows = 0;
  my $start = time;
  my $max = $FLAGS{max};
#  while (my $row = $df->get_hash()) {

  my %unusable;

  # exclude all records w/o PMID?

  my $total_records = 0;
  my $total_records_with_pmid = 0;

  my %real_pmid;

  while ($df->next_row()) {
#    my $aa = $row->{"Mutation AA"} || dump_die($row, "no AA");
    my $gene = $df->get_value("Gene name") || dump_die($df->get_hash(), "no Gene name");
    my $aa = $df->get_value("Mutation AA") || dump_die($df->get_hash(), "no AA");

    my $gws = $df->get_value("Genome-wide screen") || dump_die($df->get_hash(), "no Genome-wide-screen");

    if ($gws eq "n") {
      $unusable{"non_genome_wide_screen"}++;
      next;
    } elsif ($gws ne "y") {
      die "unknown gws value $gws";
    }

    $total_records++;

#    dump_die($df->get_hash());

    $aa =~ s/^p\.//;
    my $key = join ".", $gene, $aa;
    my $study = $df->get_value("ID_STUDY");
    my $pmid = $df->get_value("Pubmed_PMID");
    my $real_pmid;
    if ($pmid) {
      $total_records_with_pmid++;
      $real_pmid = $pmid;
      $real_pmid{$real_pmid} = 1;
    } elsif ($study) {
      $pmid = "unknown_study_" . $study;
    } else {
      die "no ID_STUDY and no pmid";
      # ever happens?
    }

    if ($study) {
      # OK
    } elsif ($pmid) {
      # sometimes no study ID, but there is a PMID
      $study = "PMID_" . $pmid;
    } else {
#      $unusable{no_study_or_pmid}++;
      die;
    }

    $study2total{$study}++;
    if ($aa =~ /^([A-Z]+)\d+\1/) {
      # silent
      $silent{$key}++;
      $pmid2silent{$pmid}{$key}++;
      $study2silent{$study}++;
      $silent2pmid{$key}{$real_pmid}++ if $real_pmid;
    }
    $study2unique{$study}{$key} = 1;
    $study2pmid{$study}{$pmid} = 1;

    my $id_sample = $df->get_value("ID_sample") || die;
    $study2sample{$study}{$id_sample} = 1;

    my $id_tumor = $df->get_value("ID_tumour") || die;
    $study2tumor{$study}{$id_tumor} = 1;

    if (++$rows % 250000 == 0) {
      printf STDERR "processed %d rows in %d sec\n", $rows, time - $start;
    }
    last if $max and $rows >= $max;
  }

  printf STDERR "total:%d with_pid:%d (%.2f)\n",
    $total_records, $total_records_with_pmid,
      ($total_records_with_pmid * 100 / $total_records);
  if (%unusable) {
    printf STDERR "unusable rows:\n";
    foreach (sort keys %unusable) {
      printf STDERR "   %s: %d\n", $_, $unusable{$_};
    }
  }

  #
  # generate PubMed abstract info:
  #
  my $pmid_file = "pmid_query.txt";
  write_simple_file([keys %real_pmid], $pmid_file);
  system "pmid2abstract.pl -file $pmid_file";
  die "pmid2abstract.pl call failed" if $?;
  my $pmid_info = read_simple_file("pmid2abstract.tab", "-hash" => "pmid");

  # - intersection of PMIDs with many silent variants and
  #   many *recurrent* silent variants?
  # - single-cell sequencing detection:
  # - what about RATIO of silent to non-silent variants in a study?
  #   single-cell should have higher rates, right?
  # - OR: measure "sameness" of variants within a study?
  #   single-cell should have highly repetitive variants, others
  #   should be more random.  Ratio of unique variants to total variants?
  #
  my $rpt = new Reporter(
			 "-file" => "cosmic_study_stats.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   study
					   pmid
					   has_pmid
					   tumor_count
					   sample_count
					   ratio_sample_to_tumor

					   total
					   silent
					   unique
					   ratio_unique_to_total
					   ratio_silent_to_total
					   recurrent_silent_unique_count
					   recurrent_silent_total_count
					   ratio_recurrent_silent_total_to_pan

					   review_flag
					   review_reasons

					   suspicious_abstract_keywords
					   recurrent_silent_genes_abstract

					   recurrent_silent_details

					   title
					   abstract
					)
				      ],
			 "-auto_qc" => 1,
			);

  my $total_recurrent_silent_pan_cosmic = 0;
  foreach my $pmid (keys %pmid2silent) {
    my @keys = keys %{$pmid2silent{$pmid}};
    foreach my $key (@keys) {
      my $count = $pmid2silent{$pmid}{$key} || 0;
      if ($count >= $MIN_COUNT_FOR_RECURRENT) {
	$total_recurrent_silent_pan_cosmic += $count;
      }
    }
  }

  foreach my $study (keys %study2total) {
    my %r;
    $r{study} = $study;

    my @pmids;
    if ($study2pmid{$study}) {
      @pmids = keys %{$study2pmid{$study}};
    }
    my @real_pmids = grep {!/^unknown/} @pmids;

    $r{pmid} = join ",", @pmids;
    $r{has_pmid} = scalar @real_pmids ? 1 : 0;

    my $total = $study2total{$study} || die;
    my $unique = scalar keys %{$study2unique{$study}};

    $r{unique} = $unique;
    $r{total} = $total;

    my $ratio_unique_to_total = $unique / $total;
    $r{ratio_unique_to_total} = $ratio_unique_to_total;
    # lower ratios mean less variety in variants, possibly single-cell?
    # HOWEVER: some studies are of just one variant, e.g.

    my $silent_total = $study2silent{$study} || 0;
    $r{silent} = $silent_total;

    my $ratio_silent_to_total = $silent_total / $total;
    $r{ratio_silent_to_total} = $ratio_silent_to_total;

    my $total_recurrent_silent = 0;
    my %over;
    foreach my $pmid (@pmids) {
      my @keys = keys %{$pmid2silent{$pmid}};
      foreach my $key (@keys) {
	my $count = $pmid2silent{$pmid}{$key} || 0;
	if ($count >= $MIN_COUNT_FOR_RECURRENT) {
	  $over{$key} = $count;
	  $total_recurrent_silent += $count;
	}
      }
    }

    my $count_unique_recurrent_silent = 0;
    my $details = "";

    if (%over) {
      $count_unique_recurrent_silent = scalar keys %over;
      $details = join ", ", map {$_ . "=" . $over{$_}} sort keys %over;
    }
    $r{recurrent_silent_total_count} = $total_recurrent_silent;
    $r{recurrent_silent_unique_count} = $count_unique_recurrent_silent;
    $r{recurrent_silent_details} = $details;

    my $ratio_recurrent_silent_total_to_pan = $total_recurrent_silent / $total_recurrent_silent_pan_cosmic;
    $r{ratio_recurrent_silent_total_to_pan} = $ratio_recurrent_silent_total_to_pan;
    # TO DO: maybe exponential notation?

    my $sc = scalar keys %{$study2sample{$study}};
    my $tc = scalar keys %{$study2tumor{$study}};
    $r{sample_count} = $sc;
    $r{tumor_count} = $tc;
    my $ratio_sample_to_tumor = $sc / $tc;
    $r{ratio_sample_to_tumor} = $ratio_sample_to_tumor;

    if (@real_pmids) {
      my @title;
      my @abstract;
      foreach my $pmid (@real_pmids) {
	my $info = $pmid_info->{$pmid} || die "no PMID cache info for $pmid";
	push @title, $info->{title};
	push @abstract, $info->{abstract};
      }

      $r{title} = join "|", @title;
      $r{abstract} = join "|", @abstract;

      my %suspicious;
      foreach my $abstract (@abstract) {
	foreach my $term (@SUSPICIOUS_ABSTRACT_TERMS) {
	  $suspicious{$term} = 1 if $abstract =~ /$term/i;
	}
      }
      $r{suspicious_abstract_keywords} = join ",", sort keys %suspicious;

      my %hits;
      foreach my $gene (sort keys %over) {
	# recurrent silent
	$hits{$gene} = 1 if grep {/$gene/} @abstract;
      }
      $r{recurrent_silent_genes_abstract} = join ",", sort keys %hits;
      # never seems to actually happen.
      # what if extended to recurrent missense?
    } else {
      $r{title} = $r{abstract} = "no PMID";
      $r{suspicious_abstract_keywords} = "";
      $r{recurrent_silent_genes_abstract} = "";
    }

    my $review_flag;
    my @review_reasons;

    if ($count_unique_recurrent_silent >= 5) {
      push @review_reasons, "absolute_recurrent_silent_count";
      # FIX ME: GLOBAL VAR
      push @review_reasons, "recurrent_variants_and_no_PMID" unless @real_pmids;
    }
    push @review_reasons, "high_ratio_sample_to_tumor" if $ratio_sample_to_tumor > 2.0;
    # allow to for e.g. WGS + WES sequencing

    push @review_reasons, "low_ratio_unique_to_total" if $ratio_unique_to_total < 0.33;
    # e.g. PMID 26672766 = intratumoral heterogeneity
    # FIX ME: GLOBAL VAR

    push @review_reasons, "high_ratio_silent_to_total" if $ratio_silent_to_total >= 0.33;
    # FIX ME: global var

    push @review_reasons, "suspicious_abstract_keyword" if $r{suspicious_abstract_keywords};

#    push @review_reasons, "high_ratio_recurrent_silent_total_to_pan" if $ratio_recurrent_silent_total_to_pan >= 0.03;
    push @review_reasons, "high_ratio_recurrent_silent_total_to_pan" if $ratio_recurrent_silent_total_to_pan >= 0.01;
    # study represent 1% or more of total recurrent variants in COSMIC

    $review_flag = @review_reasons ? 1 : 0;

    $r{review_flag} = $review_flag;
    $r{review_reasons} = join ",", @review_reasons;

    $rpt->end_row(\%r);
  }
  $rpt->finish();

  $rpt = new Reporter(
		      "-file" => "cosmic_silent_multiple_pmid.tab",
		      "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   variant
					   gene
					   AA
					   is_reviewable
					   is_reportable
					   pmid_count
					   pmids
					)
				      ],
			 "-auto_qc" => 1,
			);

  my $rr = new ReviewableReportable("-genome" => $FLAGS{genome} || die "-genome");

  foreach my $variant (sort keys %silent2pmid) {
    my @pmids = sort {$a <=> $b} keys %{$silent2pmid{$variant}};
    if (@pmids > 1) {
      my %r;
      $r{variant} = $variant;
      @r{qw(gene AA)} = split /\./, $variant;
      $r{pmid_count} = scalar @pmids;
      $r{pmids} = join ", ", @pmids;

      my $gene = $r{gene} || die;
      $r{is_reviewable} = $rr->is_reviewable($gene) || 0;
      $r{is_reportable} = $rr->is_reportable($gene) || 0;

      $rpt->end_row(\%r);
    }
  }
  $rpt->finish();



}
