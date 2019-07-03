#!/bin/env perl
# collate tumor suppressor/oncogene annotations from various sources
# MNE 4/2018
# /home/medmonso/work/jz/2017_11_02_pecan_pie_paper/ts_onco_annotations/combined
#
# TO DO:
#  - QC: do any germline-reviewable genes lack this annotation?

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list get_hash_option);
use FileUtils qw(read_simple_file);
use DelimitedFile;
use Reporter;
use TdtConfig;
use GeneSymbolMapper;
use MiscUtils qw(log_message);
use TSOncoDB;

my %FLAGS;
my @clopts = (
	      "-build",
	      "-genome=s",

	      "-oncokb=s",
	      "-tsgene=s",

	      "-yanling=s",
	      # T:\ZhangLab\JinghuiZhang\Papers\PeCanPIE\meeting_ppt\2018_04_04_gene_lists\cancer_gene_census_tier1_vs_reviewable-jz-yanling.xlsx
	      # - filtered to:
	      #   (a) usable sources (not reviewable overlap, etc
	      #   (b) LoF/hotspot annotation

	      "-debug-ram",

	      "-append=s",

	      "-qc-germline",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

if ($FLAGS{build}) {
  build_db();
} elsif ($FLAGS{"qc-germline"}) {
  qc_germline();
} else {
  die "?\n";
}

sub build_db {
  # build up database from various sources

  my @in_rows;

  #
  # parse basic info from all databases into standardized format:
  #
  my $df;

  my %excel_damage = (
		      "Dec-01" => "DEC1",
		     );
  # damaged data in TSGene

  if (my $extra = $FLAGS{append}) {
    # use existing database rather than building from scratch
    my $config_species = TdtConfig::readConfig('species', "Homo_sapiens") || die "can't find hs species config";
    my $f_db = $config_species->{TS_ONCO_DB} || die;

    #
    # existing data:
    #
    my $df = new DelimitedFile("-file" => $f_db,
			     "-headers" => 1,
			     );
    while (my $row = $df->get_hash()) {
      my $sym = $row->{sym} || die;
      foreach my $f (qw(is_ts is_onco)) {
	if (my $sources = $row->{$f}) {
	  foreach my $source (split /,/, $sources) {
	    my %r;
	    $r{sym} = $sym;
	    $r{source} = $source;
	    $r{is_ts} = "";
	    $r{is_onco} = "";

	    $r{$f} = 1;
	    push @in_rows, \%r;
	  }
	}
      }
    }

    #
    # new data:
    #
    $df = new DelimitedFile("-file" => $extra,
			     "-headers" => 1,
			     );
    foreach my $f (qw(sym source is_ts is_onco)) {
      die "where is $f" unless exists $df->headers->{$f};
    }
    while (my $row = $df->get_hash()) {
      push @in_rows, $row;
    }


  } else {
    #
    # PeCan PIE paper: Yanling/JZ reviews/additions:
    #
    my $genome = $FLAGS{genome} || die "-genome";
    my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";

    my $f_yanling = $FLAGS{yanling} || die "-yanling";
    $df = new DelimitedFile(
			    "-file" => $f_yanling,
			    "-headers" => 1,
			   );
    while (my $row = $df->get_hash()) {
      my %r;
      my $class = $row->{Yanling_comment} || die;
      $r{sym} = $row->{gene} || die;
      $r{source} = "SJ_PIE_review_2018_04";
      my $is_ts = "";
      my $is_onco = "";
      if ($class eq "Hotspot") {
	$is_onco = 1;
      } elsif ($class =~ /^loss\-of\-function/i) {
	# control character
	$is_ts = 1;
      } else {
	die $class;
      }
      $r{is_ts} = $is_ts;
      $r{is_onco} = $is_onco;
      push @in_rows, \%r;
    }

    #
    # SJ germline reportable gene annotations (SNV/indel/CNV):
    #
    my $f_gl_ts = $config_genome->{CLINCLS_GL_REPORTABLE_GENE_ANNOTATION_FILE} || die;
    $df = new DelimitedFile(
			    "-file" => $f_gl_ts,
			    "-headers" => 1,
			   );
    while (my $row = $df->get_hash()) {
      my %r;
      $r{sym} = $row->{Gene} || die;
      my $tg = $row->{TruncationGold};
      my $is_ts = "";
      my $is_onco = "";
      if ($tg eq "Y") {
	$is_ts = 1;
      } elsif ($tg eq "N") {
	$is_onco = 1;
      } elsif ($tg eq "2") {
	# both
	$is_ts = $is_onco = 1;
      } else {
	die "unhandled value $tg";
      }
      $r{is_ts} = $is_ts;
      $r{is_onco} = $is_onco;
      $r{source} = "SJ_germline";
      push @in_rows, \%r;
    }

    #
    # SJ manual somatic gene additions:
    #
    my $f_manual = $config_genome->{CLINCLS_GENES_MANUAL_FILE} || die;
    $df = new DelimitedFile(
			    "-file" => $f_manual,
			    "-headers" => 1,
			   );
    while (my $row = $df->get_hash()) {
      my %r;
      $r{sym} = $row->{Gene} || die;
      my $class = $row->{Class} || die;
      my ($is_ts, $is_recurrent) = parse_gene_class($class);
      $r{is_ts} = $is_ts;
      $r{is_onco} = $is_recurrent;
      $r{source} = "SJ_manual";
      push @in_rows, \%r;
    }

    #
    # SJ legacy medal ceremony annotations from JZ:
    #
    my %param2source = (
			"CLINCLS_CANCER_GENE_LIST_FILE" => "SJ_JZ_cancer",
			"CLINCLS_NON_CANCER_GENE_LIST_HC_FILE" => "SJ_JZ_noncancer_high_confidence",
			"CLINCLS_NON_CANCER_GENE_LIST_LC_FILE" => "SJ_JZ_noncancer_low_confidence"
			# high and low confidence
		       );

    foreach my $p (sort keys %param2source) {
      my $fn = $config_genome->{$p} || die;
      my $df = new DelimitedFile(
				 "-file" => $fn,
				 "-headers" => 1,
				);
      my $source = $param2source{$p} || die;
      while (my $row = $df->get_hash()) {
	my $gene = get_field($row, "Gene");
	my $gene_class = get_field($row, "Class");
	my ($is_ts, $is_recurrent) = parse_gene_class($gene_class);
	next if (not($gene =~ /\w/) and not($gene_class =~ /\w/));
	# blank lines

	my %r;
	$r{sym} = $gene;
	$r{source} = $source;
	$r{is_ts} = $is_ts;
	$r{is_onco} = $is_recurrent;
	push @in_rows, \%r;
      }
    }

    # TSGene:
    my $f_tsgene = $FLAGS{tsgene} || die "-tsgene";
    $df = new DelimitedFile("-file" => $f_tsgene,
			    "-headers" => 1,
			   );

    while (my $row = $df->get_hash()) {
      my %r;
      my $gene = $row->{GeneSymbol} || die;
      $gene = $excel_damage{$gene} || $gene;
      $r{sym} = $gene;
      $r{is_ts} = 1;
      # all tumor suppressors in this db
      $r{is_onco} = "";
      # no opinion
      $r{source} = "tsgene";
      push @in_rows, \%r;
    }

    # oncokb:
    my $f_oncokb = $FLAGS{oncokb} || die "-oncokb";
    $df = new DelimitedFile(
			    "-file" => $f_oncokb,
			    "-headers" => 1,
			   );

    while (my $row = $df->get_hash()) {
      my %r;
      $r{sym} = $row->{"Hugo Symbol"} || die;
      $r{is_ts} = oncokb_value_to_flag($row->{"OncoKB TSG"});
      $r{is_onco} = oncokb_value_to_flag($row->{"OncoKB Oncogene"});
      $r{source} = "oncokb";
      if ($r{is_ts} or $r{is_onco}) {
	# some genes are in list but have no ts/onco annotations, e.g. ABL2
	push @in_rows, \%r;
      }
    }
  }


  # attempt to standardize all symbols to HGNC,
  # e.g. SJ legacy files use MLL3 for KMT2C:
  my $gsm = new_gsm();
  my $idx_approved = $gsm->hgnc->index_approved();
  foreach my $g (keys %{$idx_approved}) {
    $gsm->add_gene("-gene" => $g);
  }

  foreach my $row (@in_rows) {
    my $sym = $row->{sym} || die;
    unless ($gsm->contains($sym)) {
      if (my $new = $gsm->resolve("-symbol" => $sym)) {
	printf STDERR "update %s to %s in %s\n", $sym, $new, $row->{source};
	$row->{sym} = $new;
      } else {
	printf STDERR "ERROR: can't resolve %s in %s, assuming correct\n", $sym, $row->{source};
      }
    }
  }

  # intermediate report: separate entry for each source row,
  # not unique by gene, but more easily parsable for some purposes:
  my $rpt = new Reporter(
			 "-file" => "ts_onco_intermediate.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   sym
					   source
					   is_ts
					   is_onco
					)
				      ],
			 "-auto_qc" => 1,
			);
  foreach my $r (@in_rows) {
    $rpt->end_row($r);
  }
  $rpt->finish();


  #
  # merge into one result per gene, with source list under each flag:
  #
  my %track;
  foreach my $r (@in_rows) {
    my ($sym, $source, $is_ts, $is_onco) = @{$r}{qw(sym source is_ts is_onco)};
#    printf STDERR "track %s\n", join ",", $sym, $source, $is_ts, $is_onco;
    $track{$sym}{is_ts}{$source} = 1 if $is_ts;
    $track{$sym}{is_onco}{$source} = 1 if $is_onco;
  }
  $rpt = new Reporter(
			 "-file" => "ts_onco_merged.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   sym
					   is_ts
					   is_onco
					)
				      ],
			 "-auto_qc" => 1,
			);
  foreach my $sym (sort keys %track) {
    my %r;
    $r{sym} = $sym;
    $r{is_ts} = $r{is_onco} = "";
    $r{is_ts} = join ",", sort keys %{$track{$sym}{is_ts}} if $track{$sym}{is_ts};
    $r{is_onco} = join ",", sort keys %{$track{$sym}{is_onco}} if $track{$sym}{is_onco};
    $rpt->end_row(\%r);

  }
  $rpt->finish();

}

sub oncokb_value_to_flag {
  my ($v) = @_;
  my $flag = "";
  # no opinion
  if (defined $v) {
    if ($v eq "Yes") {
      $flag = 1;
    } elsif ($v eq "") {
      $flag = "";
    } else {
      die sprintf 'unhandled value "%s"', $v;
    }
  }
  return $flag;
}

sub parse_gene_class {
  my ($string) = @_;
  my @things = split /\+/, $string;
  my $is_ts = "";
  my $is_recurrent = "";
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

sub get_field {
  my ($row, $label) = @_;
  my $value = $row->{$label};
  die "no value for $label" unless defined $value;
  return $value;
}

sub new_gsm {

  my $config_species = TdtConfig::readConfig('species', "Homo_sapiens") || die "can't find hs species config";
  my $f_hgnc = $config_species->{HGNC} || die "no HGNC config";

  ram_debug("GSM start");

  my $gsm = new GeneSymbolMapper(
				 "-hgnc_file" => $f_hgnc,
				 "-hgnc_lite" => 1,
				 "-enable_entrez_gene" => 0,
				);
  ram_debug("GSM end");
  return $gsm;
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

sub qc_germline {
  my $tsdb = new TSOncoDB();
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";

  my $reportable = read_simple_file($config_genome->{CLINCLS_GERMLINE_REPORTABLE_GENES} || die);

  my %missing;
  foreach my $gene (@{$reportable}) {
    my $row = $tsdb->find_row("-gene" => $gene);
    $missing{reportable}{$gene} = 1 unless $row;
  }

  my $f_rev = $config_genome->{CLINCLS_CANCER_RELATED_GENES_FILE} || die;
  open(IN, $f_rev) || die;
  my %rev;
  while (<IN>) {
    chomp;
    my @f = split /\t/, $_;
    die unless @f == 2;
    $rev{$f[0]} = 1;
  }
  foreach my $gene (sort keys %rev) {
    my $row = $tsdb->find_row("-gene" => $gene);
    $missing{reviewable}{$gene} = 1 unless $row;
  }

  foreach my $type (sort keys %missing) {
    printf STDERR "missing %s: %d\n", $type, scalar keys %{$missing{$type}};
    foreach my $gene (sort keys %{$missing{$type}}) {
      printf STDERR "  %s\n", $gene;
    }
  }

}
