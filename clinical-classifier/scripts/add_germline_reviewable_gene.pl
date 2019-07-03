#!/bin/env perl
# add gene(s) to germline reviewable list
# /nfs_exports/genomes/1/PCGP/BucketIntermediate/germlineSNVs/genelists/cancer_related_genes_dec19.2013.lst
#
# TO DO:

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use File::Basename;
use List::Util qw(min max);
use File::Copy;

use Cluster;
use CommandLineRebuilder;
use GenomicRangeFinder;
use TdtConfig;
use FileUtils qw(read_simple_file write_simple_file);
use MiscUtils qw(dump_die);
use RefFlatFile;
use DelimitedFile;
use GeneSymbolMapper;

my $FLANKING_BUFFER_NT = 1000;
# expand isoform-associated ranges by this # of nt on either side

my $MANUAL_GENES = init_manual_genes();

my %FLAGS;

my @clopts = (
	      "-genome=s",
	      "-ranges=s",
	      # cancer-related regions file
	      # (manual override: taken from config by default)

	      "-gene=s",
	      "-genes=s",
	      "-genes-alt=s",
	      "-genes-all=s",

	      "-refflat=s",

	      "-out=s",
	      "-keep-partial",

	      "-hack",

	      "-hgnc=s",
	      "-gene-info=s",

	      "-prune-to-new",
	      "-prune-ignore=s",
	     );

GetOptions(
	   \%FLAGS,
	   @clopts
	  );

if ($FLAGS{hack}) {
  hack();
  exit(0);
}

my $genome = $FLAGS{genome} || die "-genome";
my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";

my $f_ranges = $FLAGS{ranges};
unless ($f_ranges) {
  $f_ranges = $config_genome->{CLINCLS_CANCER_RELATED_GENES_FILE} || die;
}
printf STDERR "ranges: %s\n", $f_ranges;

if ($FLAGS{"prune-to-new"}) {
  prune_to_new();
  exit(0);
}

my $gsm = get_gsm();

my $grf_gold = new GenomicRangeFinder();
open(GLTMP, $f_ranges) || die;
while (<GLTMP>) {
  chomp;
  my @f = split /\t/, $_;
  die unless @f == 2;
  my ($gene, $loc) = @f;
  $gsm->add_gene("-gene" => $gene);
  $grf_gold->add(
		 "-range" => $loc,
		 "-value" => {
			      "gene" => $gene,
			      "location" => $loc
			     },
		);
}


my $f_refflat = $FLAGS{refflat} || $config_genome->{REFSEQ_REFFLAT} || die "no refflat";
my $rff = new RefFlatFile();
$rff->parse_file("-refflat" => $f_refflat,
		 "-type" => "refgene");
my $rows = $rff->get_rows();
my %gene2rows;

foreach my $r (@{$rows}) {
  my $gene = $r->{gene} || next;
  push @{$gene2rows{$gene}}, $r;
  # keep "sharp" symbols as-is to separate different mapping groups
}

my @rows;
if (my $one_gene = $FLAGS{gene}) {
  my %r;
  $r{gene} = $one_gene;
  push @rows, \%r;
} elsif (my $f = $FLAGS{genes}) {
  my $set = read_simple_file($f);
  foreach my $g (@{$set}) {
    my %r;
    $r{gene} = $g;
    push @rows, \%r;
  }
} elsif (my $fn = $FLAGS{"genes-alt"}) {
  my $df = new DelimitedFile("-file" => $fn,
			     "-headers" => 1,
			     );
  while (my $row = $df->get_hash()) {
    foreach my $f (qw(gene alternate)) {
      die "row missing field $f" unless $row->{$f};
    }
    push @rows, $row;
  }
} else {
  die "specify -gene GENE or -genes LISTFILE\n";
}


my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
#my $outfile = $FLAGS{out} || basename($f_ranges) . ".patched.txt";
my $outfile = $FLAGS{out} || sprintf "cancer_related_genes_%d_%02d_%02d.txt", 1900+$year, $mon + 1, $mday;

copy($f_ranges, $outfile) || die;
open(OUT, ">>" . $outfile) || die;

my $count_ok = 0;
my $count_manual = 0;
my $count_duplicate = 0;
my $count_fail = 0;

foreach my $row (@rows) {

  my @genes;
  push @genes, $row->{gene} || die "no gene";
  push @genes, $row->{alternate} if $row->{alternate} and $row->{alternate} ne $row->{gene};

  my $processed;

  #
  #  first check if any provided symbol matches an existing entry:
  #
  foreach my $gene (@genes) {
    $gene =~ s/^\s+//;
    $gene =~ s/\s$//;
    my $gene_cleaned = clean_sharp_gene_symbol($gene);
    die join ",", $gene, $gene_cleaned if $gene ne $gene_cleaned;
    if ($gsm->contains($gene)) {
      # same in both
      $count_duplicate++;
      $processed = 1;
      last;
    } elsif (my $result = $gsm->resolve("-symbol" => $gene)) {
      $count_duplicate++;
      printf STDERR "new symbol %s matches old symbol %s in current list\n",
	$gene_cleaned, $result;
      $processed = 1;
      last;
    }
  }

  my $gene_idx = 0;
  foreach my $gene (@genes) {
    last if $processed;

    my @wanted_sets;
    foreach my $suffix ("", "_locPar") {
      my $key = $gene . $suffix;
      if (my $wanted = $gene2rows{$key}) {
	push @wanted_sets, $wanted;
	last;
      }
    }
    unless (@wanted_sets) {
      # refFlat "sharp" may contain multiple mappings
      foreach my $site ("A" .. "Z") {
	my $key = sprintf '%s_loc%s', $gene, $site;
	if (my $wanted = $gene2rows{$key}) {
	  push @wanted_sets, $wanted;
	}
      }
    }


    # FIX ME: refFlat "sharp" cases w/multiple mappings!
    # jump off that bridge when we come to it

    # unless ($wanted) {
    #   printf STDERR "ERROR: no records found for \"%s\" in %s\n", $gene, $f_refflat;
    #   $count_fail++;
    #   next;
    # }

    foreach my $wanted (@wanted_sets) {
      my %chr;
      my %strand;
      my @starts;
      my @ends;
      foreach my $r (@{$wanted}) {
	$chr{$r->{chrom}}++;
	$strand{$r->{strand}}++;
	push @starts, $r->{txStart} || die;
	push @ends, $r->{txEnd} || die;
      }

      die "multiple chroms" if scalar keys %chr > 1;
      die "multiple strands" if scalar keys %strand > 1;

      my ($chr) = keys %chr;
      $chr =~ s/^chr//;
      my $start = min(@starts) - $FLANKING_BUFFER_NT;
      my $end = max(@ends) + $FLANKING_BUFFER_NT;

      my $hits = $grf_gold->find(
				 "-chrom" => $chr,
				 "-base" => $start,
				 "-end" => $end
				);

      my $usable = 1;
      if (@{$hits}) {
	my @old_genes = map {$_->{gene}} @{$hits};
	my @old_match = grep {$_ eq $gene} @old_genes;

	if (@old_match) {
	  # duplicate of existing row
	  $count_duplicate++;
	  $usable = 0;
	} else {
	  printf STDERR "WARNING: new gene %s overlaps with existing entry for %s\n", $gene, join ",", @old_genes;
	  # these should be checked closely:
	  #
	  # (A) might be legit:
	  #   1. different genes in same region:
	  #      - GP1BB/SEPT5
	  #      - VREPB1/IGL
	  #   2. genes where flanking regions happen to overlap
	  # (B) OTOH: maybe new list prefers a different symbol than
	  #   found in old list?
	}
      }
      if ($usable) {
	printf OUT "%s\t%s:%d-%d\n",
	  $gene, $chr, $start, $end;
	$count_ok++;
      }

      $processed = 1;
      # lookup resolved using this symbol

    } # $wanted

    $gene_idx++;
  }  # $gene

  unless ($processed) {
    my $manual;
    foreach my $g (@genes) {
      $manual = $MANUAL_GENES->{$g};
      last if $manual;
    }

    if ($manual) {
      printf STDERR "manual entry for %s\n", join ",", @genes;

      my $start = ($manual->{start} || die) - $FLANKING_BUFFER_NT;
      my $end = ($manual->{end} || die) + $FLANKING_BUFFER_NT;

      printf OUT "%s\t%s:%d-%d\n",
	$manual->{gene},
	  $manual->{chr},
	    $start, $end;
	$count_manual++;
    } else {
      printf STDERR "ERROR: no records found for %s in %s\n", join(",", @genes), $f_refflat;
      $count_fail++;
    }
  }


}

close OUT;

printf STDERR "skipped %d duplicates\n", $count_duplicate if $count_duplicate;
printf STDERR "added %d new manual records\n", $count_manual if $count_manual;
my $fatal;

if ($count_ok or $count_manual) {
  printf STDERR "added %d new auto entries\n", $count_ok;
} else {
  printf STDERR "ERROR: no records successfully added\n";
  $fatal = 1;
}

if ($count_fail) {
  printf STDERR "ERROR: %d failed record(s)\n", $count_fail;
  $fatal = 1;
}

if ($fatal) {
  if ($FLAGS{"keep-partial"}) {
    printf STDERR "fatal error, but keeping output anyway\n";
  } else {
    printf STDERR "fatal error, deleting output\n";
    unlink $outfile;
  }
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
    printf STDERR "stripping %s => %s\n", $before, $gene;
  }
  return $gene;
}

sub init_manual_genes {
  #
  # raw coordinates for manual entries (e.g. NG_)
  # - add NCBI coordinates, flanking sequence will be added later

  my %manual;
  my $row;

  $row = {
	  "gene" => "IGHA1",
	  "chr" => 14,
	  "start" => 106173505,
	  "end" => 106175001,
	 };
  $manual{$row->{gene}} = $row;
  # http://www.ncbi.nlm.nih.gov/gene/3493

  $row = {
	  "gene" => "IGHA2",
	  "chr" => 14,
	  "start" => 106053274,
	  "end" => 106054731
	 };
  $manual{$row->{gene}} = $row;
  # http://www.ncbi.nlm.nih.gov/gene/3494

  $row = {
	  "gene" => "IGHE",
	  "chr" => 14,
	  "start" => 106066403,
	  "end" => 106068064
	 };
  $manual{$row->{gene}} = $row;
  # http://www.ncbi.nlm.nih.gov/gene/3497

  $row = {
	  "gene" => "IGHG1",
	  "chr" => 14,
	  "start" => 106207810,
	  "end" => 106209407
	 };
  $manual{$row->{gene}} = $row;
  # http://www.ncbi.nlm.nih.gov/gene/3500

  $row = {
	  "gene" => "IGHG2",
	  "chr" => 14,
	  "start" => 106109540,
	  "end" => 106111126
	 };
  $manual{$row->{gene}} = $row;
  # http://www.ncbi.nlm.nih.gov/gene/3501

  $row = {
	  "gene" => "IGHG3",
	  "chr" => 14,
	  "start" => 106232251,
	  "end" => 106237742
	 };
  $manual{$row->{gene}} = $row;
  # http://www.ncbi.nlm.nih.gov/gene/3502

  $row = {
	  "gene" => "IGHG4",
	  "chr" => 14,
	  "start" => 106090813,
	  "end" => 106092402
	 };
  $manual{$row->{gene}} = $row;
  # http://www.ncbi.nlm.nih.gov/gene/3503

  $row = {
	  "gene" => "IGHM",
	  "chr" => 14,
	  "start" => 106318298,
	  "end" => 106322322
	 };
  $manual{$row->{gene}} = $row;
  # http://www.ncbi.nlm.nih.gov/gene/3507

  $row = {
	  "gene" => "IGKC",
	  "chr" => 2,
	  "start" => 89156874,
	  "end" => 89157196
	 };
  $manual{$row->{gene}} = $row;
  # http://www.ncbi.nlm.nih.gov/gene/3514

  $row = {
	  "gene" => "TRAC",
	  "chr" => 14,
	  "start" => 23016447,
	  "end" => 23021075
	 };
  $manual{$row->{gene}} = $row;
  # http://www.ncbi.nlm.nih.gov/gene/28755

  return \%manual;
}

sub get_gsm {
  my $config_species = TdtConfig::readConfig('species', "Homo_sapiens") || die "can't find hs species config";
  my $hgnc = $FLAGS{hgnc};
  my $gene_info = $FLAGS{"gene-info"};

  unless ($hgnc) {
    $hgnc = $config_species->{HGNC} || die "no HGNC config";
  }
  unless ($gene_info) {
    $gene_info = $config_species->{ENTREZ_GENE_GENEINFO} || die "no config ENTREZ_GENE_GENEINFO";
  }

  my $gsm = new GeneSymbolMapper(
				 "-hgnc_file" => $hgnc,
				 "-eg_file" => $gene_info,
#				 "-hgnc_synonym_enable" => 0,
				 # disable for e.g. FAH/FANCA
				);

  return $gsm;
}

sub hack {
  my $gsm = get_gsm();
  die;
}

sub prune_to_new {
  my $new_genes = read_simple_file($FLAGS{genes} || die "-genes");
  my $all_genes = read_simple_file($FLAGS{"genes-all"} || die "-genes-all");
  my $gsm_main = get_gsm();
  my %ignore;
  if (my $ignore = $FLAGS{"prune-ignore"}) {
    %ignore = map {$_, 1} split /,/, $ignore;
  }

  my $gsm_new = new GeneSymbolMapper("-clone" => $gsm_main);
  foreach my $g (@{$new_genes}) {
    $gsm_new->add_gene("-gene" => $g);
  }

  my $gsm_all = new GeneSymbolMapper("-clone" => $gsm_main);
  foreach my $g (@{$all_genes}) {
    $gsm_all->add_gene("-gene" => $g);
  }

  my $gsm_orig = new GeneSymbolMapper("-clone" => $gsm_main);
  #
  # for each entry on the ORIGINAL list, keep if present on the new list,
  # otherwise discard:
  #
  open(IN, $f_ranges) || die;
  my $outfile = basename($f_ranges) . ".trimmed";
  my $wf = new WorkingFile($outfile);
  my $fh = $wf->output_filehandle();
  my $fatal;
  while (my $line = <IN>) {
    my ($gene) = split /\t/, $line;
    if ($gsm_new->find($gene)) {
      $gsm_orig->add_gene("-gene" => $gene);
      print $fh $line;
    } else {
      if ($gsm_all->find($gene) or $ignore{$gene}) {
	# ok to drop
	printf STDERR "dropping gene %s\n", $gene;
      } else {
	# gene not confirmed reviewed before dropping
	printf STDERR "FATAL ERROR: gene %s not in genes-all\n", $gene;
	$fatal = 1;
      }
    }
  }
  die "fatal error" if $fatal;
  $wf->finish();

  #
  #  make a list of NEW genes not already present in the old list:
  #
  my @new_genes;
  foreach my $g (@{$new_genes}) {
    unless ($gsm_orig->find($g)) {
      push @new_genes, $g;
    }
  }

  my $f_new_to_add = "new_genes_to_add.txt";
  write_simple_file(\@new_genes, $f_new_to_add);

  die scalar @new_genes;


}
