#!/bin/env perl
#
# parse COSMIC variants, attempting to map details to genomic space
#
# MNE 3/2015

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;

use MiscUtils qw(dump_die log_message);
use DelimitedFile;
use NucleotideSubstitutionParser;
use FAI;
use TdtConfig;
use TabixPrep;
use TARTANUtils qw(tartan_genome_put_helper);

use File::Basename;

my %FLAGS;
GetOptions(\%FLAGS,
	   "-genome=s",

	   "-cosmic=s",
	   # optional (taken from config)
	   # /nfs_exports/genomes/1/projects/ClinicalSeq/medal_config/2013_09_20/CosmicMutantExport_v66_250713_cleaned_2013_12_20.tab

	   "-tabix",
	   "-max=i",
	   "-keep-tempfiles",
	   "-ping=i",

	   "-tartan-index=s",
	   # output dir
	  );

my $genome = $FLAGS{genome} || die "-genome";
my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
#dump_die($config_genome);

if (my $out_dir = $FLAGS{"tartan-index"}) {
  tartan_genome_put_helper(
			   "-genome" => $genome,
			   "-subdir" => "CLINICAL/COSMIC/COSMIC_tabix",
			   "-out" => $out_dir,
			  );
  exit(0);
}

my $genome_fasta = $config_genome->{FASTA} || die;
my $fai = new FAI("-fasta" => $genome_fasta);

my $cosmic = $FLAGS{cosmic};
#$cosmic = $config_genome->{CLINCLS_COSMIC_MUTANT_FILE} unless $cosmic;
# require user to be specific:
# for Protein Paint, JZ wants this version:
# /nfs_exports/genomes/1/projects/ClinicalSeq/COSMIC_confirmed_somatic_gene_mod.out

die "-cosmic" unless $cosmic;

printf STDERR "COSMIC: %s\n", $cosmic;

my $TABIX_MODE = $FLAGS{tabix};

my $df = new DelimitedFile(
			   "-file" => $cosmic,
			   "-headers" => 1,
			  );

my $outfile = basename($cosmic) . ".genomic.tab";

my $outfile_fail = basename($cosmic) . ".genomic_failed.tab";

my @extra = qw(
		Chr
		WU_HG19_Pos
		ReferenceAllele
		MutantAllele
		genomic_parsable
		genomic_parsable_type
		genomic_parsable_notes
	     );
# put position before comments as blank fields break sort call in tabix mode

my $rpt;
my $tp;
if ($TABIX_MODE) {
  my @headers = @extra;
  push @headers, @{$df->headers_raw};
  # TO DO:
  # - maybe strip out most columns for smallest possible file size?
  # - remove duplicate rows for the same variant?

  $tp = new TabixPrep(
		      "-outfile" => $outfile . ".gz",
		      "-headers" => \@headers,
		      "-header_chr" => "Chr",
		      "-header_start" => "WU_HG19_Pos",
		      "-delete_tempfiles" => ($FLAGS{"keep-tempfiles"} ? 0 : 1),
		     );


  #$tp->add_row("-row" => $row);
  #$tp->finish();
} else {
  $rpt = $df->get_reporter(
			   "-file" => $outfile,
			   "-extra" => \@extra,
			  );
}

my $rpt_fail = $df->get_reporter(
				 "-file" => $outfile_fail,
				 "-extra" => \@extra,
				);

my $nsp = new NucleotideSubstitutionParser();

my $missing_genomic_pos = 0;

my ($F_POS, $F_STRAND);
my $row_count = 0;
my $max = $FLAGS{max};

my $ping = $FLAGS{ping} || 500;

while (my $row = $df->get_hash()) {
#  dump_die($row, "debug", 1);
  $row_count++;

  log_message("processed $row_count") if $ping and $row_count % $ping == 0;
  last if $max and $row_count > $max;

  my $parsable = 0;
  my $ref_base = "";
  my $var_base = "";
  my $genomic_parsable_type = "";
  my $chr = "";
  my $start = "";
  my @notes;

  unless ($F_POS) {
    $F_POS = get_position_column($row);
    $F_STRAND = get_strand_column($row);
  }
  my $pos_raw = $row->{$F_POS};

  if ($pos_raw) {
    #
    #  genomic position available
    #
    my @p1 = split /:/, $pos_raw;
    die unless @p1 == 2;
    my $range;
    ($chr, $range) = @p1;

    if ($chr >= 1 and $chr <= 22) {
      # normal
    } elsif ($chr = 23) {
      $chr = "X";
      # bleah
    } else {
      die "unhandled chr $chr";
    }

    my @r = split /\-/, $range;
    die unless @r == 2;
    my $end;
    ($start, $end) = @r;
#    my $cds_raw = $row->{"Mutation CDS"} || dump_die($row, "no Mutation CDS");
    my $cds_raw = $row->{"Mutation CDS"} || "";
    # can be blank in v72
    my $strand = $row->{$F_STRAND} || die;

    $nsp->auto_strand_fix($strand);
    if ($nsp->parse($cds_raw)) {
      $parsable = 1;

      if ($nsp->is_substitution()) {
	$genomic_parsable_type = "substitution";
	$ref_base = $nsp->reference_sequence();
	$var_base = $nsp->variant_sequence();
      } elsif ($nsp->is_deletion()) {
	$genomic_parsable_type = "deletion";
	$ref_base = $nsp->reference_sequence();

	unless ($ref_base =~ /^[acgt]+$/i) {
	  # c.643_657del15
	  my $len = $end - $start + 1;
	  #	$ref_base = "N" x $len;
	  $ref_base = $fai->get_chunk(
				      "-id" => $chr,
				      "-start" => $start,
				      "-length" => $len
				     );
	  push @notes, "generated_reference_sequence";
	}
	$var_base = "-";
      } elsif ($nsp->is_insertion()) {
	$genomic_parsable_type = "insertion";
	$ref_base = "-";
	$var_base = $nsp->variant_sequence();
	if ($var_base) {
	  unless ($var_base =~ /^[acgt]+$/i) {
	    $parsable = 0;
	    push @notes, "invalid_inserted_bases";
	  }
	} else {
	  $parsable = 0;
	  push @notes, "no_inserted_sequence";
	}
      } elsif ($nsp->is_complex_indel()) {
	$genomic_parsable_type = "complex";
	push @notes, "complex";
	$ref_base = $nsp->reference_sequence();
	$var_base = $nsp->variant_sequence();

	if ($ref_base eq $nsp->REF_SEQ_UNSPECIFIED) {
	  my $len = $end - $start + 1;
	  $ref_base = $fai->get_chunk(
				      "-id" => $chr,
				      "-start" => $start,
				      "-length" => $len
				     );
	  push @notes, "generated_reference_sequence";
	}
      } else {
	die;
      }
    }

    if ($ref_base and $ref_base ne "-") {
      #
      # sanity-check reference sequence
      #
      my $sanity = $fai->sanity_check(
				      "-sequence" => $ref_base,
				      "-id" => $chr,
				      "-start" => $start,
				     );
      unless ($sanity) {
	$parsable = 0;
	push @notes, "reference_sequence_mismatch";
      }


      unless ($ref_base =~ /^[acgt]+$/i) {
	$parsable = 0;
	push @notes, "reference_sequence_invalid_nucleotide";
      }
    }

  } else {
    $missing_genomic_pos++;
    push @notes, "no_genomic_coordinates";
  }

  if ($parsable) {
    foreach my $seq ($ref_base, $var_base) {
      next if $seq eq "-";
      # insertion/deletion: OK
      if ($seq =~ /^[ACGT]+$/) {
	# OK
      } else {
	$parsable = 0;
	push @notes, "non_ACGT_sequence";
      }
    }
  }

  $row->{genomic_parsable} = $parsable;
  $row->{genomic_parsable_type} = $genomic_parsable_type;
  $row->{genomic_parsable_notes} = join ",", @notes;
  $row->{Chr} = $chr;
  $row->{WU_HG19_Pos} = $start;
  $row->{ReferenceAllele} = $ref_base;
  $row->{MutantAllele} = $var_base;

  if ($parsable) {
    if ($TABIX_MODE) {
      $tp->add_row("-row" => $row);
    } else {
      $rpt->end_row($row);
    }
  } else {
    $rpt_fail->end_row($row);
  }
}
if ($TABIX_MODE) {
  $tp->finish();
} else {
  $rpt->finish();
}
$rpt_fail->finish();

printf STDERR "missing genomic position: %d\n", $missing_genomic_pos;


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
