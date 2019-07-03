package dbSNP_DBS;
# fast search interface for dbSNP flatfiles
# MNE 2/2014

use strict;
use Carp qw(confess);

use Configurable;
use Exporter;

use DelimitedBinarySearch;
use GenomeUtils qw(cook_chromosome_name reverse_complement);
use MiscUtils qw(dump_die);

use constant CLASS_SNV => "single";
use constant CLASS_MNV => "mnp";

@dbSNP_DBS::ISA = qw(Configurable Exporter);
@dbSNP_DBS::EXPORT_OK = qw();

my @UCSC_SNP_HEADERS = (
			"bin",
			"chrom",
			"chromStart",
			"chromEnd",
			"name",
			"score",
			"strand",
			"refNCBI",
			"refUCSC",
			"observed",
			"molType",
			"class",
			"valid",
			"avHet",
			"avHetSE",
			"func",
			"locType",
			"weight",
			"exceptions",
			"submitterCount",
			"submitters",
			"alleleFreqCount",
			"alleles",
			"alleleNs",
			"alleleFreqs",
			"bitfields",
		       );
# HACK: exact fields will depend on version parsed!
# this is dangerous if fields are reordered, etc.
# maybe parse associated .sql file??
# e.g.
# http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/snp137.sql

use MethodMaker qw(
	dbsnp
        dbs
        snvs_only
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->snvs_only(1);
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self, %options) = @_;
  my $dbsnp_file = $self->dbsnp || die "-dbsnp";

  my $dbs = new DelimitedBinarySearch(
    "-file" => $dbsnp_file,
    "-is_headered" => 0,
      );
  $dbs->max_line_length(13000);
  # rs80338793: large deletion

  $self->dbs($dbs);
}

sub find {
  my ($self, %options) = @_;
  my @filtered;

  if (my $r = $options{"-sj"}) {
    $options{"-chrom"} = $r->{Chr} || die;
    $options{"-base-number"} = $r->{WU_HG19_Pos} || $r->{WU_HG18_Pos} || confess "can't find base number";
    $options{"-base-reference"} = $r->{ReferenceAllele} || die "need ReferenceAllele";
    $options{"-base-variant"} = $r->{MutantAllele} || die "need MutantAllele";
  }

  my $chrom = cook_chromosome_name(($options{"-chrom"} || die "-chrom"),
				   "-ucsc" => 1);
  # UCSC-format chrom name
  my $bn = $options{"-base-number"} || die "-base-number";

  my $dbs = $self->dbs() || die;
  $dbs->hashify_headers(\@UCSC_SNP_HEADERS);
  # hack: will likely 'splode for later version

  my $hits = $dbs->find(
		     "-comparators" => [
					{
					"column_number" => 2,
					"type" => "string",
					"value" => $chrom
					 # chromosome
					},

					{
					"type" => "number",
#					"column_number" => 4,
#					"value" => $bn
					"column_number" => 3,
					"value" => $bn - 1,
					# base #

# critically import to search by START field (interbase / 0-based)
# rather than end.  File is sorted by start, NOT END!  if we search by
# end we may miss some entries, e.g.  rs80359151 @ chr13.32953932

					},
				       ],

		       "-value" => $bn,
		       "-verbose" => $options{"-verbose"}
		      );
  my $snvs_only = $self->snvs_only();
  my $ra = $options{"-base-reference"};
  my $va = $options{"-base-variant"};

  if ($snvs_only) {
    # SNV mode
    my $bad;
    $bad = 1 if ($ra and $ra eq "-") or ($va and $va eq "-");
    # user variant is an indel, nothing will match
    return [] if $bad;
  }

  my $query_is_mnv;
  if ($ra) {
    $query_is_mnv = length($ra) > 1;
  }

  foreach my $hit (@{$hits}) {
#    dump_die($hit, "checking hit", 1);
    my $usable = 1;

    my $db_ref = $hit->{refUCSC};
    my $class = $hit->{class};

    if ($snvs_only) {
      if ($ra) {
	if ($query_is_mnv) {
	  # user query type must match database entry type
	  $usable = 0 unless $class eq CLASS_MNV;
	} else {
	  $usable = 0 unless $class eq CLASS_SNV;
	}
      }

      $usable = 0 if $hit->{chromStart} == $hit->{chromEnd};
      # rs41293467: zero span
    }

    if ($usable and $ra) {
      if ($class eq CLASS_SNV) {
	# for SNVs, expect the reference batch to match perfectly
	if (length($db_ref) != length($ra)) {
	  # e.g. chr2.233243796.T.C hits rs1130338, annotated
	  # as an SNV but with multi-nucleotide reference and variant bases
	  # (feh)
	  $usable = 0;
	} elsif ($db_ref ne uc($ra)) {
	  dump_die($hit, "ERROR: reference base mismatch, user=$ra db=$db_ref", 1); 
	  $usable = 0;
	}
      } elsif ($class eq CLASS_MNV) {
	if (uc($ra) eq $db_ref) {
	  # good: db entry has reference sequence entry that's the
	  # same length as the MNV sequence, and it matches.
	  # (e.g. rs34714427)
	} else {
#	  dump_die($hit, "MNV reference mismatch, user=$ra db=$db_ref");
	  $usable = 0;
	  # some of these might salvagable if the desired allele
	  # appears in the allele set even if the given reference
	  # base is only 1 nt.
	}
      } else {
	dump_die($hit, "unhandled class $class, user=$ra db=$db_ref");
      }
    }

    if ($usable and $va) {
      # user specified variant allele
      my $vb = $self->get_variant_bases($hit);
      if (grep {$va eq $_} @{$vb}) {
	# ok
      } else {
#	dump_die($options{"-sj"}, "SJ query:", 1);
#	dump_die($hit, "dbSNP mismatch!");
	# different variant at site
	$usable = 0;
      }
    }

#    dump_die($hit, "debug $usable", 1);

    push @filtered, $hit if $usable;
  }

  return \@filtered;
}

sub get_variant_bases {
  # return genomic variant bases for a SNV
  # (correcting for strand)
  my ($self, $hit) = @_;
  my $observed = $hit->{observed} || die;
  my $strand = $hit->{strand} || die;
  my @entries = split /\//, $observed;
  my $reference_base = uc($hit->{refUCSC} || die);
  my $saw_reference;
  my @variant_bases;
  foreach my $entry (@entries) {
    if ($strand eq "-") {
      $entry = reverse_complement($entry);
    } elsif ($strand ne "+") {
      die;
    }

    $entry = uc($entry);
    if ($entry eq $reference_base) {
      $saw_reference = 1;
    } else {
      push @variant_bases, $entry;
    }
  }

#  rs144571919

  dump_die($hit, sprintf("WARNING: didn't find reference %s in %s!",
			 $reference_base, $hit->{name}), 1)
      unless $saw_reference;
  # rs144571919:
  # reference is a C, observed bases are A/C, strand is +...WTF?
  return \@variant_bases;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
