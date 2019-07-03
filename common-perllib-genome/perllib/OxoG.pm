package OxoG;
# get data for OxoG artifact, possibly for use with
# https://www.broadinstitute.org/cancer/cga/dtoxog
# (or standalone fix)
# MNE 1/2016

use strict;
use Exporter;

use Bio::DB::Sam;
use Bio::DB::Sam::Constants;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use Variant;

use constant CIGAR_MATCH => "M";
use constant CIGAR_INSERTION => "I";
use constant CIGAR_DELETION => "D";
use constant CIGAR_SOFT_CLIP => "S";
use constant CIGAR_SKIP => "N";
# FIX ME: Bio::DB::Sam::Constants
# see SAM specification.  Maybe these are in Bio modules anywhere?
# TO DO: utility code to fetch aligned blocks a-la picard?

@OxoG::ISA = qw(Configurable Exporter);
@OxoG::EXPORT_OK = qw();

use MethodMaker qw(
	is_possible_artifact

	count_ref_F1R2
	count_ref_F2R1
	count_alt_F2R1
	count_alt_F1R2
        count_alt_clonal

count_alt_fwd_read1
count_alt_fwd_read2
count_alt_rev_read1
count_alt_rev_read2

        require_properly_paired
        allow_optical_pcr_duplicates
        min_base_quality
        min_mapq
verbose
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->require_properly_paired(1);
  $self->allow_optical_pcr_duplicates(0);
  $self->min_base_quality(20);
  $self->min_mapq(0);
  $self->verbose(0);
  $self->configure(%options);
  return $self;
}

sub query {
  my ($self, %options) = @_;
  my $snv4 = get_hash_option(\%options, "-snv4");
  my $bams = get_hash_option(\%options, "-bams");
  my $v = new Variant();
  $v->import_snv4($snv4);
  die "must be SNV" unless $v->is_snv;

  my $ra = uc($v->reference_allele);
  my $va = uc($v->variant_allele);
  my $target_base = $v->start;

  my $possible_artifact = 0;
  if ($ra eq "C" and $va eq "A") {
    # artifact reads align to +
    $possible_artifact = 1;
  } elsif ($ra eq "G" and $va eq "T") {
    # artifact reads align to -
    $possible_artifact = 1;
  }
  $self->is_possible_artifact($possible_artifact);

  # for DtoxoG: https://www.broadinstitute.org/cancer/cga/dtoxog
  my $alt_F1R2 = 0;
  my $alt_F2R1 = 0;
  my $ref_F1R2 = 0;
  my $ref_F2R1 = 0;
  # NOTE: ref counts do NOT go through clonal filtering like alt counts,
  # not implemented as I'm not sure these will ever be used
  # (we're not planning to use dtoxog)

  my $alt_clonal = 0;
  my $alt_fwd_read1 = 0;
  my $alt_fwd_read2 = 0;
  my $alt_rev_read1 = 0;
  my $alt_rev_read2 = 0;

  my $VERBOSE = $self->verbose;
  my $allow_optical_pcr_duplicates = $self->allow_optical_pcr_duplicates();
  my $min_mapq = $self->min_mapq();
  my $min_base_quality = $self->min_base_quality();
  my $require_properly_paired = $self->require_properly_paired();

  foreach my $f_bam (@{$bams}) {
    # gather counts for each BAM

    my $bam = Bio::DB::Sam->new(
      "-bam" => $f_bam,
      "-expand_flags" => 1,
	);
    my $rnm = new ReferenceNameMapper();
    foreach my $chrom ($bam->seq_ids()) {
      # standardize reference names
      $rnm->add_name($chrom);
    }

    my $q_chr = $rnm->find_name($v->reference_name) || die;
    my @alignments_raw = $bam->get_features_by_location(
						    -seq_id => $q_chr,
						    -start  => $target_base,
						    -end    => $target_base
						   );

    #
    #  pre-filter reads to those passing minimum QC:
    #
    my @alignments;
    foreach my $aln (@alignments_raw) {
      my $usable = 1;
      $usable = 0 if $aln->qual() < $min_mapq;
      # minimum mapping quality

      unless ($allow_optical_pcr_duplicates) {
	my $optical_pcr_dup = $aln->get_tag_values("DUPLICATE");
	die "optical/pcr dup not defined" unless defined $optical_pcr_dup;
	$usable = 0 if $optical_pcr_dup;
      }

      if ($require_properly_paired) {
	$usable = 0 unless $aln->get_tag_values("MAP_PAIR");
      }
      push @alignments, $aln if $usable;
    }

    #
    #  first pass: check for clonal alignments
    #
    my %saw;
    my %saw_f;
    my %saw_r;
    my %clonal;
    foreach my $aln (@alignments) {
      my $seq_id = $aln->query->seq_id();
      $saw{$seq_id}++;
      if ($aln->strand > 0) {
	$saw_f{$seq_id}++;
      } elsif ($aln->strand < 0) {
	$saw_r{$seq_id}++;
      } else {
	die;
      }
    }
    foreach my $seq_id (sort keys %saw) {
      if ($saw{$seq_id} > 1) {
	if ($saw_f{$seq_id} and $saw_r{$seq_id}) {
	  $clonal{$seq_id} = 1;
	} else {
	  die "ERROR, sanity check fail: saw read $seq_id multiple times but not on opposite strands";
	  # need example to check/test:
	  # I think this might happen e.g. for some RNA mappers
	}
      }
    }

    #
    #  second pass: count
    #
    my %alt_ids;
    foreach my $aln (@alignments) {
      my $seq_id = $aln->query->seq_id();

      my $first_read = $aln->get_tag_values("FIRST_MATE");
      my $second_read = $aln->get_tag_values("SECOND_MATE");

      die unless $first_read or $second_read;
      die if $first_read and $second_read;

      my $strand_plus = $aln->strand > 0 ? 1 : 0;
      my $cigar_a = $aln->cigar_array();

      #    printf STDERR "id:%s usable:%d\n", $aln->query->seq_id(), $usable;

      #
      #  Walk through CIGAR and track.
      #
      my $ref_pos = $aln->start;
      # reference base alignment start (1-based)
      my $query_idx = $aln->query->start - 1;
      # query sequence index (0-based)

      my $query_dna = $aln->query->dna();
      my $query_qual = $aln->qscore();
      die "ERROR: seq/quality length mismatch" unless length($query_dna) == @{$query_qual};

      printf STDERR "pos:%d read:%s strand:%s CIGAR:%s\n",
      $ref_pos,
      $aln->query->seq_id,
      $aln->strand,
      $aln->cigar_str, $ref_pos if $VERBOSE;

      for (my $i = 0; $i < @{$cigar_a}; $i++) {
	my $c = $cigar_a->[$i];
	# walk through CIGAR
	my ($ctype, $clen) = @{$c};

	printf STDERR "  cigar:%s %d\n", $ctype, $clen if $VERBOSE;

	if ($ctype eq CIGAR_SOFT_CLIP) {
	  # ignore:
	  # - leading soft-clips are ignored in query start
	  # - we're not using trailing clips
	} elsif ($ctype eq CIGAR_MATCH) {
	  # match or mismatch: affects both query and reference
	  my ($q_nt, $q_qual);

	  for (my $j = 0; $j < $clen; $j++) {
	    $q_nt = uc(substr($query_dna, $query_idx, 1));
	    $q_qual = $query_qual->[$query_idx];
	    printf STDERR "  %d: %s %d\n", $ref_pos, $q_nt, $q_qual if $VERBOSE;

	    my $usable = 1;
	    $usable = 0 unless $ref_pos == $target_base;
	    $usable = 0 unless $q_qual >= $min_base_quality;

	    # in the basic counts, should we count reads not flagged as properly paired??

	    if ($ref_pos == $target_base and $q_qual >= $min_base_quality) {
	      #
	      #  good enough to track
	      #

	      my ($is_F1R2, $is_F2R1);
	      if ($first_read and $strand_plus) {
		#
		# first read of pair is on +, 2nd read is on -
		#
		$is_F1R2 = 1;
	      } elsif ($second_read and !$strand_plus) {
		$is_F1R2 = 1;
	      } elsif ($second_read and $strand_plus) {
		#
		# second read is on +, 1st read is on -
		#
		$is_F2R1 = 1;
	      } elsif ($first_read and !$strand_plus) {
		$is_F2R1 = 1;
	      } else {
		die sprintf "ignore %s %s strand:%s\n", $q_nt, $aln->query->seq_id, $aln->strand;
	      }

	      if ($usable) {
		if ($q_nt eq $ra) {
		  # read matches reference
		  $ref_F1R2++ if $is_F1R2;
		  $ref_F2R1++ if $is_F2R1;
		} elsif ($q_nt eq $va) {
		  # read matches variant
		  $alt_ids{$seq_id} = 1;

		  my $clonal = $clonal{$seq_id};

		  unless ($clonal) {
		    $alt_F1R2++ if $is_F1R2;
		    $alt_F2R1++ if $is_F2R1;

		    if ($strand_plus) {
		      if ($first_read) {
			$alt_fwd_read1++;
		      } else {
			$alt_fwd_read2++;
		      }
		    } else {
		      if ($first_read) {
			$alt_rev_read1++;
		      } else {
			$alt_rev_read2++;
		      }
		    }
		  }

		  my $tag;
		  if ($clonal) {
		    $tag = "clonal";
		  } elsif ($is_F1R2) {
		    $tag = "F1R2";
		  } elsif ($is_F2R1) {
		    $tag = "F2R1";
		  } else {
		    $tag = "ignore";
		  }

		  printf STDERR "variant:%s class:%s read:%s strand:%s\n", $q_nt, $tag, $seq_id, $aln->strand if $VERBOSE;

		}
	      }
	    }

	    $ref_pos++;
	    $query_idx++;
	  }
	} elsif ($ctype eq CIGAR_SKIP) {
	  # skip: affects the reference but not the query
	  $ref_pos += $clen;
	} elsif ($ctype eq CIGAR_INSERTION) {
	  # insertion: affects query but not reference
	  $query_idx += $clen;
	} elsif ($ctype eq CIGAR_DELETION) {
	  # deletion: affects reference but not query
	  $ref_pos += $clen;
	} else {
	  die "unhandled CIGAR entry $ctype";
	}
      }
    }

    foreach my $k (keys %clonal) {
      delete $clonal{$k} unless $alt_ids{$k};
    }
    $alt_clonal += scalar keys %clonal;

  }

  $self->count_alt_clonal($alt_clonal);


  if ($VERBOSE) {
    printf STDERR "alt_F2R1: %d\n", $alt_F2R1;
    printf STDERR "alt_F1R2: %d\n", $alt_F1R2;
    printf STDERR "ref_F1R2: %d\n", $ref_F1R2;
    printf STDERR "ref_F2R1: %d\n", $ref_F2R1;
    printf STDERR "alt_fwd_read1: %d\n", $alt_fwd_read1;
    printf STDERR "alt_fwd_read2: %d\n", $alt_fwd_read2;
    printf STDERR "alt_rev_read1: %d\n", $alt_rev_read1;
    printf STDERR "alt_rev_read2: %d\n", $alt_rev_read2;
  }

  $self->count_alt_F2R1($alt_F2R1);
  $self->count_alt_F1R2($alt_F1R2);
  $self->count_ref_F1R2($ref_F1R2);
  $self->count_ref_F2R1($ref_F2R1);
  $self->count_alt_fwd_read1($alt_fwd_read1);
  $self->count_alt_fwd_read2($alt_fwd_read2);
  $self->count_alt_rev_read1($alt_rev_read1);
  $self->count_alt_rev_read2($alt_rev_read2);

}



1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
