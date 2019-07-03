package IndelRemap;
# remap reads to reference sequence reflecting an indel for
# cleaner mapping
# MNE 11/2015
#
# - good for WGS
# - might epic fail in RNA!
#

use strict;

use Bio::SeqIO;
use Bio::DB::Sam;
use Bio::DB::Sam::Constants;
use Bio::Seq::Quality;

use Configurable;
use Exporter;
use TdtConfig;
use FAI;
use GenomeUtils qw(reverse_complement);
use BWA;
use TemporaryFileWrangler;
use MiscUtils qw(dump_die);

@IndelRemap::ISA = qw(Configurable Exporter);
@IndelRemap::EXPORT_OK = qw();

use constant CIGAR_MATCH => "M";
use constant CIGAR_INSERTION => "I";
use constant CIGAR_DELETION => "D";
use constant CIGAR_SOFT_CLIP => "S";
# see SAM specification.  Maybe these are in Bio modules anywhere?

use constant SAM_QUERY_BUFFER => 10;

my $REFERENCE_FLANK = 200;
my $REFERENCE_ID = "ref_indel";

use MethodMaker qw(
		    bam_realigned
                    modified_reference_fasta

		    variant
		    tfw

		    verbose

		    min_flanking_sequence
		    min_flanking_sequence_softclip
		    min_flanking_sequence_inverse_indel

                    min_quality
                    min_mapq
                    allow_optical_pcr_duplicates
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->tfw(new TemporaryFileWrangler());
#  $self->min_flanking_sequence(10);
  $self->min_flanking_sequence(7);
  $self->min_flanking_sequence_softclip(12);
  # be much stricter if adjacent CIGAR entry is a soft clip
  $self->min_flanking_sequence_inverse_indel(12);
  # also beware of non-supporting reads that introduce an indel
  # of the opposite type nearby.
  # e.g.
  # http://bamviewer-rt:8080/BAMViewer/aceview/splash?tumorname=/nfs_exports/genomes/1/projects/EXCAP/PanTARGET/BucketRaw/xmaSJAML/SJAML040538_D2.bam&ref=hg19&region=1&center=64624758&fullPath=true
  # deletion of 3 nt.  Realignment of non-supporting reads may align
  # perfectly at the new reference site, but introduce INSERTIONS of
  # the deleted sequence nearby!

  $self->min_quality(10);
  # min quality this should really be > 0, especially for small indels
  $self->min_mapq(0);
  # default 0 (in case BAM is SJ RNASeq)
  $self->allow_optical_pcr_duplicates(1);
  # default 1 (in case of amplicon or certain RNA sequencing projects)

  $self->verbose(0);
  $self->configure(%options);
  return $self;
}

sub remap {
  my ($self, %options) = @_;
  my $v = $options{"-variant"} || die;
  my $f_bam = $options{"-bam"} || die;
  $self->variant($v);
  my $tfw = $self->tfw();
#  $tfw->verbose(1);

  my $fasta = $options{"-fasta"};
  unless ($fasta) {
    my $genome = $options{"-genome"} || die;
    my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
    $fasta = $config_genome->{FASTA} || die "no fasta";
  }
  die unless $fasta;
  my $fai = new FAI("-fasta" => $fasta);

  my $new_reference;
  #
  #  build modified reference sequence reflecting the target deletion
  #  and create FASTA:
  #
  my $up = $fai->get_chunk(
			   "-id" => $v->reference_name(),
			   "-start" => $v->start() - $REFERENCE_FLANK,
			   "-length" => $REFERENCE_FLANK
			  );
  if ($v->is_deletion()) {
    my $down = $fai->get_chunk(
			       "-id" => $v->reference_name(),
			       "-start" => $v->end() + 1,
			       "-length" => $REFERENCE_FLANK
			    );
    $new_reference = $up . $down;
  } elsif ($v->is_insertion()) {
    #
    # the given position is the base number BEFORE the variant
    #
    my $up = $fai->get_chunk(
			     "-id" => $v->reference_name(),
			     "-start" => $v->start() - $REFERENCE_FLANK,
			     "-length" => $REFERENCE_FLANK + 1,
			  );
    my $down = $fai->get_chunk(
			       "-id" => $v->reference_name(),
			       "-start" => $v->start() + 1,
			       "-length" => $REFERENCE_FLANK
			    );

    printf STDERR "modified insertion ref: %s%s%s\n",
      lc($up), uc($v->variant_allele), lc($down) if $self->verbose();
    $new_reference = join "", $up, $v->variant_allele, $down;
  } elsif ($v->is_complex) {
    my $ra = $v->reference_allele();
    my $va = $v->variant_allele();
    foreach ($ra, $va) {
      die $_ if /\W/;
    }
    my $down = $fai->get_chunk(
			       "-id" => $v->reference_name(),
			       "-start" => $v->end() + 1,
			       "-length" => $REFERENCE_FLANK
			    );
    $new_reference = join "", lc($up), $va, lc($down);
  } else {
    die "not ins/del/complex";
  }

#  die $new_reference;

  my $fasta_file = $tfw->get_tempfile("-append" => ".temp_new_ref.fa", "-glob" => 1);
  my $fa = Bio::SeqIO->new(
			   "-format" => "fasta",
			   "-file"   => ">" . $fasta_file
			  );

  my $bs = Bio::Seq->new(
			 "-display_id" => $REFERENCE_ID,
			 "-seq" => $new_reference
			);

  $fa->write_seq($bs);
  $fa->close();
  $self->modified_reference_fasta($fasta_file);

  #
  #  find reads overlapping target site:
  #
  my $reference_event_length;
  if ($v->is_insertion()) {
    $reference_event_length = 1;
  } elsif ($v->is_deletion()) {
    $reference_event_length = $v->event_length();
  } elsif ($v->is_complex) {
    $reference_event_length = length $v->reference_allele();
  } else {
    die;
  }

  my $buf = SAM_QUERY_BUFFER + $reference_event_length;
  # for deletions, ensure we query overlapping the end site as well
  my $bam = Bio::DB::Sam->new(
			      "-bam" => $f_bam,
			     );
  my $qs = $v->start - SAM_QUERY_BUFFER;
  my $qe = $v->end + SAM_QUERY_BUFFER;

  my @alignments;
  if (my $aln = $options{"-alignments"}) {
    @alignments = @{$aln};
  } else {
    @alignments = $bam->get_features_by_location(
						  -seq_id => $v->reference_name,
						  -start  => $qs,
						  -end    => $qe,
						 );
  }

  #
  #  write as FASTQ, undoing BAM alignment strand adjustments:
  #
  my $fastq_file = $tfw->get_tempfile("-append" => ".temp_realign.fq", "-glob" => 1);

  my $fq = Bio::SeqIO->new(
			   "-format"    => "fastq",
			   "-file"      => ">" . $fastq_file
			  );
#  die join "\n", map {$_->query->seq_id} @alignments;

  foreach my $a (@alignments) {
    my $query_dna = $a->query->dna();
    my $query_qual = $a->qscore();
    # TO DO: bioperl fastq
    my $strand = $a->strand;
    if ($strand == 1) {
      # mapped natively to +
    } elsif ($strand == -1) {
      # undo BAM strand compensation (will be reapplied in realignment)
      $query_dna = reverse_complement($query_dna);
      $query_qual = [ reverse(@{$query_qual}) ];
    } else {
      die "strand";
    }
    my $bsq = Bio::Seq::Quality->new(
				     "-qual" => $query_qual,
				     "-id" => $a->query->seq_id,
				     "-seq" => $query_dna
				    );
    $fq->write_seq($bsq);
  }
  $fq->close();

  #
  #  realign
  #
  my $bwa = new BWA("-bwasw" => 1);
  my $new_bam = $bwa->bwa(
			  "-fastq" => $fastq_file,
			  "-reference" => $fasta_file
			 );
  # TO DO: cleanup
  $self->bam_realigned($new_bam);
}

sub get_indel_coverage {
  #
  # get count of reads supporting target indel.  In the remapped BAM,
  # these are reads with ungapped alignments covering the target site
  # (plus some minimum level of flanking sequence)
  #
  my ($self) = @_;
  my $f_bam = $self->bam_realigned();
  my $v = $self->variant();

  my $tstart = $REFERENCE_FLANK;
  # last base before event
  my $tend;
  # first base after event

  my $inverse_cigar;

  if ($v->is_insertion()) {
    # updated reference contains insertion
    $tend = $tstart + 1 + $v->event_length();
    $inverse_cigar = CIGAR_DELETION;
  } elsif ($v->is_deletion) {
    # deletion:
    # updated reference no longer contains deleted region
    $tend = $tstart + 1;
    $inverse_cigar = CIGAR_INSERTION;
  } elsif ($v->is_complex) {
    # complex:
    # updated reference no longer contains deleted region
    my $ra = $v->reference_allele;
    my $va = $v->variant_allele;
    if (length($ra) > length($va)) {
      # net deletion
      $tend = $tstart + length($va);
    } else {
      # net insertion
      die "TEST ME: net insertion";
    }
    $inverse_cigar = "";
  }

  my $tstart0 = $tstart;
  my $tend0 = $tend;
  # narrow site to check for sequence quality

  my $mf = $self->min_flanking_sequence - 1;
  # already have 1
  if ($mf > 0) {
    $tstart -= $mf;
    $tend += $mf;
    # adjust region for minimum flanking sequence
  }

  my $min_flanking_sequence_softclip = $self->min_flanking_sequence_softclip() - 1;
  # subtract 1 since generated coordinates we're using already contain
  # 1 flanking base
  my $min_flanking_sequence_inverse = $self->min_flanking_sequence_inverse_indel() - 1;

  my $bam = Bio::DB::Sam->new(
			      "-bam" => $f_bam,
			      "-expand_flags" => 1,
			     );
  my @alignments = $bam->get_features_by_location(
						  -seq_id => $REFERENCE_ID,
						  -start  => $REFERENCE_FLANK,
						  -end    => $REFERENCE_FLANK
						 );
  my $supporting_reads = 0;
  my $min_quality = $self->min_quality();
  my $min_mapq = $self->min_mapq();
  my $allow_optical_pcr_duplicates = $self->allow_optical_pcr_duplicates();

  foreach my $a (@alignments) {
    my $cigar_a = $a->cigar_array();

    my $usable = 1;
    $usable = 0 if $a->qual() < $min_mapq;
    # minimum mapping quality

    unless ($allow_optical_pcr_duplicates) {
      my $optical_pcr_dup = $a->get_tag_values("DUPLICATE");
      die "optical/pcr dup not defined" unless defined $optical_pcr_dup;
      $usable = 0 if $optical_pcr_dup;
    }

    if ($usable) {
      $usable = 0;

      my $ref_pos = $a->start;
      # reference base alignment start (1-based)
      my $query_idx = $a->query->start - 1;
      # query sequence index (0-based)

      my $query_dna = $a->query->dna();
      my $query_qual = $a->qscore();
      die "ERROR: seq/quality length mismatch" unless length($query_dna) == @{$query_qual};

      printf STDERR "pos:%d read:%s strand:%s CIGAR:%s\n",
	$ref_pos,
	  $a->query->seq_id,
	    $a->strand,
	      $a->cigar_str if $self->verbose;

      for (my $i = 0; $i < @{$cigar_a}; $i++) {
	my $c = $cigar_a->[$i];
	# walk through CIGAR
	my ($ctype, $clen) = @{$c};

#	printf STDERR "%s\n", join " ", $ctype, $clen;

	if ($ctype eq CIGAR_MATCH) {
	  # match or mismatch: affects both query and reference

	  my $ref_end = $ref_pos + $clen - 1;
	  if ($tstart >= $ref_pos and $tend <= $ref_end) {
            # target falls entirely within
#	    $usable = 1;
	    $usable = "$ref_pos $ref_end $tstart $tend";
#	    die $a->cigar_str . " $ref_pos $ref_end";

	    # flanking sequence checks for soft clips:
	    #
	    # when remapping to a tandem repeat region aligner
	    # may simply shift the mismatch to an adjacent soft clip.
	    # for example in chr1.64624751.TCT.- reads that DON'T
	    # support the deletion overlap the adjusted reference on
	    # remap because of the repeat region, and the deleted
	    # sequence is soft-clipped later.  This may be an argument
	    # for a longer flanking sequence, or maybe even only
	    # attempting to remap either reads with nearby soft clips
	    # and/or reads whose align start/end is near the site.
	    my $flank_left = $tstart0 - $ref_pos;
	    my $flank_right = $ref_end - $tend0;
	    if ($i > 0) {
	      my ($c2) = @{$cigar_a->[$i - 1]};
	      if ($c2 eq CIGAR_SOFT_CLIP) {
		# previous entry is a softclip
		if ($flank_left < $min_flanking_sequence_softclip) {
		  # soft clipping is too close
#		  dump_die($v, "left clip flank fail " . $a->query->seq_id);
		  $usable = "";
		}
	      }

	      if (
		($v->is_complex and ($c2 eq CIGAR_INSERTION or 
				     $c2 eq CIGAR_DELETION)) or
		$c2 eq $inverse_cigar
		# maybe more sophisticated handling for complex?
		  ) {
		if ($flank_left < $min_flanking_sequence_inverse) {
#		  dump_die($v, "invert leftflank fail " . $a->query->seq_id);
		  $usable = "";
		} else {
#		  dump_die($v, "invert leftflank OK " . $a->query->seq_id);
		}
	      }
	    }

	    if ($i < @{$cigar_a} - 1) {
	      my ($c2) = @{$cigar_a->[$i + 1]};
	      if ($c2 eq CIGAR_SOFT_CLIP) {
		if ($flank_right < $min_flanking_sequence_softclip) {
#		  die "right clip flank fail " . $a->query->seq_id;
		  $usable = "";
		} else {
#		  die "right clip flank OK " . $a->query->seq_id;
		}
	      }

	      if ($c2 eq $inverse_cigar) {
		if ($flank_right < $min_flanking_sequence_inverse) {
		  $usable = "";
#		  dump_die($v, "invert rightflank fail " . $a->query->seq_id);
		} else {
		  dump_die($v, "invert rightflank OK " . $a->query->seq_id);
		}
	      }

	    }


#	    die "get site quality";
	    my $ri = $ref_pos - 1;
	    for (my $i = 0; $i < $clen; $i++, $ri++) {
	      my $q = $query_qual->[$i];
#	      printf STDERR "%s\n", join " " , $i, $ri, substr($query_dna, $i, 1), $q;
	      if ($ri >= $tstart0 and $ri <= $tend0 and $q < $min_quality) {
		$usable = "";
	      }
	    }

	  }
	  $ref_pos += $clen;
	  $query_idx += $clen;
	} elsif ($ctype eq CIGAR_INSERTION) {
	  # insertion: affects query but not reference
	  $query_idx += $clen;
	} elsif ($ctype eq CIGAR_DELETION) {
	  # deletion: affects reference but not query
	  $ref_pos += $clen;
	} elsif ($ctype eq CIGAR_SOFT_CLIP) {
	  # ignore, since the alignment start position
	  # returned by the API is *net* of leading soft clipping.
	} else {
	  die "unhandled CIGAR entry $ctype";
	}
      }
    }
    $supporting_reads++ if $usable;
    printf STDERR "call for %s: remap_usable=%s\n", $a->query->seq_id, $usable if $self->verbose;
  }

  return $supporting_reads;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
