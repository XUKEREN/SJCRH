package MNVCounter;
# get count of reads supporting a MNV or complex variant from a BAM
# MNE 1/2016
#
# rewrite idea: find read position AFTER the event, then look
# upstream for full/partial evidence of variant allele.
# current code expects specific soft clipping behavior which
# may not be reliable.  HOWEVER: won't work for short reads that
# don't fully show variant allele (there is no anchoring after the event)
#

use strict;
use Exporter;
use Bio::DB::Sam;
use MNVCounter;
use ReferenceNameMapper;
use Variant;
use UniqueAlignmentID;

use Configurable;
use MiscUtils qw(dump_die get_hash_option);

@MNVCounter::ISA = qw(Configurable Exporter);
@MNVCounter::EXPORT_OK = qw();

use constant CIGAR_MATCH => "M";
use constant CIGAR_INSERTION => "I";
use constant CIGAR_DELETION => "D";
use constant CIGAR_SOFT_CLIP => "S";
use constant CIGAR_HARD_CLIP => "H";
use constant CIGAR_SKIP => "N";
# see SAM specification.  Maybe these are in Bio modules anywhere?

my $SAM_QUERY_BUFFER = 15;
# extend query to also capture soft clipped regions

my $COMPLEX_MAX_PRECEDING_INDEL_FORGIVE_LENGTH = 1;
# there may be spurious deletions BEFORE the site.
# http://bamviewer-rt:8080/BAMViewer/aceview/splash?tumorname=/rgs01/resgen/prod/tartan/runs/RNA_mapping/QP5rSpP9/output/SJAML040612_D1/finish/SJAML040612_D1.bam&ref=hg19&region=7&center=50444478&fullPath=true

my $MAX_EDGE_RESCUE_END_DISTANCE = 10;

use MethodMaker qw(
		    bam
		    bam_ref
		    rnm
		    supporting_reads
		    supporting_reads_indel
		    supporting_reads_soft_clip
verbose
allow_optical_pcr_duplicates
min_mapq
min_base_quality_coverage
coverage

complex_leading_deletion

partial_match_enable
partial_match_min_variant_allele_length
partial_match_min_variant_match_length

filter_name

		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->allow_optical_pcr_duplicates(1);
  $self->min_mapq(0);
  $self->min_base_quality_coverage(15);

  $self->partial_match_enable(1);
  if (0) {
    printf STDERR "WARNING: partial matching disabled\n";
    $self->partial_match_enable(0);
  }

  $self->partial_match_min_variant_allele_length(6);
  $self->partial_match_min_variant_match_length(4);
  # xiaotu 4/4/2016:
  # how about we use 6nt (of alternative allele) to define an allele
  # to be long?  If you use the left-most (or right-most) basepairs of
  # the alternative allele, can you do 3? Or 4?

  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self, %options) = @_;

  my $f_bam = $self->bam || die "-bam";

  my $bam = Bio::DB::Sam->new(
			      "-bam" => $f_bam,
			      "-expand_flags" => 1,
			     );
  my $rnm = new ReferenceNameMapper();
  foreach my $chrom ($bam->seq_ids()) {
    # standardize reference names
    $rnm->add_name($chrom);
  }
  $self->rnm($rnm);
  $self->bam_ref($bam);
}

sub query {
  my ($self, %options) = @_;
  my $snv4 = $options{"-snv4"} || die "-snv4";
  my ($chr_raw, $q_pos, $q_ra, $q_va) = split /\./, $snv4;
  $q_va = uc($q_va);
  my $q_chr = $self->rnm->find_name($chr_raw) || die "no BAM chrom for $chr_raw
";
  my $v = new Variant();
  $v->import_snv4($snv4);
  my $complex_mode;
  if ($v->is_mnv) {
  } elsif ($v->is_complex) {
    $complex_mode = 1;
  } else {
    dump_die($v, "variant is not a MNV/complex");
  }

  my $ref_pos_after = $v->end + 1;
  # first reference base AFTER the event 

  my $partial_match_enable = $self->partial_match_enable();
  my $partial_match_min_variant_allele_length = $self->partial_match_min_variant_allele_length();
  my $partial_match_min_variant_match_length = $self->partial_match_min_variant_match_length();

  my $complex_leading_deletion;

  my $ref_start = $v->start;
  my $ref_end = $v->end;

  my $tmp = $q_pos;

  my $qs = $q_pos - $SAM_QUERY_BUFFER;
  my $qe = $v->end + $SAM_QUERY_BUFFER;

  my $allow_optical_pcr_duplicates = $self->allow_optical_pcr_duplicates();
  my $min_mapq = $self->min_mapq();
  my $min_base_quality_coverage = $self->min_base_quality_coverage();
  my $VERBOSE = $self->verbose;

  my @alignments = $self->bam_ref->get_features_by_location(
						  -seq_id => $q_chr,
						  -start  => $qs,
						  -end    => $qe
						 );
  # TO DO: passed-in alignments

  my %total_coverage;
  $self->coverage(\%total_coverage);

  my $supporting_reads = 0;
  my $supporting_reads_indel = 0;
  my $supporting_reads_soft_clip = 0;

  if (0) {
    printf STDERR "DEBUG READ FILTER!!\n";
#    @alignments = grep {$_->query->seq_id() =~ /77517/} @alignments;
#    @alignments = grep {$_->query->seq_id() =~ /30450/} @alignments;
#    @alignments = grep {$_->query->seq_id() =~ /74172/} @alignments;
    @alignments = grep {$_->query->seq_id() =~ /80380/} @alignments;
  }

  my $filter_name = $self->filter_name();

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

    if ($filter_name) {
      next unless $a->query->seq_id() eq $filter_name;
    }

    my $downstream_chunk;
    my $has_deletion;
    my $has_insertion;
    my $has_soft_clip;
    my $seq_idx_after;
    # first reference base AFTER the event, in read index space
    my $cigar_match_i;
    my $cigar_match_qi;
    my $disqualified;

    printf STDERR "id:%s:%s check:%d\n", $a->query->seq_id(), $a->strand(), $usable if $VERBOSE;
    if ($usable) {
      #
      #  Walk through CIGAR and track.
      #
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
      $a->cigar_str, $ref_pos if $VERBOSE;

      my $last_del = 0;

      for (my $i = 0; $i < @{$cigar_a}; $i++) {
	my $c = $cigar_a->[$i];
	# walk through CIGAR
	my ($ctype, $clen) = @{$c};

	printf STDERR "  cigar:%s %d\n", $ctype, $clen if $VERBOSE;
	my $usable = 1;

	my $seq_ref = "NULL";
	my $seq_var = "NULL";

	if ($ctype eq CIGAR_MATCH or $ctype eq CIGAR_SOFT_CLIP) {
	  # match or mismatch: affects both query and reference
	  my $q;

	  my $is_leading_softclip;
	  if ($i == 0 and $ctype eq CIGAR_SOFT_CLIP) {

	    $ref_pos -= $clen;
	    $query_idx -= $clen;
	    # alignment start SKIPS any leading soft-clipped sequence.
	    # move alignment start back so we can include soft-clipped
	    # region in total coverage count without breaking sync.
	    # (this is fair since we will count MNVs in trailing clips)
	    $is_leading_softclip = 1;
	  }
	  
	  for (my $j = 0; $j < $clen; $j++) {
#	    printf STDERR "debug2 %s\n", join ",", $ref_pos, $ref_pos_after;
	    $seq_idx_after = $query_idx if $ref_pos == $ref_pos_after;
	    
	    $q = $query_qual->[$query_idx];
	    printf STDERR "  %d: %s %d\n", $ref_pos, substr($query_dna, $query_idx, 1), $q if $VERBOSE;

	    $total_coverage{$ref_pos}++ if $q >= $min_base_quality_coverage;

	    if ($ref_pos >= $q_pos and !defined($downstream_chunk)) {
	      my $dist = $ref_pos - $q_pos;
	      my $ok;

#	      $downstream_chunk = uc(substr($query_dna, $query_idx));
#	      printf STDERR "debug: %s\n", join ",", $a->query->seq_id, $ref_pos, $dist, $downstream_chunk, $q_pos;

	      if ($is_leading_softclip) {
		# soft-clipped region at start of the read,
		# special handling required
		$ok = 1;
	      } elsif ($dist == 0) {
		# site encountered with aligned block
		if ($ctype eq CIGAR_SOFT_CLIP and $j > 0) {
		  # trailing soft clipping is OK if the clipping starts
		  # immediately at the event site.  NOT OK if clipping
		  # starts before the event starts, because that portion
		  # should match.  e.g.
		  # http://bamviewer-rt:8080/BAMViewer/aceview/splash?normalname=/rgs01/resgen/prod/tartan/index/data/SJLIFE/SJLIFE/SJALL041247_G1/WHOLE_GENOME/bam/SJALL041247_G1.bam&center=162107943&fullPath=true&ref=hg20&region=chr1
		  # E00332:134:HK22NCCXX:4:2111:11627:32004 +
		  # event is ATTGACATTTGAACAATAAAATT>G, for this read
		  # there is a G at the expected position, but the soft clip
		  # earlier and the downstream sequence doesn't match
		  #
		  # ugh, there may be some nightmare cases where the
		  # read is soft clipped BUT IT STILL MATCHES THE REFERENCE!
		  # This filter would discard these reads.
		  $disqualified = 1;
		} else {
		  $ok = 1;
		}
	      } elsif ($dist == $last_del) {
		# a deletion touched the start of the MNV site.
		# The read may still support the MNV, but have a
		# spurious pattern of indels that interferes
		# with direct counting.
		$ok = 1;
	      } else {
		# unusable:
		# - read starts after the MNV site begins
		# - aligned sequence is too far away
		#   (large deletion or CIGAR N skip)
	      }

	      if ($ok) {
		if ($complex_mode) {
		  #
		  # complex indel
		  #
		  if ($ctype eq CIGAR_SOFT_CLIP) {
		    if ($is_leading_softclip) {
		      #
		      #  leading soft clip
		      #
		      if (length($q_ra) > length($q_va)) {
			# net deletion
#			my $idx = $query_idx + (length($q_ra) - length($q_va));
			if ($clen >= length($q_va)) {
			  # soft clip is long enough to contain entire
			  # variant allele
			  my $idx = ($query_idx + $clen) - length($q_va);
			  if ($idx < length($query_dna)) {
			    # clipping not always compatible, separate
			    # rescue later
			    $downstream_chunk = uc(substr($query_dna, $idx));
			    # expect inserted bases at END of leading soft clip
			  }
			}
#			printf STDERR "downstream: raw:%s cooked:%s\n", substr($query_dna, $query_idx), $downstream_chunk;
		      } else {
			# net insertion
			my $idx = $query_idx + length($q_ra) - length($q_va);
			$downstream_chunk = uc(substr($query_dna, $idx));
#			die "net insertion in leading softclip, test me $downstream_chunk";
		      }
		    } else {
		      # trailing soft clip, sometimes this aligns perfectly:
		      # http://bamviewer-rt:8080/BAMViewer/aceview/splash?ref=hg19&fullPath=FALSE&tumorname=projects/brT/SJBALL/SJBALL002231_D1-TARGET-10-PARGUZ-09A-01.bam&normalname=projects/brT/SJBALL/SJBALL002231_G1-TARGET-10-PARGUZ-10A-01.bam&region=chr19&center=34945243
		      $downstream_chunk = uc(substr($query_dna, $query_idx));
		    }
		  } else {
		    my $idx = $query_idx;
#		    die "hey now $last_del" if $last_del;
		    $complex_leading_deletion = 1 if $last_del;
		    if (0) {
		      # DISABLED: breaks this case:
		      # http://bamviewer-rt:8080/BAMViewer/aceview/splash?normalname=/nfs_exports/genomes/1/projects/EXCAP/TARGET_ALL/BucketRaw/SJALL/SJALL015916_G1-TARGET-10-PARDFB-10A-01.bam&tumorname=/nfs_exports/genomes/1/projects/EXCAP/TARGET_ALL/BucketRaw/SJBALL/SJBALL002179_D1-TARGET-10-PARDFB-09A-01.bam&center=68697869&fullPath=true&ref=hg19&region=1
		      #
		      # ...counterargument for the reason this code was
		      # added in the first place???
		      #

		      if ($last_del and $last_del <= $COMPLEX_MAX_PRECEDING_INDEL_FORGIVE_LENGTH) {
			$idx += $last_del;
			# compensate for interfering preceding deletions
		      }
		    }

		    $downstream_chunk = uc(substr($query_dna, $idx));
		    if ($ctype eq CIGAR_MATCH) {
		      $cigar_match_i = $i;
		      $cigar_match_qi = $idx;
		    }

		  }
		} else {
		  #
		  # MNV
		  #
		  $downstream_chunk = uc(substr($query_dna, $query_idx));
		}
		printf STDERR "  downstream: %s\n", join ",", $a->query->seq_id, $ref_pos, $dist, $downstream_chunk, $q_pos if $VERBOSE;
	      } else {
		$downstream_chunk = "";
	      }

	      if ($ctype eq CIGAR_SOFT_CLIP and
		  $ref_pos >= $ref_start and
		  $ref_pos <= $ref_end) {
		$has_soft_clip = $a;
	      }
	    }
	    $ref_pos++;
	    $query_idx++;
	  }
	} elsif ($ctype eq CIGAR_SKIP) {
	  # skip: affects the reference but not the query
	  $ref_pos += $clen;
	  $last_del = 0;
	} elsif ($ctype eq CIGAR_INSERTION) {
	  # insertion: affects query but not reference

	  if ($complex_mode and 
	      $ref_pos >= $q_pos and
	      !defined($downstream_chunk)) {
	    my $dist = $ref_pos - $q_pos;
	    my $ok = $dist == 0;

	    if ($ok) {
	      $downstream_chunk = uc(substr($query_dna, $query_idx));
	    }
	  }

	  $query_idx += $clen;
	  $last_del = 0;
	  $has_insertion = $a if $ref_pos >= $ref_start and
	      $ref_pos <= $ref_end;
	} elsif ($ctype eq CIGAR_DELETION) {
	  # deletion: affects reference but not query
	  for (my $j = 0; $j < $clen; $j++) {
	    $has_deletion = $a if $ref_pos >= $ref_start and $ref_pos <= $ref_end;
	    $total_coverage{$ref_pos++}++;
	  }
	  $last_del = $clen;
	} elsif ($ctype eq CIGAR_HARD_CLIP) {
	  # hard clipping: sequence is not present in read
	  for (my $j = 0; $j < $clen; $j++) {
#	    printf STDERR "hard clip pos:%d read:%s strand:%s CIGAR:%s qidx:%d down:%s\n", $ref_pos, $a->query->seq_id, $a->strand, $a->cigar_str, $query_idx, substr($query_dna, $query_idx);
	    $ref_pos++;
	  }
	} else {
	  die "unhandled CIGAR entry $ctype";
	}
      }
    }

    #
    #  main checks requiring evidence of entire variant allele:
    #
    my $usable2;
    if ($downstream_chunk) {
      $usable2 = index($downstream_chunk, $q_va) == 0 ? 1 : 0;
#      printf STDERR "downstream for %s = %s, usable2=%d\n", $a->query->seq_id, $downstream_chunk, $usable2;
    }

#    printf STDERR "debug %s\n", join " ", $a->query->seq_id, $q_va, ($seq_idx_after || "n/a");

    if ($usable and not($usable2) and not($disqualified) and $seq_idx_after) {
      # if the base number AFTER the event is known, separate check looking
      # upstream from that site.  This may rescue some complex variants.
      my $upstream = substr($a->query->dna(), 0, $seq_idx_after);

      if (length($upstream) >= length($q_va)) {
	my $chunk = substr($upstream, - length($q_va));
	if ($chunk eq $q_va) {
	  #	  dump_die($v, "rescued! " . $a->query->seq_id);
	  $usable2 = 1;
	  printf STDERR "full upstream rescue for %s\n", join " ", $a->query->seq_id, ($has_soft_clip ? "soft" : "not_soft");
	}
      }
    }


    if ($usable2 and defined($cigar_match_i) and
	length($q_ra) != length($q_va)) {
      # for complex indels where anchor was found in a cigar M,
      # require an adjacent CIGAR indel, otherwise a SNV matching
      # the variant allele can trigger a false positive.

      # link2bambino.pl -url 'http://bamviewer-rt:8080/BAMViewer/aceview/splash?tumorname=/rgs01/resgen/prod/tartan/index/data/SCMC/SCMC/SJALL040462_O2/VALCAP/bam/SJALL040462_O2.bam&ref=hg19&region=chr11&center=36619739&fullPath=true' -extra "-restrict ST-E00209:293:HFMF3ALXX:3:2103:16092:6689 -no-db"
      # 10.36619743.CC.T
      # this read has a T and so superficially matches the downstream
      # sequence, however it doesn't contain a deletion!
      #
      # this is more likely the fewer bases there are in the variant
      # allele, as a SNV might match.
      my $ok;
      foreach my $dir (-1, 1) {
	my $ci = $cigar_match_i + $dir;
	if ($ci >=0 and $ci < @{$cigar_a}) {
	  # check for an adjacent event
	  my $t = $cigar_a->[$ci]->[0];
	  $ok = 1 if $t eq CIGAR_INSERTION or $t eq CIGAR_DELETION;
	  # typically we'd expect to find an insertion or a deletion
	  # depending on whether the variant allele is longer or
	  # shorter than the reference allele.  However some rare mapping
	  # cases stymie this, e.g. 
	  # mchack -bam /rgs01/resgen/prod/tartan/index/data/SCMC/SCMC/SJALL018372_O3/VALCAP/bam/SJALL018372_O3.bam -snv4 chr10.104850749.G.ATGAGGA
	  # here we'd expect an insertion, however the read
	  # ST-E00244:210:H57MYALXX:4:2104:25083:28312 shows a D
	  # of the target base and a subsequent cigar M, 3 bases of
	  # which don't match the reference (near read end).
	}
      }

      unless ($ok) {
	# no flanking indel event
	if ($downstream_chunk) {
	  # downstream rescue, e.g.
	  # http://bamviewer-rt:8080/BAMViewer/aceview/splash?normalname=/rgs01/resgen/prod/tartan/index/data/SJLIFE/SJLIFE/SJALL041247_G1/WHOLE_GENOME/bam/SJALL041247_G1.bam&center=162107943&fullPath=true&ref=hg20&region=chr1
	# E00332:134:HK22NCCXX:4:2111:30107:32514 on -
	  my $dsl = length($downstream_chunk);
	  if ($dsl <= length($q_va)) {
	    # variant evidence is at very end of read: OK
	    $ok = 1;
	  } elsif ($dsl <= $MAX_EDGE_RESCUE_END_DISTANCE) {
	    # might be a small mapped portion that shows evidence of the event,
	    # (typically with mismatches), similar to soft clipping
	    $ok = 1;
	  }
	} elsif ($seq_idx_after) {
	  # upstream rescue
	  my $upstream = substr($a->query->dna(), 0, $seq_idx_after);
	  if (length($upstream) <= length($q_va)) {
	    printf STDERR "TEST ME: upstream read-end rescue\n";
	    $ok = 1;
	  }
	}
      }


      if ($ok) {
#	printf STDERR "%s usable\n", $a->query->seq_id;
      } else {
	printf STDERR "NOTE: %s strand:%s disqualified for complex %s in %s, CIGAR qi=%d\n", $a->query->seq_id, $a->strand, $snv4, $self->bam, $cigar_match_qi;
	$usable2 = 0;


      }
    }



    #
    #  optional checking for partial matches to non-variant alleles:
    #
    if (not($usable2) and
	$partial_match_enable) {

      if ($downstream_chunk and
	  length($q_va) > length($downstream_chunk) and
	  # only allow if downstream chunk is shorter than variant sequence!
	  length($q_va) >= $partial_match_min_variant_allele_length and
	  length($downstream_chunk) >= $partial_match_min_variant_match_length) {
	# rescue reads with long variant alleles where the read sequence
	# is not long enough to show the full variant sequence,
	# but is long enough to show a significant portion of it.
	#
	# this section deals with TRAILING matches, i.e. we know
	# what downstream sequence to look for, the read is just too
	# short to contain the whole thing
	my $search = substr($q_va, 0, $partial_match_min_variant_match_length);
	if (index($downstream_chunk, $search) == 0) {
	  $usable2 = 1;
	  printf STDERR "short rescue downstream for %s\n", join " ", $a->query->seq_id, $downstream_chunk;
	}
      }

      if ($seq_idx_after and
	  length($q_va) >= $partial_match_min_variant_allele_length) {
	# partial variant allele rescue for LEADING matches.
	# here we know the base number AFTER the event, so look upstream
	# for a partial match.
	my $upstream = substr($a->query->dna(), 0, $seq_idx_after);

	if (length($upstream) >= $partial_match_min_variant_match_length and
	    length($q_va) > length($upstream)
	    # only if upstream sequence is not long enough to observe variant
	    ) {
	  my $search = substr($q_va, - $partial_match_min_variant_match_length);
	  my $chunk = substr($upstream, - $partial_match_min_variant_match_length);
	  if ($chunk eq $search) {
	    $usable2 = 1;
	    printf STDERR "short rescue upstream for %s\n", join " ", $a->query->seq_id, $chunk;
	  }
	}
      }
    }

    printf STDERR "debug: %s\n", join " ", (defined $usable2 ? $usable2 : "no_chunk"), "has_softclip:" . ($has_soft_clip ? 1 : 0), $a->query->seq_id, $a->query->strand, $downstream_chunk if $VERBOSE;
    if ($usable2) {
      die "shouldn't happen: unusable read" unless $usable;
      $supporting_reads++;
      if ($has_insertion or $has_deletion) {
	$supporting_reads_indel++;
      } elsif ($has_soft_clip) {
	$supporting_reads_soft_clip++;
      }
    }
  }

  $self->complex_leading_deletion($complex_leading_deletion);
  $self->supporting_reads($supporting_reads);
  $self->supporting_reads_indel($supporting_reads_indel);
  $self->supporting_reads_soft_clip($supporting_reads_soft_clip);

  return $supporting_reads;
}

sub get_coverage {
  my ($self, %options) = @_;
  my $base = $options{"-base"} || die "-base";
  return $self->coverage->{$base};
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
