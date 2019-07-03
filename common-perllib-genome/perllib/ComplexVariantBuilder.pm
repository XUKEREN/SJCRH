package ComplexVariantBuilder;
# scan a BAM interval and construct putative complex variants
# MNE 1/2016

use strict;
use Exporter;
use Bio::DB::Sam;
use ComplexVariantBuilder;
use ReferenceNameMapper;
use Variant;
use Counter;

use Configurable;
use MiscUtils qw(dump_die get_hash_option);

@ComplexVariantBuilder::ISA = qw(Configurable Exporter);
@ComplexVariantBuilder::EXPORT_OK = qw();

use constant CIGAR_MATCH => "M";
use constant CIGAR_INSERTION => "I";
use constant CIGAR_DELETION => "D";
use constant CIGAR_SOFT_CLIP => "S";
use constant CIGAR_HARD_CLIP => "H";
use constant CIGAR_SKIP => "N";
# see SAM specification.  Maybe these are in Bio modules anywhere?

my $SAM_QUERY_BUFFER = 0;
# based on formally-aligned regions only so not needed

my $MAX_NEIGHBOR_DISTANCE = 7;
# 6 at least
#my $MAX_NEIGHBOR_DISTANCE = 18;
#
# this example requires 18:
#
# http://bamviewer-rt:8080/BAMViewer/aceview/splash?tumorname=projects/WHOLEGENOME/PanTARGET/BucketRaw/SJAMLIF/SJAML040688_D1.bam&normalname=projects/WHOLEGENOME/PanTARGET/BucketRaw/SJAMLIF/SJAML040688_G1.bam&ref=hg19&region=17&center=42288184&fullPath=false
# 
my $MAX_NEIGHBOR_DISTANCE_SNV = 5;
# be stricter about SNV-to-SNV joining to avoid constructing large MNVs
# from nearby SNPs, e.g.
# 
# http://bamviewer-rt:8080/BAMViewer/aceview/splash?tumorname=/nfs_exports/genomes/1/projects/VALCAP/TARGET_ALL/BucketRaw/SJCOGALL/SJCOGALL010887_R1-PANKAK_R.bam&normalname=/nfs_exports/genomes/1/projects/VALCAP/TARGET_ALL/BucketRaw/SJCOGALL/SJCOGALL010887_G1-PANKAK_G.bam&ref=hg19&region=chr16&center=55169853&fullPath=true


my $MIN_FLANKING_SEQUENCE = 12;
# complex variants may look like MNVs near read edges

use MethodMaker qw(
		    bam
		    bam_ref
		    rnm
                    fai
nt_cache_single
chr_raw
verbose
min_good_reads
called
called_simple
raw_alignment_count
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->verbose(0);
  $self->min_good_reads(2);
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

sub reset {
  my ($self) = @_;
  $self->nt_cache_single({});
}

sub query {
  my ($self, %options) = @_;
  my ($chr_raw, $qs, $qe);
  my $read_regexp = $options{"-read-regexp"};

  $self->reset();

  my $gate_pos;

  if (my $snv4 = $options{"-snv4"}) {
    my ($q_pos, $q_ra, $q_va);
    ($chr_raw, $q_pos, $q_ra, $q_va) = split /\./, $snv4;
    my $v = new Variant();
    $v->import_snv4($snv4);
    $qs = $v->start;
    $qe = $v->end;
    $gate_pos = $q_pos;
  } elsif ($chr_raw = $options{"-chr"}) {
    $qs = $qe = $gate_pos = $options{"-position"} || die "-chr requires -position";
    if ($options{"-interval-mode"}) {
      $qe = $options{"-end"} || die "-end";
      $gate_pos = undef;
    }
  } else {
    die "specify -snv4 or -chr and -position";
  }

  $self->chr_raw($chr_raw);
  my $q_chr = $self->rnm->find_name($chr_raw) || die "no BAM chrom for $chr_raw
";

  my $allow_optical_pcr_duplicates = 1;
  my $min_mapq = 0;
  my $min_base_quality = 20;
#  my $min_base_quality = 10;
  # base on high-quality events only

  my $VERBOSE = $self->verbose();
  my %good;
  my %simple;

  my @alignments = $self->bam_ref->get_features_by_location(
						  -seq_id => $q_chr,
						  -start  => $qs,
						  -end    => $qe
						 );
  $self->raw_alignment_count(scalar @alignments);
  printf STDERR "alignments: %d\n", scalar @alignments if $VERBOSE;

  my $c = new Counter(\@alignments, "-mod" => 250);
  foreach my $aln (@alignments) {
    $c->next();
    if ($read_regexp) {
      printf STDERR "DEBUG: read filter\n";
#      next unless $aln->query->seq_id() =~ /51626/ and $aln->strand == 1;
#      next unless $aln->query->seq_id() =~ /93045/;
      next unless $aln->query->seq_id() =~ /$read_regexp/;
    }

    my $cigar_a = $aln->cigar_array();
    my $usable = 1;
    $usable = 0 if $aln->qual() < $min_mapq;
    # minimum mapping quality

    unless ($allow_optical_pcr_duplicates) {
      my $optical_pcr_dup = $aln->get_tag_values("DUPLICATE");
      die "optical/pcr dup not defined" unless defined $optical_pcr_dup;
      $usable = 0 if $optical_pcr_dup;
    }

#    printf STDERR "id:%s usable:%d\n", $aln->query->seq_id(), $usable;
    my %events;
    my %skip_sites;

    if ($usable) {
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

      printf STDERR "pos:%d read:%s strand:%s length:%d qi:%d CIGAR:%s read:%s\n",
	$ref_pos,
	  $aln->query->seq_id,
	    $aln->strand,
	      length($query_dna),
		$query_idx,
		  $aln->cigar_str,
		    $query_dna,
		      if $VERBOSE;

      for (my $i = 0; $i < @{$cigar_a}; $i++) {
	my $c = $cigar_a->[$i];
	# walk through CIGAR
	my ($ctype, $clen) = @{$c};

	printf STDERR "  cigar:%s %d\n", $ctype, $clen if $VERBOSE;

	if ($ctype eq CIGAR_SOFT_CLIP) {
	  # ignore:
	  # - query alignment start is AFTER leading clips
	  # - we're not using trailing clips
	} elsif ($ctype eq CIGAR_HARD_CLIP) {
	  $query_idx -= $clen if $i == 0;
	  # leading hard clip: query pointer incorrect in these cases,
	  # seems to imply hard-clipped sequence is present in read, it isn't
	} elsif ($ctype eq CIGAR_MATCH) {
	  # match or mismatch: affects both query and reference
	  my ($q_nt, $q_qual);

	  for (my $j = 0; $j < $clen; $j++) {
	    $q_nt = uc(substr($query_dna, $query_idx, 1));
	    die sprintf "q index %d beyond array of %d", $query_idx, scalar(@{$query_qual}) if $query_idx >= @{$query_qual};
	    $q_qual = $query_qual->[$query_idx];
	    printf STDERR "  %d: %s %d\n", $ref_pos, $q_nt, $q_qual if $VERBOSE;
	    my $ref_nt = $self->get_ref_single($ref_pos);
	    if ($q_nt ne $ref_nt) {
	      # mismatch
	      my %event;
	      $event{qc} = $q_qual >= $min_base_quality ? 1 : 0;
	      $event{type} = CIGAR_MATCH();
	      $event{pos} = $ref_pos;
	      $event{length} = 1;
	      $event{sequence} = $q_nt;
	      $events{$ref_pos}{CIGAR_MATCH()} = \%event;
	    }

	    $ref_pos++;
	    $query_idx++;
	  }
	} elsif ($ctype eq CIGAR_SKIP) {
	  # skip: affects the reference but not the query
	  $skip_sites{$ref_pos} = 1;
	  $ref_pos += $clen;
	  $skip_sites{$ref_pos} = 1;
	} elsif ($ctype eq CIGAR_INSERTION) {
	  # insertion: affects query but not reference
	  my $usable = 1;
	  my ($q_nt, $q_qual);
	  my $inserted = "";
	  my $qc = 1;
	  for (my $j = 0; $j < $clen; $j++) {
	    $q_nt = uc(substr($query_dna, $query_idx, 1));
	    $q_qual = $query_qual->[$query_idx];
	    printf STDERR "  %d: %s %d (ins)\n", $ref_pos, $q_nt, $q_qual if $VERBOSE;

	    $qc = 0 if $q_qual < $min_base_quality;
	    $inserted .= $q_nt;
	    $query_idx++;
	  }
	  
	  my %event;
	  $event{type} = CIGAR_INSERTION;
	  $event{pos} = $ref_pos;
	  $event{length} = $clen;
	  $event{sequence} = $inserted;
	  $event{qc} = $qc;
	  $events{$ref_pos}{CIGAR_INSERTION()} = \%event;
	} elsif ($ctype eq CIGAR_DELETION) {
	  # deletion: affects reference but not query
	  my %event;
	  $event{qc} = 1;
	  $event{type} = CIGAR_DELETION;
	  $event{pos} = $ref_pos;
	  $event{length} = $clen;
	  $events{$ref_pos}{CIGAR_DELETION()} = \%event;

	  $ref_pos += $clen;
	} else {
	  die "unhandled CIGAR entry $ctype";
	}
      }
    }

    #
    # finished parsing CIGAR
    #

    my @starts = sort {$a <=> $b} keys %events;

    my @ranges;

    for (my $si = 0; $si < @starts; $si++) {
      my $last_pos = $starts[$si];
      my $ni = $si + 1;
      while (1) {
	last if $ni >= @starts;
	my $np = $starts[$ni];
	my $dist = $np - $last_pos;

	my $ok = 0;

	if ($events{$last_pos}{CIGAR_MATCH()} and
	    $events{$np}{CIGAR_MATCH()}) {
	  if ($dist <= $MAX_NEIGHBOR_DISTANCE_SNV) {
	    $ok = 1;
	  } else {
#	    printf STDERR "reject bad merge from $last_pos to $np dist=$dist\n";
	  }
	} else {
	  # mixed event
	  $ok = 1 if $dist <= $MAX_NEIGHBOR_DISTANCE;
	}

	if (not($ok) and $events{$last_pos}{CIGAR_DELETION()}) {
	  my $pos2 = $last_pos + $events{$last_pos}{CIGAR_DELETION()}{length} - 1;
	  $dist = $np - $pos2;
	  $ok = 1 if $dist <= $MAX_NEIGHBOR_DISTANCE;
	}

	if ($ok) {
	  # close enough to combine, continue
	  $ni++;
	  $last_pos = $np;
	} else {
	  last;
	}
      }

      if ($last_pos >= $starts[$si]) {
	# usable set (even if only 1, can be 2 events at site)
	push @ranges, [ $starts[$si], $last_pos ];
	#	  printf STDERR "new range: %d-%d\n", $starts[$si], $last_pos;
      }
    }

    my $ranges_final = prune_ranges(\@ranges);
    # prune ranges that are subsets of larger ranges
    #      printf STDERR "ranges:%d final:%d %s\n", scalar(@ranges), scalar @{$ranges_final}, join ", ", map {$_->[0] . "-" . $_->[1]} @{$ranges_final};

    my $aln_start = $aln->start;
    my $aln_end = $aln->end;

    foreach my $range (@{$ranges_final}) {
      # combinable
      my ($start_pos, $end_pos) = @{$range};

      my $seq_ref = "";
      my $seq_var = "";
      my %deleted;

      my $qc_problem = 0;

      for (my $ref_pos = $start_pos; $ref_pos <= $end_pos; $ref_pos++) {
	$seq_ref .= $self->get_ref_single($ref_pos);

	my %e = map {$_, 1} keys %{$events{$ref_pos}};
	if (scalar keys %e > 1) {
	  die unless $e{CIGAR_MATCH()} and $e{CIGAR_INSERTION()};
	  # test me
	}

	my ($first, @rest);
	foreach my $type (keys %{$events{$ref_pos}}) {
	  if ($type eq CIGAR_INSERTION) {
	    $first = $type;
	  } else {
	    push @rest, $type;
	  }
	}
	my @all;
	push @all, $first if $first;
	push @all, @rest;

	my $has_snv;

	foreach my $type (@all) {
	  # process events at this site in mapping order:
	  # insertions processed first as they occur BEFORE this base
	  my $event = $events{$ref_pos}{$type};
	  my $qc = $event->{qc};
	  die "qc not defined" unless defined $qc;
	  $qc_problem = "poor_quality_inserted_sequence" unless $qc;
	  # poor-quality inserted sequence

	  my $dist_from_as = $ref_pos - $aln_start;
	  my $dist_from_ae = $aln_end - $ref_pos;
	  if ($dist_from_as < $MIN_FLANKING_SEQUENCE) {
	    $qc_problem = "too_close_to_align_start";
	  } elsif ($dist_from_ae < $MIN_FLANKING_SEQUENCE) {
	    $qc_problem = "too_close_to_align_end";
	  }

	  foreach my $ss (keys %skip_sites) {
	    if (abs($ref_pos - $ss) < $MIN_FLANKING_SEQUENCE) {
	      $qc_problem = "too_close_to_skip";
	    }
	  }

	  my $add_reference;

	  if ($type eq CIGAR_DELETION) {
	    my $end = $ref_pos + $event->{length} - 1;
	    for (my $i = $ref_pos; $i <= $end; $i++) {
	      $deleted{$i} = 1;
	      $end_pos = $i if $i > $end_pos;
	      # if the final event is a deletion, ensure we include
	      # the reference including the entire event
	    }
	  } elsif ($type eq CIGAR_INSERTION) {
	    $seq_var .= $event->{sequence} || die;
	  } elsif ($type eq CIGAR_MATCH) {
	    $seq_var .= $event->{sequence} || die;
	    $has_snv = 1;
	  } else {
	    die $type;
	  }
	}

	$seq_var .= $self->get_ref_single($ref_pos) unless $has_snv or $deleted{$ref_pos};
      }

      my $raw_key = join ".", $chr_raw, $start_pos, $seq_ref, $seq_var;
      printf STDERR "raw %d-%d => %s in %s\n", $start_pos, $end_pos, $raw_key, $aln->query->seq_id if $VERBOSE;

      my $final_start = $start_pos;
      #	die join ",", $final_start, $seq_ref, $seq_var;

      if (0) {
	# DISABLE:
	# adjusting the start position can interfere with counting
	# because insertions may appear in the alignment BEFORE
	# the adjusted site.
	while (substr($seq_ref, 0, 1) eq substr($seq_var, 0, 1)) {
	  # trim extraneous leading bases
	  $seq_ref = substr($seq_ref, 1);
	  $seq_var = substr($seq_var, 1);
	  $final_start++;
	}
      }

      while (substr($seq_ref, -1) eq substr($seq_var, -1)) {
	# can happen if trailing insertion
	$seq_ref = substr($seq_ref, 0, length($seq_ref) - 1);
	$seq_var = substr($seq_var, 0, length($seq_var) - 1);
      }

      my $usable;
      if ($gate_pos) {
	if ($gate_pos >= $final_start and $gate_pos <= $end_pos) {
	  $usable = 1;
	} else {
	  my $dist;
	  if ($gate_pos < $final_start) {
	    # site is near start gate position
	    $dist = $final_start - $gate_pos;
	  } elsif ($gate_pos >= $end_pos) {
	    $dist = $gate_pos - $end_pos;
	    # UGH: maybe need another threshold for this situation??
	    printf STDERR "read:%s start:%d end:%d gate:%d dist:%d\n",
	      $aln->query->seq_id(), $final_start, $end_pos, $gate_pos, $dist if $VERBOSE;
	  } else {
	    die "fix me";
	  }
	  printf STDERR "read:%s start:%d end:%d gate:%d dist:%d\n",
	    $aln->query->seq_id(), $final_start, $end_pos, $gate_pos, $dist if $VERBOSE;
	  $usable = $dist <= $MAX_NEIGHBOR_DISTANCE;
	  # maybe need another threshold for this situation?
	}
      } else {
	$usable = 1;
      }
      $usable = 0 if $qc_problem;

      my $key = join ".", $chr_raw, $final_start, $seq_ref, $seq_var;
      my $simple_indel;
      my $simple_usable;
      if (length($seq_ref) == 1 and length($seq_var) == 1) {
#	printf STDERR "skipping SNV $key\n";
	$usable = 0;
	# final usable data is only a SNV, ignore
      } elsif (length($seq_ref) == 0 or length($seq_var) == 0) {
	# final data is a simple indel, ignore
	foreach ($seq_ref, $seq_var) {
	  $_ = "-" unless $_;
	}
	$key = join ".", $chr_raw, $final_start, $seq_ref, $seq_var;
	$usable = 0;

	$simple_indel = 1;

	if ($gate_pos) {
	  my $dist = abs($final_start - $gate_pos);
	  $simple_usable = $dist <= $MAX_NEIGHBOR_DISTANCE;
	}
#	printf STDERR "skipping indel $key\n";
      }

      if ($usable) {
	printf STDERR "final %d-%d => %s in %s %s\n", $start_pos, $end_pos, $key, $aln->query->seq_id, $aln->query->strand if $VERBOSE;
	push @{$good{$key}}, $a;
      } elsif ($simple_indel) {
	push @{$simple{$key}}, $a if $simple_usable;
      } elsif ($VERBOSE and $qc_problem) {
	printf STDERR "QC problem in %s: %s\n", $aln->query->seq_id, $qc_problem;
      }
    }				# $ranges
  }				# $a

  my $min_good = $self->min_good_reads();
  foreach my $snv4 (keys %good) {
    my $count = @{$good{$snv4}};
#    printf STDERR "var:%s count:%d\n", $snv4, $count;
    delete $good{$snv4} if $count < $min_good;
  }
  foreach my $snv4 (keys %simple) {
    my $count = @{$simple{$snv4}};
#    printf STDERR "var:%s count:%d\n", $snv4, $count;
    delete $simple{$snv4} if $count < $min_good;
  }

  $self->called(\%good);
  $self->called_simple(\%simple);

  return scalar keys %good;
}

sub get_ref_single {
  my ($self, $pos) = @_;

  my $nt;
  if (1) {
    # fastest: load entire chrom into memory
    my $seq = $self->fai->get_sequence("-id" => $self->chr_raw());
    return uc(substr($$seq, $pos - 1, 1));
  } elsif (0) {
    # actually a bit slower than the below  :(
    my $cache = $self->nt_cache_single;
    $nt = $cache->{$pos};
    unless ($nt) {
      $nt = uc($self->fai->get_chunk_buffered(
		 "-id" => $self->chr_raw,
		 "-start" => $pos,
		 "-length" => 1
	       ));
      $cache->{$pos} = $nt;
    }
  } else {
    my $cache = $self->nt_cache_single;
    $nt = $cache->{$pos};
    unless ($nt) {
      # reliable but SLOW
      $nt = uc($self->fai->get_chunk(
		 "-id" => $self->chr_raw,
		 "-start" => $pos,
		 "-length" => 1
	       ));
      $cache->{$pos} = $nt;
    }
  }

  return $nt;
}

sub prune_ranges {
  # STATIC
  my ($ranges_all) = @_;
  my @passed;
  foreach my $r (@{$ranges_all}) {
    my $ok = 1;
    foreach my $other (@{$ranges_all}) {
      if ($r ne $other) {
	if ($r->[0] >= $other->[0] and $r->[1] <= $other->[1]) {
	  $ok = 0;
	  last;
	}
      }
    }
    push @passed, $r if $ok;
  }
  return \@passed;
}

sub found_simple {
  my ($self) = @_;
  return scalar keys %{$self->called_simple};
}

sub simple_has_best_evidence {
  my ($self) = @_;
  my $called = $self->called();
  my $called_simple = $self->called_simple();

  my @called = sort {$b <=> $a} map {scalar @{$called->{$_}}} keys %{$called};
  my @simple = sort {$b <=> $a} map {scalar @{$called_simple->{$_}}} keys %{$called_simple};

  my $best_called = @called ? $called[0] : 0;
  my $best_simple = @simple ? $simple[0] : 0;

  return $best_simple > $best_called;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
