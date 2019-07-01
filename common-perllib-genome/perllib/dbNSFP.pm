package dbNSFP;
# binary search dbNSFP flatfiles
# (e.g. /nfs_exports/genomes/1/Homo_sapiens/dbNSFP/)
# Michael Edmonson 12/2013
#
# TO DO:
# - NHLBI v3: filter to target transcript by translating ENST to refGene
#   (PITA as ENST IDs have historically aged/obsoleted rapidly)

use strict;
use Carp qw(confess cluck);

use FileHandle;
use POSIX qw(SEEK_SET);
use Time::HiRes qw(time);

use Configurable;
use GenomeUtils qw(cook_chromosome_name);
use AAParser;
use FileUtils qw(universal_open);
use DelimitedBinarySearch;

my $NSFP_DELIM = ";";
my $NSFP_NULL = ".";

my $VERBOSE = 0;

my @PRED_CLEAN_FIELDS = qw(
		     SIFT_pred
		     Polyphen2_HDIV_pred
		     Polyphen2_HVAR_pred
		     MutationAssessor_pred
		  );

@dbNSFP::ISA = qw(Configurable);

use MethodMaker qw(
	directory
chr2file
fh

headers
blankify

hits
raw_hits

uniprot_idmapping
unparsable_aa_ok
require_codon_number_match

dbs
f_pos

clean_pred_calls
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->blankify(1);
  # replace empty "." values with blanks
  $self->require_codon_number_match(0);
  # by default don't require codon number to match for AA filtering.
  # This USUALLY works but not always
  $self->f_pos("pos(1-coor)");
  # version 2.1
  $self->clean_pred_calls(1);
  # call field cleanups: remove null entries and duplicates
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self, %options) = @_;
  my $dir = $self->directory;
  if ($dir) {
    # not specified in tabix mode
    my @files = glob($dir . "/dbNSFP*_variant.chr*");
    die unless @files;
    my %chr2file;
    foreach my $fn (@files) {
      $fn =~ /chr(\w+)$/ || die;
      die "duplicate" if $chr2file{$1};
      my $chr = cook_chromosome_name($1);
      $chr2file{$chr} = $fn;
    }
    $self->chr2file(\%chr2file);
  }
}

sub get_fh_nearby {
  # approximate search: return filehandle upstream of target site
  # for manual searching e.g. of a region mapped to a gene
  my ($self, %options) = @_;
  my $chr = cook_chromosome_name($options{"-chr"} || die "-chr");
  my $wanted_pos = $options{"-pos"} || die "-pos";

  my $fn = $self->chr2file()->{$chr} || die "no file for $chr";

  my $dbs = new DelimitedBinarySearch(
				      "-file" => $fn,
				      "-is_headered" => 1,
				      "-max_line_length" => 2500
				      # chr.5140890570 = 2027
				     );
  $self->dbs($dbs);

  my $dbs_raw = $dbs->find(
			"-comparators" => [
					   {
					    "column_name" => "pos(1-coor)",
					    "type" => "number",
					    "value" => $wanted_pos
					   }
					  ],
			"-verbose" => $options{"-verbose"},
		       );

  return $dbs->get_fh_nearby();
}

sub find {
  # find hits associated with a genomic position.
  # Multiple results are possible.
  # if -aa is specified tries to winnow results to that event
  # likewise for -nm
  my ($self, %options) = @_;
  my $chr = cook_chromosome_name($options{"-chr"} || die "-chr");
  my $wanted_pos = $options{"-pos"} || die "-pos";
  my $ref_base = $options{"-reference-base"} || die "-reference-base";
  my $var_base = $options{"-variant-base"} || die "-variant-base";
  my $aa = $options{"-aa"};
  my $nm = $options{"-nm"};
  if ($nm and not($nm =~ /^NM_/)) {
    # sometimes bad data in SJ reports
    cluck sprintf "dbNSFP: can't use -nm \"%s\", not in NM_ format!", $nm if $VERBOSE;
    $nm = undef;
  }
  my $pre_hits = $options{"-hits"};

  my $start_time = time();

  my $is_indel;
  foreach ($ref_base, $var_base) {
    $is_indel = 1 if length($_) > 1 or /\-/;
  }
  if ($is_indel) {
    printf STDERR "dbNSFP: can't look up indels/MNVs: ref=%s var=%s\n", $ref_base, $var_base if $VERBOSE;
    my @hits;
    $self->hits(\@hits);
    $self->raw_hits(\@hits);
    return \@hits;
  }

#  printf STDERR "new search: %s\n", join ".", $chr, $wanted_pos, $ref_base, $var_base;

  my $uniprot2nm;
  if ($nm) {
    confess "broken NM_ format: $nm" unless $nm =~ /^NM_/;
    if (0) {
      print STDERR "\n ****  DEBUG, skip uniprot load **** \n\n";
      $uniprot2nm = {};
      $uniprot2nm->{"P54278"}{"NM_000535"} = 1;
    } else {
      $uniprot2nm = $self->get_uniprot2nm();
    }
  }

  my $dbs_raw;
  if ($pre_hits) {
    $dbs_raw = $pre_hits;
  } else {
    my $fn = $self->chr2file()->{$chr} || confess "no file for $chr";
    die unless -s $fn;

    my $dbs = new DelimitedBinarySearch(
					"-file" => $fn,
					"-is_headered" => 1,
					"-max_line_length" => 2500
					# chr.5140890570 = 2027
				       );

    $dbs_raw = $dbs->find(
			  "-comparators" => [
						{
						 "column_name" => "pos(1-coor)",
						 "type" => "number",
						 "value" => $wanted_pos
						}
					       ],
			     "-verbose" => $options{"-verbose"},
			    );
    $self->headers($dbs->headers);
  }

  my @hits;
  my $blankify = $self->blankify();
  my $f_pos = $self->f_pos() || die;

  # detect v2/v3 format and field names here if @{$dbs_raw}

  my $dbnsfp_version;
  my $f_uniprot_acc;

  foreach my $hit (@{$dbs_raw}) {
    die unless $hit->{$f_pos} == $wanted_pos;

    unless ($dbnsfp_version) {
      if (exists $hit->{"Uniprot_acc"}) {
	$f_uniprot_acc = "Uniprot_acc";
	$dbnsfp_version = 2;
      } elsif (exists $hit->{"Uniprot_acc_Polyphen2"}) {
	$f_uniprot_acc = "Uniprot_acc_Polyphen2";
	$dbnsfp_version = 3;
      } else {
	dump_die([$hit], "unknown dbNSFP data version");
      }
    }

    if ($blankify) {
      foreach (values %{$hit}) {
	$_ = "" if $_ eq ".";
      }
    }
    push @hits, $hit if uc($ref_base) eq $hit->{ref} and
      uc($var_base) eq $hit->{alt};
  }

  my $end_time = time();
#  printf STDERR "dbNSFP search time for %s: %f\n", join(".", $chr, $wanted_pos, $ref_base, $var_base), $end_time - $start_time;

  my $aa_parsable;
  my $aa_silent;
  my @raw_hits = @hits;
  $self->raw_hits([ @raw_hits ]);

  my %saw_nm;
  if ($nm) {
    #
    #  filter results to specified NM_ record via UniProt lookup.
    #
    my $user_nm = $nm;
    $user_nm =~ s/\.\d+$//;
    # don't require sub-version match

    my @filtered;
    foreach my $r (@hits) {
      my $acc_list = $r->{$f_uniprot_acc};
      my %wanted_indices;
      if ($acc_list) {
	my @accs = split /$NSFP_DELIM/, $acc_list;

	#
	#  first pass: try looking up UniProt accessions
	#  *EXACTLY* as they appear in dbNSFP.  These
	#  may or may not have a sub-isoform number,
	#  or even a mix of unversioned and versioned
	#  for the SAME UniProt ID (e.g. "O95197;O95197-2").
	#
	for (my $i = 0; $i < @accs; $i++) {
	  my $acc = $accs[$i];
	  if (my $map = $uniprot2nm->{$acc}) {
	    # some identifiers don't map to NM_, e.g. C9J167
	    foreach my $nm (keys %{$map}) {
	      printf STDERR "hit for %s: %s\n", $acc, $nm if $VERBOSE;
	      if ($nm eq $user_nm) {
		$wanted_indices{$i} = 1;
		printf STDERR "  hit user NM $nm\n" if $VERBOSE;
	      }
	      $saw_nm{$nm} = 1;
	    }
	  }
	}

	unless (%wanted_indices) {
	  #
	  #  try alternate lookups, but ONLY if the primary lookup fails.
	  #
	  for (my $i = 0; $i < @accs; $i++) {
	    my $acc_raw = $accs[$i];

	    my @try;
	    if ($acc_raw =~ /\-\d+$/) {
	      # dbNSFP ID contains a sub-version number, try without.
	      my $stripped = $acc_raw;
	      $stripped =~ s/\-\d+$//;
	      push @try, $stripped;
	    } else {
	      # dbNSFP does not mention the sub-version number.
	      # sometimes the index contains ONLY versioned entries,
	      # e.g. for Q2M3G0 index has Q2M3G0-1 through Q2M3G0-4.
	      #
	      # ALSO "-1" is not guaranteed to be present.
	      # e.g. for Q8N398, there is an entry for Q8N398-2 only:
	      # gzip -dc HUMAN_9606_idmapping.dat.gz |egrep 'E9PF42|B9EGN7|Q8N398'|grep -i refseq
	      # Q8N398-2        RefSeq  NP_612354.1
	      # Q8N398-2        RefSeq_NT       NM_138345.1
	      # Q8N398-2        RefSeq  XP_006713882.1
	      # Q8N398-2        RefSeq_NT       XM_006713819.1
	      # E9PF42  RefSeq  XP_006713881.1
	      # E9PF42  RefSeq_NT       XM_006713818.1
	      #
	      # The example earlier of "O95197;O95197-2"
	      # demonstrates why this must only be done AFTER the
	      # primary lookup fails.  In that case the 2nd raw
	      # entry is the one we want.  Trying this alternate
	      # lookup first would have failed: O95197 would have
	      # been extended to "O95197-1" which in this case
	      # does NOT match the accession in O95197-2.
	      for (my $num = 1; $num < 40; $num++) {
		# observed 37
		push @try, $acc_raw . "-" . $num;
	      }
	    }

	    foreach my $acc (@try) {
	      if (my $map = $uniprot2nm->{$acc}) {
		# some identifiers don't map to NM_, e.g. C9J167
		foreach my $nm (keys %{$map}) {
		  printf STDERR "secondary hit for %s: %s\n", $acc, $nm if $VERBOSE;
		  if ($nm eq $user_nm) {
		    $wanted_indices{$i} = 1;
		    printf STDERR "  secondary hit user NM $nm\n" if $VERBOSE;
		  }
		  $saw_nm{$nm} = 1;
		}
	      }
	      last if %wanted_indices;
	      # if a hit found to desired lookup, stop.
	    }
	  }
	}

	if (%wanted_indices) {
	  # this record contains data for the desired accession

	  if (scalar keys %wanted_indices > 1) {
	    #
	    # even with filtering there are multiple usable indexes
	    #
	    my @accs = split /$NSFP_DELIM/, $r->{$f_uniprot_acc};
	    my @canonical_i;
	    my @secondary_i;
	    foreach my $i (sort {$a <=> $b} keys %wanted_indices) {
	      if ($accs[$i] =~ /\-\d+$/) {
		# secondary isoform, e.g. P54278-2
		push @secondary_i, $i;
	      } else {
		# canonical isoform, e.g P54278
		push @canonical_i, $i;
	      }
	    }

	    if (@canonical_i and @secondary_i) {
	      # if record contains both canonical and secondary IDs,
	      # remove secondary ones
	      delete @wanted_indices{@secondary_i};
	    }
	  }

	  #
	  # additional resolution ideas:
	  #  - choose the entry where Uniprot_aapos best matches "aapos"
	  #  - choose the entry where the UniProtKB-ID is NOT named after
	  #    the uniprot ID:
	  # P25054  UniProtKB-ID    APC_HUMAN (reviewed)
	  # Q4LE70  UniProtKB-ID    Q4LE70_HUMAN (not reviewed)
	  #  - choose the entry having the most unique database tags
	  #    44 tags_P25054
	  #    23 tags_Q4LE70
	  #  - choose the entry having the most data rows:
	  #    P25054 = 139
	  #    Q4LE70 = 30
	  #

	  if (scalar keys %wanted_indices > 1) {
	    my @accs = split /$NSFP_DELIM/, $r->{$f_uniprot_acc};

#	    dump_die([ $r ], "debug", 1);

	    printf STDERR "attempting to resolve ambiguous UniProt IDs for %s: %s...",
	    join(".", $chr, $wanted_pos, $ref_base, $var_base, $nm),
	    join(", ", map {$accs[$_] . "=" . $dbNSFP::UNIPROT2KBID{$accs[$_]}} sort {$a <=> $b} keys %wanted_indices) if $VERBOSE;

	    my %score;
	    foreach my $i (keys %wanted_indices) {
	      my $acc = $accs[$i];
	      $acc =~ s/\-\d+$//;
              # don't look up sub-indexes
	      my $score = $dbNSFP::UNIPROT_COUNTS{$acc} || die "no counts for $acc!";
	      if (my $kbid = $dbNSFP::UNIPROT2KBID{$acc}) {
		$score = 0 if $kbid =~ /$acc/;
		# attempting to resolve ambiguous UniProt IDs for 12.121176083.G.A.NM_000017: E5KSD5=E5KSD5_HUMAN, P16219=ACADS_HUMAN
		# winner by counts alone is E5KSD5, FAIL (unreviewed)
		# reviewed entries seem to contain gene name rather 
		# than accession
	      }
	      $score{$i} = $score;
	    }

	    my ($best, @others) = sort {$score{$b} <=> $score{$a}} keys %wanted_indices;

	    foreach ($best, @others) {
	      printf STDERR "  index %d: score %d\n", $_, $score{$_} if $VERBOSE;
	    }

	    if ($score{$best} > $score{$others[0]}) {
	      # logic worked
	      printf STDERR "winner chosen by best score!\n" if $VERBOSE;
	    } else {
	      my @candidates = grep {$score{$_} == $score{$best}} keys %wanted_indices;
	      # might have already narrowed the list down by score method
	      #
	      # e.g. chr13.36413279.G.A   NM_004734  A357V

	      $best = tiebreak_pph2($r, \@candidates);
	      if (defined $best) {
		# picked a winner by highest PolyPhen2 damage prediction
		printf STDERR "yay: PPH2 tiebreak\n" if $VERBOSE;
	      } else {
		$best = (sort {$a <=> $b} @candidates)[0];
		# just choose the first ID for consistency
		printf STDERR "PPH2 tiebreak failed, bailing out.  candidates=%s winner=%d\n", join(",", @candidates), $best if $VERBOSE;
	      }
	      @others = grep {$_ ne $best} keys %wanted_indices;
	    }

	    delete @wanted_indices{@others};
	    my $winner = $accs[$best];
	    printf STDERR "winner: %s\n", $winner if $VERBOSE;
	    $winner =~ s/\-\d+$//;
	    my $kbid = $dbNSFP::UNIPROT2KBID{$winner} || die "no kbid for $winner";

#	    die "uh-oh: resolved to possibly dubious ID" if $kbid =~ /$winner/;
	    # sometimes this happens despite our best efforts, e.g.:
	    #
	    # "-chr" => 18,
	    # "-pos" => 13068878,
	    # "-reference-base" => "C",
	    # "-variant-base" => "T",
	    # "-nm" => "NM_032142",
	    #
	    # ...hits E9PF99 and Q9HCK3 (both unreviewed).
	  }

#	  dump_die([$r], "debug, before splice", 1);

	  my @splice_fields = (
			       $f_uniprot_acc,
			       "Polyphen2_HDIV_score",
			       "Polyphen2_HDIV_pred",
			       "Polyphen2_HVAR_score",
			       "Polyphen2_HVAR_pred",
			      );
	  if ($dbnsfp_version == 2) {
	    push @splice_fields, (
	      "Uniprot_id",
	      "Uniprot_aapos",
	    );
	  } elsif ($dbnsfp_version == 3) {
	    push @splice_fields, (
				  "Uniprot_id_Polyphen2",
				  "Uniprot_aapos_Polyphen2",
				 );
	  } else {
	    die "unknown split fields";
	  }

	  foreach my $f (@splice_fields) {
	    splice_field_data($r, $f, \%wanted_indices);
	  }

	  push @filtered, $r;
	}
      }
    }

    if (@filtered) {
      @hits = @filtered;
    } else {
      # nothing made it through
      if (%saw_nm) {
	# NM_ accessions could be found, but none mapped.
	# Note these are NOT the raw UniProt accessions:
	# sometimes these are present but don't map to NM, e.g.
	# 11.1265752.A.G has A7Y9J9;E9PBJ0 but neither has an NM_
	my $msg = sprintf "dbNSFP lookup failed: aa=%s wanted nm=%s, found=%s; keeping anyway\n", $aa, $nm, join ",", sort keys %saw_nm;
	dump_die(\@hits, $msg, 1) if $VERBOSE;
	# might need further investigation
	# examples of why we should raw keep results in any case:
	#
	# 1. aa=Q6765L wanted NM_001164507, only found NM_004543
	#    appears to be a different but compatible accession
	#   (ref/alt codon and # are the same)
	#
	# 2. aa=D116_E2splice wanted NM_004048, found XM_005254549
	#    - fuzzy codon match successful (D116)
	#    - hits F5H6I0 instead (predicted)
      } else {
	# since no NM_ accessions could be found, filtering is not possible
      }
    }
  }

#  printf STDERR "hits after phase 1: %d\n", scalar @hits;

#  dump_die(\@hits, "phase 1", 1);

  if ($aa) {
    #
    #  filter to records matching user-specified AA annotation
    #
    my $p = new AAParser();
    if (my $cooked = $p->parse_substitution($aa)) {
      # well-formed annotation
      $aa_parsable = 1;
      my $codon_reference = $p->codon_reference;
      my $codon_number = $p->codon_number;
      my $codon_variant = $p->codon_variant;
#      dump_die(\@hits, "reference is stop codon") if $codon_reference eq "*";

      foreach ($codon_reference, $codon_variant) {
	$_ = "X" if $_ eq "*";
	# dbNSFP encodes these as an X  :/
      }

      $aa_silent = 1 if $codon_reference eq $codon_variant;
      # 18.47431159.C.T = V818V

      my @filtered;
      my @filtered_unknown;
      my $comparable;
      foreach my $r (@hits) {
	if ($r->{aaref} and
	    $r->{aaalt} and
	    $r->{aapos} and
	    ($r->{aapos} || "") ne -1) {
	  # entry has info we can compare against
	  $comparable = 1;
	  my $ok = 1;
	  $ok = 0 unless $r->{aaref} eq $codon_reference;
	  $ok = 0 unless $r->{aaalt} eq $codon_variant;

	  if ($self->require_codon_number_match()) {
	    my $codon_number_ok;
	    foreach my $pos_field (qw(aapos Uniprot_aapos)) {
	      if (my $v = $r->{$pos_field}) {
		foreach my $cn (split /$NSFP_DELIM/, $v) {
		  $codon_number_ok = 1 if $cn == $codon_number;
		}
	      }
	    }
	    $ok = 0 unless $codon_number_ok;
	  }
	  push @filtered, $r if $ok;
	} else {
	  push @filtered_unknown, $r;
	  # there is a record, but it's unannotated
	}
      }
      @hits = @filtered if $comparable;
      # - if annotations are present, require that they match,
      #   (even if that means discarding all results).
      # - if annotations are missing, allow results through untouched
      #   QUESTION: should we do this??
    } elsif ($p->parse($aa)) {
      # non-standard AA annotation, e.g.
      # 15.45007900.G.A D116_E2splice
      if ($p->codon_start == $p->codon_end) {
	my $codon_number = $p->codon_start || die;
	$aa_parsable = 2;
	
	my @filtered;
	foreach my $r (@hits) {
	  my $ok;
	  foreach my $cn (split /$NSFP_DELIM/, $r->{aapos}) {
	    $ok = 1 if $cn == $codon_number;
	  }
	  if ($ok) {
	    push @filtered, $r;
	    printf STDERR "WARNING: dbNSFP fuzzy match, CHECK ME: aa=%s aaref=%s aaalt=%s\n", $aa, @{$r}{qw(aaref aaalt)} if $VERBOSE;
	  }
	}

	if (not(@filtered) and
	    @hits == 1 and 
	    $hits[0]->{aapos} == -1 and
	    $aa =~ /splice/i
	    ) {
	  @filtered = @hits;
	  # many splice-annotated variants seem to have a single entry
	  # but without AA annotations.  Maybe these are not always
	  # available/complete?  Pass through for now.
	  printf STDERR "dbNSFP: passing through single non-AA genomic lookup for %s.%s.%s.%s AA %s\n", $chr, $wanted_pos, $ref_base, $var_base, $aa if $VERBOSE;
	}
	@hits = @filtered;
      } else {
	printf STDERR "dbNSFP: AA event %s spans codon, can't look up!\n", $aa if $VERBOSE;
	@hits = ();
      }
    } else {
      printf STDERR "WARNING: can't parse AA %s\n", $aa if $VERBOSE;
      # might not be fatal, only a problem if results are ambiguous
      # and AA is needed for disambiguation
    }

#    printf STDERR "hits after phase 2: %d\n", scalar @hits;

    #
    #  debug/error checks:
    #
    my $error_msg = sprintf "dbNSFP: lookup for %s.%s.%s.%s AA=%s NM=%s, result_count=%d, raw_count=%d raw_aaref=%s raw_aapos=%s raw_aaalt=%s saw_nm=%s", $chr, $wanted_pos, $ref_base, $var_base, $aa, $nm, scalar(@hits), scalar(@raw_hits),
    join(",", map {$_->{aaref}} @raw_hits),
    join(",", map {$_->{aapos}} @raw_hits),
    join(",", map {$_->{aaalt}} @raw_hits),
    join(",", sort keys %saw_nm);
    if (@hits > 1) {
      # ambiguity is always trouble
      dump_die(\@hits, $error_msg, 1);
      # however still some sticky cases where not sure what to do,
      # e.g.
      # 20.30956815.C.A AA E3_exon
    } elsif (@hits == 0) {
      # no results
      if ($aa_parsable) {
	# expect results in this case
	if ($aa_silent) {
	  # dbNSFP doesn't contain silent
	} else {
	  # confess($error_msg);
	  if (@raw_hits) {
	    # check to make sure not some kind of synchronization error,
	    # e.g. for stop codons
	    dump_die(\@raw_hits, $error_msg . ", all hits removed", 1) if $VERBOSE;
	  } else {
	    printf STDERR "%s\n", $error_msg if $VERBOSE;
	  }
	  
	  # sometimes results are not present, not sure why, e.g.
	  # 1.16915439.G.A AA P133L, result_count=0, raw_count=0 
	  # NBPF1 / NM_017940
	  # so now warn rather than die
	}
      } elsif ($aa and not($self->unparsable_aa_ok)) {
	# not sure here: AA might be relevant but parsing is not
	# yet handled properly
	if ($aa =~ /utr/i or
	    $aa =~ /intron/i
	    ) {
	  # safe to ignore
	} elsif ($aa =~ /splice_region/i) {
	  # splice region: no data in db?
	} elsif ($aa =~ /_exon/i or
		 # E1_exon
		 $aa =~ /E\d+_splice/i or
		 # E5_splice
		 $aa =~ /splice_E\d/i or
		 # splice_E3
		 $aa =~ /_E\d+SPLICE$/i
		 # V80_E4SPLICE
	    ) {
	  # might be uppercased after munging by ?SITH/GEDI?
	} elsif (lc($aa) eq "unknown") {
	  # temporary hack classification for NHLBI common variants
	} else {
	  confess "FATAL ERROR: $error_msg";
	}
      } else {
	# no AA annotation, ignore
      }
    }
  }

  if ($self->clean_pred_calls and @hits) {
    # remove duplicates from results, e.g. replace P;.;P with P.
    # also condenses SIFT predictions for dbNSFP v3, which are indexed
    # by ENSEMBL transcript ID, disambiguation not currently supported
    foreach my $hit (@hits) {
      foreach my $f (@PRED_CLEAN_FIELDS) {
	my $fv = $hit->{$f};
	if ($fv and $fv =~ /$NSFP_DELIM/) {
	  my @v = split /$NSFP_DELIM/, $fv;
	  my @unique;
	  my %saw;
	  foreach my $v (@v) {
	    next if $v eq $NSFP_NULL or $saw{$v};
	    push @unique, $v;
	    $saw{$v} = 1;
	  }
	  $hit->{$f} = join $NSFP_DELIM, @unique;
	}
      }
    }
  }


#  printf STDERR "dbNSFP: found results\n" if @hits;
  $self->hits(\@hits);

  return \@hits;
}

sub dump_die {
  my ($rows, $msg, $warn) = @_;
  printf STDERR "dump_die(): dumping %d db hits:\n", scalar @{$rows};
  foreach my $row (@{$rows}) {
    printf STDERR "NEW ROW:\n";
    foreach (sort keys %{$row}) {
      printf STDERR "  %s: %s\n", $_, $row->{$_};
    }
  }
  if ($warn) {
    printf STDERR "WARNING: %s\n", $msg;
  } else {
    confess($msg);
  }
}

sub get_uniprot2nm {
  my ($self) = @_;
  unless (%dbNSFP::UNIPROT2NM) {
    # singleton
    my $map_file = $self->uniprot_idmapping() || die "need -uniprot_idmapping";
    # ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
    printf STDERR "parsing %s...", $map_file;
    my $fh = universal_open($map_file) || die;
    my $saw_subisoforms;
    while (<$fh>) {
      chomp;
      my ($uniprot, $type, $id) = split /\t/, $_;
      $saw_subisoforms = 1 if $uniprot =~ /\-\d+$/;

      $dbNSFP::UNIPROT_COUNTS{$uniprot}++;
      # hack to help resolve hits to multiple uniprot entries, e.g.
      # "-chr" => 5,
      # "-pos" => 112173767,
      # "-reference-base" => "T",
      # "-variant-base" => "G",
      # "-nm" => "NM_000038",

      $id =~ s/\.\d+$//;
      # don't require sub-version match
      if ($type eq "RefSeq_NT") {
	# may map to multiple entries, e.g.
	# P31946  RefSeq_NT       NM_003404.4
	# P31946  RefSeq_NT       NM_139323.3
#	printf STDERR "map uniprot %s => %s\n", $uniprot, $id;
	$dbNSFP::UNIPROT2NM{$uniprot}{$id} = 1;
      } elsif ($type eq "UniProtKB-ID") {
	die "duplicate $uniprot" if $dbNSFP::UNIPROT2KBID{$uniprot};
	$dbNSFP::UNIPROT2KBID{$uniprot} = $id;
      }
    }
    close $fh || die "i/o error $! $?";
    die "id parsing error" unless scalar keys %dbNSFP::UNIPROT2NM;

    die "ERROR: didn't find any sub-isoform identifiers (required for refSeq mapping).  Newer UniProt ID mapping file needed?" unless $saw_subisoforms;
    print STDERR "\n";
  }
  return \%dbNSFP::UNIPROT2NM;
}

sub splice_field_data {
  # static
  my ($row, $field, $wanted) = @_;
  die unless exists $row->{$field};
  $row->{$field . "_raw"} = $row->{$field};
  my @raw = split /$NSFP_DELIM/, $row->{$field};
  if (@raw > 1) {
    # in v3, some columns may contain only one value.
    # may be a compensation for pointlessly repetitive entries, e.g.
    # B;B;B is now B?
    my @filtered;
    for (my $i=0; $i < @raw; $i++) {
      push @filtered, $raw[$i] if $wanted->{$i};
    }
    $row->{$field} = join $NSFP_DELIM, @filtered;
  }
}

sub get_uniprot_aachange {
  # return AA change from (possibly filtered) Uniprot list
  #
  # reannotation of meta variants list, unique by position:
  # - traditional (requires single UniProt AA pos):
  #   - 8400 AA annotations (raw)
  #   - 6424 (removing Term, i.e. passthrough)
  #   - 519 discrepancies
  # - allowing "aapos" if UniProt pos not available:
  #   - 8499 AA annotations (raw)
  #   - 7036 (removing Term)
  #   - 542 discrepancies
  # - if both of the above fail but we have an ambiguous list of
  #   codon #s, use the first one:
  #   - 8761 AA annotations (Raw)
  #   - 8645 (removing Term)
  #   - 1012 discrepancies
  my ($self, %options) = @_;
  my $row = $options{"-row"} || die "-row";
  my $strict = $options{"-strict"} || 0;
  my $relax_if_nonsense = $options{"-relax-if-nonsense"};

  my $aa;

  my $aa_ref = $row->{aaref};
  my $aa_alt = $row->{aaalt};

  my @multi;

  my @fields_to_search = "Uniprot_aapos";
  # primary, often pre-filtered to just our transcript of interest
  push @fields_to_search, "aapos" unless $strict;
  # sometimes however not available, use ordinary codon # field.
  # more dangerous: not necessarily the transcript we want!

  if ($aa_ref and $aa_alt) {
    $aa_alt = "*" if $aa_alt eq "X";
    # hack: always true?
    foreach my $field (@fields_to_search) {
      if (my $pos = $row->{$field}) {
	my @p = split /$NSFP_DELIM/, $pos;
	my %v = map {$_, 1} @p;
	if (scalar keys %v == 1) {
	  $aa = join "", $aa_ref, $p[0], $aa_alt;
#	  dump_die([$row], "hey now $aa") if $field eq "aapos";
	  last;
	} elsif (scalar keys %v > 1) {
	  push @multi, \@p;
	}
      }
    }
  }

  if (not($aa) and not($strict) and @multi) {
    # ugly: no single-codon annotation is available, but we do
    # have a list of codon numbers.
    # even more dangerous: known ambiguity!!
    # (last ditch)
    $aa = join "", $aa_ref, $multi[0]->[0], $aa_alt;
#    dump_die([ $row ], "hey now $aa");
  }

  if (not($aa) and $relax_if_nonsense) {
    my $aap = new AAParser();
    if ($aap->parse_substitution($relax_if_nonsense)) {
      # nonsense mutations often seem to be missing uniprot info
      if ($aap->is_nonsense()) {
	dump_die([ $row ]);
      }
    }

  }

  return $aa;
}

sub tiebreak_pph2 {
  # attempt tie-breaking based on HDIV
  my ($row, $indices) = @_;
  my $field = "Polyphen2_HDIV_pred";
  die unless exists $row->{$field};
  my @values = split /$NSFP_DELIM/, $row->{$field};
  my %scores;
  my %value2score = (
		     "." => 0,
		     "B" => 1,
		     "P" => 2,
		     "D" => 3
      );
  my $verbose = 0;

  my %wanted = map {$_ => 1} @{$indices};
  printf STDERR "start: %s\n", join ",", @{$indices} if $verbose;

  for (my $i = 0; $i < @values; $i++) {
    next unless $wanted{$i};
    my $score = 0;
    my $v = $values[$i];
    if ($v) {
      $score = $value2score{$v};
      printf STDERR "index %d => %s => %d\n", $i, $v, $score if $verbose;
#      die sprintf "no score for value %s, row=%s", $v, join ",", @{$row} unless defined $score;
      dump_die([$row], "no score for value $v") unless defined $score;
      $scores{$i} = $score;
    }
  }

  my ($best_score) = sort {$b <=> $a} values %scores;
  my $winner;
  if ($best_score) {
    my @nominees = grep {$scores{$_} == $best_score} @{$indices};
    ($winner) = sort {$a <=> $b} @nominees;
    # if multiple entries with the same score, choose the first index
  }
  return $winner;
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
