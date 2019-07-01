package MapSNVToGenome;
# translate base position in a short DNA sequence to genomic

use strict;
use Configurable;
use Exporter;
use Carp qw(confess);

@MapSNVToGenome::ISA = qw(Configurable Exporter);
@MapSNVToGenome::EXPORT_OK = qw();

use BLATClient;
use GenomeUtils qw(cook_chromosome_name);
use FileUtils qw(find_binary);

#my $FLANK = 100;
my $FLANK = 25;
my $BLAST_EXPECT = .0001;

use MethodMaker qw(
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  return $self;
}

sub find {
  my ($self, %options) = @_;
  my $sequence = $options{"-sequence"} || die;
  # source sequence
  my $base_number = $options{"-base-number"} || die;
  my $base = $options{"-base"} ||  die;
  my $bi = $base_number - 1;

  my $min_start = $options{"-min-start"} || die "-min-start";
  my $max_end = $options{"-max-end"} || die "-max-end";
  # gate results to known start/end

  my $si = $bi - $FLANK;
  my $ei = $bi + $FLANK;
  die if $si < 0;
  die if $ei >= length($sequence);

  my $chunk = substr($sequence, $si, ($ei - $si));
  my $chunk_len = length($chunk);
  
  my $ci = $FLANK;
  # index of target base in excerpt sequence
#  die unless $base, substr($chunk, $ci, 1) eq $base;
  die unless substr($chunk, $ci, 1) eq $base;

  my $qbn = $ci + 1;
  # base number of target base in excerpt sequence

  my $target_chr = cook_chromosome_name($options{"-chrom"} || die "-chrom");

  my $id = sprintf '%s_%d_%s_idx_%d', $target_chr, $base_number, $base, $ci;
  my $blat_fa = sprintf 'blat_%s.fa', $id;
  unless (-s $blat_fa) {
    open(BP, ">" . $blat_fa) || die;
    printf BP ">%s\n%s\n", $id, $chunk;
    close BP;
  }

  my $p;
  if (1) {
    # BLAST
    # blastall -p blastn -d 13.fa -i blat_13_2281_C_idx_25.fa -o outfile
    # TEMPORARY HACK, FIX ME
    my $db_file = sprintf '%s.fa', $target_chr;
    die "where is $db_file" unless -s $db_file;
    my $outfile = $blat_fa . ".blast";
    unless (-s $outfile) {
      find_binary("blastall", "-die" => 1) || die "blastall not on path";
      my $cmd = sprintf 'blastall -p blastn -e %f -d %s -i %s -o %s -F 0',
      $BLAST_EXPECT, $db_file, $blat_fa, $outfile;
      system($cmd);
      die "$cmd exited with $?" if $?;
      die unless -s $outfile;
    }
    $p = Bio::SearchIO->new(-file => $outfile, -format => 'blast');
  } else {
    my $bc = new BLATClient();
    $bc->format("pslx");
    $p = $bc->query("-fasta" => $blat_fa);
    # unfortunately it appears bioperl will only parse .psl, 
    # so alignment data in .pslx not used?
  }

  my $mapped_base_number;

 OUTER:
  while (my $result = $p->next_result) {
    # one result object per query sequence
    my $qname = $result->query_name();
    printf STDERR "NEW query: %s\n", $qname;
    my @hsp_user_match;
    my @hsp_other;
    my $raw_hsp_count = 0;

    while (my $hit = $result->next_hit()) {
      # hits from this query to a database sequence
      # (Bio::Search::Hit::HitI)
      my $hit_refname = cook_chromosome_name($hit->name(), "-return-unknown" => 1);
      # accept non-standard hits, e.g. GL000228.1
      next unless $hit_refname eq $target_chr;

      printf STDERR "  hit to %s\n", $hit->name();
      # target name (chrom)

      my @all_hsp;
      # this is SUPER-RICKETY:
      # - only works with +
      while (my $hsp = $hit->next_hsp) {
	# high-scoring pairs within the hit
	# (Bio::Search::HSP::HSPI)
	printf STDERR "hsp %s\n", join ",", $hsp->range("query"), $hsp->num_identical, $hsp->frac_identical("query");

	my ($ss, $se) = $hsp->range("hit");
	if ($ss >= $min_start and $se <= $max_end) {
	  push @all_hsp, $hsp;
	}
      }

      my @sorted = sort {$b->score <=> $a->score} @all_hsp;
      my $best_hsp = $sorted[0];
      if (@sorted > 1) {
	if ($sorted[0]->score == $sorted[1]->score) {
	  # just in case
	  confess "FAIL: ambiguous hit!!";
	}
      }

      my $hsp = $best_hsp;

      if ($hsp) {
	$raw_hsp_count++;
	printf STDERR "    score:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d\n",
	$hsp->score,
	$hsp->strand,
	$hsp->range("query"),
	$hsp->range("hit"),
	$hsp->num_identical(),
	$hsp->frac_identical("query"),
	$hsp->length("query"),
	$hsp->length("hit"),
	$hsp->length("total");

	my ($qs, $qe) = $hsp->range("query");
	my ($ss, $se) = $hsp->range("hit");

	if ($qbn >= $qs and $qbn <= $qe) {
	  # query base is within this HSP
	  if ($hsp->strand <= 0) {
	    printf STDERR "ERROR, can't parse - strand\n";
	  } elsif ($hsp->gaps("query")) {
	    printf STDERR "ERROR: gap in query sequence\n";
	    printf STDERR "query=%s\n", $hsp->query_string();
	    printf STDERR "hit=%s\n", $hsp->hit_string();
	    printf STDERR "homology=%s\n", $hsp->homology_string();
	    die;
	  } elsif ($hsp->gaps("hit")) {
	    printf STDERR "ERROR: gap in hit sequence\n";
	    die;
	  } else {
	    my ($query_base_num, $subject_base_num);
	    for ($query_base_num = $qs, $subject_base_num = $ss;
		 $query_base_num < $qe;
		 $query_base_num++, $subject_base_num++) {
	      if ($query_base_num == $qbn) {
		print "mapped $query_base_num => $subject_base_num\n";
		$mapped_base_number = $subject_base_num;
		last OUTER;
	      }
	    }
	  }
	}
      }
    }
  }

  return $mapped_base_number;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
