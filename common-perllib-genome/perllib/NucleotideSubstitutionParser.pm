package NucleotideSubstitutionParser;
# parse nucleotide substitutions in c.* format
# found in COSMIC, etc.
# see http://www.hgvs.org/
# This code also handles a lot of dirty real-world data.
# MNE 2013-

use strict;
use Configurable;
use Exporter;
use Carp qw(cluck confess);

use GenomeUtils qw(reverse_complement);

use constant REF_SEQ_UNSPECIFIED => "unspecified";
use constant INDEL_CHAR => "-";

@NucleotideSubstitutionParser::ISA = qw(Configurable Exporter);

use MethodMaker qw(
	auto_strand_fix
is_parsable
is_substitution
is_insertion
is_deletion
is_complex_indel
is_coding_intronic
reference_sequence
variant_sequence
event_length
complex_insertion_length
complex_deletion_length
start
end
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  return $self;
}

sub parse {
  my ($self, $nt, %options) = @_;
  my $parsable = 0;

  $nt =~ s/^[cgm]\.//;
  # m.: NC_012920.1:m.3308T>C (Clinvar, mitochondrion)

  my $is_substitution = 0;
  my $is_insertion = 0;
  my $is_deletion = 0;
  my $is_complex_indel = 0;
  my $is_coding_intronic = 0;
  my $start = 0;
  my $end = 0;
  my $reference_sequence = "";
  my $variant_sequence = "";
  my $event_length = 0;
  my $complex_insertion_length = 0;
  my $complex_deletion_length = 0;
  
  if ($nt =~ /^(\d+)([A-Z])\>([A-Z])$/) {
    # SNV, e.g. COSMIC c.733C>T
    $parsable = 1;
    $is_substitution = 1;
    ($start, $end, $reference_sequence, $variant_sequence) = ($1, $1, $2, $3);
    $event_length = 1;
  } elsif ($nt =~ /^(\d+)_(\d+)del(\d+)$/) {
    # range deletion, e.g. COSMIC c.4409_4422del14
    $parsable = 1;
    $is_deletion = 1;
    $start = $1;
    $end = $2;
    $event_length = $3;
    $reference_sequence = "N" x $event_length;
    cluck  "length mismatch $nt" unless ($start + $event_length) == ($end + 1);
    # IARC TP53: 7579401_7579401del3
    # (*sigh*)
    # COSMIC? (via erin): 7_378del1016
    $variant_sequence = INDEL_CHAR;
  } elsif ($nt =~ /^(\d+)_(\d+)del$/) {
    # range deletion, e.g. NHGRI BRCA1 g.41232401_41236235del
    $parsable = 1;
    $is_deletion = 1;
    $start = $1;
    $end = $2;
    $event_length = ($end - $start) + 1;
    $variant_sequence = INDEL_CHAR;
  } elsif ($nt =~ /(\d+)del([A-Z])$/) {
    # single-base deletion, e.g. COSMIC c.750delG
    $parsable = 1;
    $is_deletion = 1;
    $start = $end = $1;
    $reference_sequence = $2;
    $variant_sequence = INDEL_CHAR;
    $event_length = 1;
  } elsif ($nt =~ /^(\d+)_(\d+)del([A-Z]+)$/) {
    # range deletion, e.g. COSMIC c.4479_4482delGGAA
    $parsable = 1;
    $is_deletion = 1;
    $start = $1;
    $end = $2;
    $reference_sequence = $3;
    $variant_sequence = INDEL_CHAR;
    $event_length = length($3);
    cluck sprintf "ERROR: length mismatch in $nt\n" unless ($start + $event_length) == ($end + 1);
    # RB1: 156712_156715delCCG
  } elsif ($nt =~ /^(\d+)_(\d+)ins([A-Z]+)$/) {
    # insertion of one or more bases, e.g. COSMIC c.3860_3861insGATGAAAT
    $parsable = 1;
    $is_insertion = 1;
    ($start, $end, $variant_sequence) = ($1, $2, $3);
    $event_length = length $variant_sequence;
    $reference_sequence = INDEL_CHAR;
    cluck sprintf "ERROR: range mismatch in %s", $nt unless $end - $start == 1;
    # RB1: 170380_170380insT
#    die unless $end - $start == 1;
  } elsif ($nt =~ /(\d+)_(\d+)ins(\d+)$/) {
    # insertion of a length of bases (not specified)
    # e.g. COSMIC c.1020_1021ins10
    $parsable = 1;
    $is_insertion = 1;
    ($start, $end, $event_length) = ($1, $2, $3);
    $reference_sequence = INDEL_CHAR;
    cluck "range mismatch $nt" unless $end - $start == 1;
    $variant_sequence = "N" x $event_length;
    # TP53: 7578554_7579312ins3
    # actually a complex event? i.e. large deletion and small insertion???
    # (*sigh*)
  } elsif ($nt =~ /^(\d+)([\+\-]\d+)([A-Z]+)\>([A-Z]+)$/) {
    # intronic variant in coding DNA, e.g.:
    # c.88+2T>G
    # c.89-1G>T
    # http://www.hgvs.org/mutnomen/examplesDNA.html
    # - does NOT fit well in this package's model!!!
#    die "hey now $nt $1 $2 $3 $4";
    $parsable = 1;
    # SO-SO: intronic offset is problematic.
    # however, if we are retrieving genomic coordinates elsewhere
    # base parsing components might still be useful.
    $start = $end = $1;
    my $offset = $2;
    $reference_sequence = $3;
    $variant_sequence = $4;
    $is_coding_intronic = 1;
    if (length($reference_sequence) == length($variant_sequence)) {
      # substitution
      $is_substitution = 1;
      $event_length = length $reference_sequence;
    } else {
      if (length($reference_sequence) < length($variant_sequence) and
	  index($variant_sequence, $reference_sequence) == 0) {
	# e.g. c.3123+1G>GA
	# net, appears to be an insertion of A rather than a complex event
	my $inserted = substr($variant_sequence, length($reference_sequence));
	$is_insertion = 1;
	$event_length = length($inserted);
	$reference_sequence = "";
	$variant_sequence = $inserted;
      } else {
	# complex?, e.g. COSMIC 151-16CTGTT>TGAGGCCCC
	$is_complex_indel = 1;
	$event_length = -1;
	$complex_deletion_length = length($reference_sequence);
	$complex_insertion_length = length($variant_sequence);
      }
    }
  } elsif ($nt =~ /^(\d+)([\+\-]\d+)del([acgt]+|\d+)$/i) {
    # intronic COSMIC deletion (single)
    # intronic variant in coding DNA, e.g.
    # c.1608-3delt
    # c.669+1delg
    # c.250+2delGT
    # 
    # or a size:
    # c.1643-20del12

    $parsable = 1;
    $is_deletion = 1;
    $is_coding_intronic = 1;
    $start = $end = $1;
    my $offset = $2;
    $reference_sequence = uc($3);
    if ($reference_sequence =~ /^\d+$/) {
      # c.1643-20del12
      $reference_sequence = "N" x $reference_sequence;
    }

    $event_length = length($reference_sequence);
    $variant_sequence = INDEL_CHAR;
#  } elsif ($nt =~ /^(\d+)([\+\-]\d+)_(\d+)([\+\-]\d+)del([acgt]+)$/) {
#  } elsif ($nt =~ /^(\d+)([\+\-]\d+)_(\d+)([\+\-]\d+)del([acgt]+|\d+)$/) {
  } elsif ($nt =~ /^(\d+)([\+\-]\d+)?_(\d+)([\+\-]\d+)?del([acgt]+|\d+)$/i and
	   ($2 or $4)) {
    # intronic COSMIC deletion, range version
    # - at least one side must have a +/- offset
    # - sequence may be either nucleotides or length

    # c.655-12_655-9delaacc
    # c.1114+1_1115-1delaggacacatcaactg

    # c.1325_1335+12del23 
    #  start=no offset, end has +12 offset
    #
    # c.368_368+1delTg
    #
    # sometimes a count, e.g.:
    # c.375+1_375+18del18
    my ($pos1, $offset1, $pos2, $offset2, $sequence) = ($1, $2, $3, $4, $5);
    die "grabbed too much for $nt" unless $offset1 or $offset2;

#    die "parse error in $nt, $pos1 != $pos2" unless $pos1 == $pos2;
    # not so fast: 1114+1_1115-1delaggacacatcaactg
    my $size;
    if ($sequence =~ /^\d+$/) {
      # a count rather than the actual sequence, e.g. 
      # c.375+1_375+18del18
      $size = $sequence;
      $sequence = "N" x $size;
#      printf STDERR "debug %s\n", join " ", $nt, $size, $sequence;
    } elsif ($offset1 and $offset2) {
      # 2 offsets specified
      my $abs = abs($offset1 - $offset2);
      $size = $abs + 1;
    } else {
      # only 1 offset specified, e.g. 1198-8_1199delcactccagAA
      $size = length($sequence);
    }
    die "size mismatch in $nt" unless $size = length($sequence);
    $parsable = 1;
    $is_deletion = 1;
    $is_coding_intronic = 1;
    $start = $pos1;
    $end = $pos2;
    $reference_sequence = uc($sequence);
    $variant_sequence = INDEL_CHAR;
    $event_length = length($reference_sequence);
#    printf STDERR "  debug2 $reference_sequence\n";
  } elsif ($nt =~ /^(\d+)([\+\-]\d+)ins([acgt]+)$/i) {
    # intronic COSMIC insertion (single site)
    # c.600-2insA
    $parsable = 1;
    $is_insertion = 1;
    $is_coding_intronic = 1;
    $start = $end = $1;
    my $offset = $2;
    $variant_sequence = uc($3);
    $reference_sequence = INDEL_CHAR;
    $event_length = length($variant_sequence);
#  } elsif ($nt =~ /^(\d+)([\+\-]\d+)_(\d+)([\+\-]\d+)ins([acgt]+|\d+)$/i) {
  } elsif ($nt =~ /^(\d+)([\+\-]\d+)?_(\d+)([\+\-]\d+)?ins([acgt]+|\d+)$/i and
	   ($2 or $4)) {
    # intronic COSMIC insertion (range)
    # c.5383-4_5383-3insT
    # c.977+1_977+2insT
    # c.2737-3_2737-2insTT
    # 
    # or cases with an offset on only one edge:
    # c.2225_2225+1insA 
    my ($pos1, $offset1, $pos2, $offset2, $sequence) = ($1, $2, $3, $4, $5);
    # exceptions?
    if ($sequence =~ /^\d+$/) {
      # c.453+1_453+2ins10
      $sequence = "N" x $sequence;
    }
    $parsable = 1;
    $is_insertion = 1;
    $is_coding_intronic = 1;
    if ($pos1 == $pos2) {
      $start = $pos1;
      $end = $pos2;
    } else {
      printf STDERR "WARNING: possible parse error in $nt, pos1 $pos1 != pos2 $pos2, using $pos1\n";
      $start = $end = $pos1;
    }
    $reference_sequence = INDEL_CHAR;
    $variant_sequence = uc($sequence);
    $event_length = length($sequence);
  } elsif ($nt =~ /^(\d+)([\+\-]\d+)?_(\d+)([\+\-]\d+)?([acgt]+)>([acgt]+)$/i
	   and ($2 or $4)) {
    # c.3158-1_3158GA>CT
    my ($pos1, $offset1, $pos2, $offset2, $seq_ref, $seq_var) = ($1, $2, $3, $4, $5, $6);
    $is_coding_intronic = 1;
    $parsable = 1;
    $start = $pos1;
    $end = $pos2;
    $reference_sequence = uc($seq_ref);
    $variant_sequence = uc($seq_var);
    if (length($seq_ref) == length($seq_var)) {
      # substitution, e.g. c.3158-1_3158GA>CT
      $is_substitution = 1;
      $event_length = length $reference_sequence;
    } else {
      # e.g. 
      # 489-2_495AGAAGGAAG>AGAAG 
      # is this maybe just a deletion?
      # AGAAGG survives, GAAG deleted?
      $is_complex_indel = 1;
      $event_length = -1;
      $complex_deletion_length = length($reference_sequence);
      $complex_insertion_length = length($variant_sequence);
    }
#    die join ",", $nt, $start, $end, $reference_sequence, $variant_sequence;
  } elsif ($nt =~ /^(\d+)_(\d+)([A-Z]+)>([A-Z]+)$/) {
    # substitution: either di/tri etc. or complex
    # c.549_550CC>TT
    ($start, $end, $reference_sequence, $variant_sequence) = ($1, $2, $3, $4);
    if (length($reference_sequence) == length($variant_sequence)) {
      # multi-base substitution
      $parsable = 1;
      $is_substitution = 1;
      $event_length = length($reference_sequence);
    } else {
      # complex, e.g. c.3919_3921ATA>TT (3 deleted, 2 inserted)
      $parsable = 1;
      $is_complex_indel = 1;
      $event_length = -1;
      $complex_deletion_length = $end - $start + 1;
#      confess "ref deletion length mismatch in $nt" unless $complex_deletion_length == length($reference_sequence);
      cluck "ref deletion length mismatch in $nt" unless $complex_deletion_length == length($reference_sequence);
      # sigh: 168_172AGCT>TA
      $complex_insertion_length = length($variant_sequence);
    }
  } elsif ($nt =~ /^(\d+)_(\d+)>([A-Z]+)$/) {
    # complex indel, e.g. COSMIC c.341_376>TC
    ($start, $end, $variant_sequence) = ($1, $2, $3);
    $reference_sequence = REF_SEQ_UNSPECIFIED;
    $parsable = 1;
    $is_complex_indel = 1;
    $event_length = -1;
    $complex_deletion_length = $end - $start + 1;
    $complex_insertion_length = length($variant_sequence);
  } elsif ($nt =~ /(\d+)([A-Z])>([A-Z]+)$/) {
    # complex indel, e.g. COSMIC c.616A>GG
    ($start, $end, $reference_sequence, $variant_sequence) = ($1, $1, $2, $3);
    $parsable = 1;
    $is_complex_indel = 1;
    $event_length = -1;
    $complex_insertion_length = length($variant_sequence);
    $complex_deletion_length = length($reference_sequence);
    die "length mismatch $nt $start $end $complex_insertion_length $complex_deletion_length" unless ($start + $complex_deletion_length) == ($end + 1);
  } elsif ($nt =~ /^(\d+)del(\d+)$/) {
    # TP53: g.7578206del1
    $start = $end = $1;
    $event_length = $2;
    $is_deletion = 1;
    confess "multi-base event in $nt" unless $event_length == 1;
    # verify start/end for for multi-base
    $parsable = 1;
  } elsif ($nt =~ /^(\d+)_(\d+)delins([AaCcGgTt]+)$/) {
    # TP53: g.7578531_7578532delinsG
    # 7577641_7577651delinstgggctctggg
    # complex, inserted bases specified
    ($start, $end, $variant_sequence) = ($1, $2, $3);
    $parsable = 1;
    my $el = $end - $start + 1;
    if ($el == length($variant_sequence)) {
      # MNV, e.g. COSMIC g.156085052_156085053delinsTT
      $is_substitution = 1;
      $event_length = $el;
      $reference_sequence = "N" x $event_length;
    } else {
      $is_complex_indel = 1;
      $event_length = -1;
      $complex_insertion_length = length($variant_sequence);
      $complex_deletion_length = ($end - $start) + 1;
    }
  } elsif ($nt =~ /^(\d+)_(\d+)delins(\d+)$/) {
    # TP53: g.7579366_7579371delins4, g.7578204_7578212delins3
    # complex, just length of inserted bases specified
    ($start, $end, $complex_insertion_length) = ($1, $2, $3);
    $parsable = 1;
    $is_complex_indel = 1;
    $event_length = -1;
    $complex_deletion_length = ($end - $start) + 1;
    $reference_sequence = "N" x $complex_deletion_length;
    $variant_sequence = "N" x $complex_insertion_length;
  } elsif ($nt =~ /^(\d+)delins(\d+)$/) {
    # TP53: g.7578230delins2
    # complex, single base deletion, just count of inserted bases specified
    ($start, $end, $complex_insertion_length) = ($1, $1, $2);
    $parsable = 1;
    $is_complex_indel = 1;
    $event_length = -1;
    $complex_deletion_length = 1;
    $reference_sequence = "N" x $complex_deletion_length;
    $variant_sequence = "N" x $complex_insertion_length;
  } elsif ($nt =~ /^(\d+)delins([ACGT]+)$/) {
    # TP53: g.7577093delinsGGG
    # complex, single base deletion, inserted bases specified
    ($start, $end, $variant_sequence) = ($1, $1, $2);
    $complex_insertion_length = length($variant_sequence);
    $parsable = 1;
    $is_complex_indel = 1;
    $event_length = -1;
    $complex_deletion_length = 1;
  } elsif ($nt =~ /^(\d+)_(\d+)dup$/) {
    # TP53: g.7579398_7579439dup
    # duplicated region, see http://www.hgvs.org/mutnomen/disc.html
#    ($start, $end) = ($1, $2);
    ($start, $end) = ($2, $2);
    # i.e. an insertion at the end of the range
    $parsable = 1;
    $is_insertion = 1;
    $event_length = ($end - $start) + 1;
  } elsif ($nt =~ /^(\d+)dup([ACGT]+)$/) {
    # ClinVar: g.43606600dupG
    $parsable = 1;
    $is_insertion = 1;
    $event_length = length $2;
    $start = $end = $1;
    $variant_sequence = $2;
  } elsif ($nt =~ /^(\d+)_(\d+)dup([ACGT]+)$/) {
    # ClinVar: NC_000016.9:g.2129049_2129056dupCTCACCAG
    # duplicate region
    $parsable = 1;
    $is_insertion = 1;
#    ($start, $end, $variant_sequence) = ($1, $2, $3);
    ($start, $end, $variant_sequence) = ($2, $2, $3);
    $event_length = length $variant_sequence;
  } elsif ($nt =~ /^(\d+)_(\d+)del([ACGT]+)ins([ACGT]+)$/) {
    # ClinVar: g.43609944_43609945delGCinsCG
    # complex indel where both sequences are specified
    # note: if event lengths are the same this is a 
    # multi-nucleotide substitution rather than a complex event
    # (as above)
    my ($ref_del, $ref_ins);
    ($start, $end, $ref_del, $ref_ins) = ($1, $2, $3, $4);
    $parsable = 1;
    $reference_sequence = $ref_del;
    $variant_sequence = $ref_ins;
    if (length($reference_sequence) == length($variant_sequence)) {
      # MNV
      $is_substitution = 1;
      $event_length = length($reference_sequence);
    } else {
      # complex
      $is_complex_indel = 1;
      $event_length = -1;
      $complex_deletion_length = length($ref_del);
      $complex_insertion_length = length($ref_ins);
    }
  } elsif ($nt =~ /^(\d+)del([ACGT]+)ins([ACGT]+)$/) {
    # ClinVar:
    # complex indel where both sequences are specified, single base version
    # NC_000013.10:g.32907469delCinsAA
    my ($ref_del, $ref_ins);
    ($start, $end, $ref_del, $ref_ins) = ($1, $1, $2, $3);
    $parsable = 1;
    $is_complex_indel = 1;
    $event_length = -1;
    $complex_deletion_length = length($ref_del);
    $complex_insertion_length = length($ref_ins);
    $reference_sequence = $ref_del;
    $variant_sequence = $ref_ins;
  } elsif ($nt =~ /^(\d+)_(\d+)del(\d+)ins([ACGT]+)$/) {
    # ClinVar: g.41209098_41209137del40insGA
    # complex indel where deletion size is specified and
    # inserted bases are specified
    my $inserted_not_saved;
    ($start, $end, $complex_deletion_length, $inserted_not_saved) = ($1, $2, $3, $4);
    $parsable = 1;
    my $el = $end - $start + 1;
    if ($el == length($inserted_not_saved)) {
      # MNV, e.g. COSMIC g.153296095_153296115del21insTGCTCAAGTCCTGGGGCTCAG
      $is_substitution = 1;
      $event_length = $el;
      $reference_sequence = "N" x $event_length;
      $variant_sequence = $inserted_not_saved;
      $complex_deletion_length = 0;
    } else {
      $is_complex_indel = 1;
      $event_length = -1;
      # FIX ME: 
      # might be able to "rescue" these by pulling in genomic sequence.
      $complex_insertion_length = length $inserted_not_saved;
    }
  } elsif ($nt =~ /^(\d+)_(\d+)dup(\d+)$/) {
    # duplicated region
    # e.g. ClinVar 88476172_88476195dup24
    $parsable = 1;
    $is_insertion = 1;
    $start = $end = $2;
    # duplicated region, so consider site after the first copy
    $event_length = $3;
    $variant_sequence = "N" x $event_length;
    my $elen = $2 - $1 + 1;
    confess "event length mismatch parsing $nt" unless $elen == $event_length;
  } elsif ($nt =~ /^(\d+)_(\d+)([ACGT]+)$/) {
    # range replaced with specified sequence
    $parsable = 1;
    $is_complex_indel = 1;
    $start = $1;
    $end = $2;
    $variant_sequence = $3;
    $complex_deletion_length = $2 - $1 + 1;
    $complex_insertion_length = length($variant_sequence);
  } elsif ($nt =~ /^(\d+)del(\d+)ins(\d+)$/) {
    # g.7577070del1ins8
    $parsable = 1;
    $is_complex_indel = 1;
    ($start, $complex_deletion_length, $complex_insertion_length) = ($1, $2, $3);
    $end = $start + ($complex_deletion_length - 1);
    $reference_sequence = "N" x $complex_deletion_length;
    $variant_sequence = "N" x $complex_insertion_length;
  } elsif ($nt =~ /^(\d+)_(\d+)del(\d+)ins(\d+)$/) {
    # g.7578390_7578392del3ins1
    $parsable = 1;
    $is_complex_indel = 1;
    ($start, $end, $complex_deletion_length, $complex_insertion_length) = ($1, $2, $3, $4);
    $reference_sequence = "N" x $complex_deletion_length;
    $variant_sequence = "N" x $complex_insertion_length;
  } elsif ($nt =~ /^(\d+)ins([ACGT]+)$/) {
    # simple insertion, e.g. Alsod: c.1527insTGGC
    $is_insertion = 1;
    $start = $end = $1;
    $reference_sequence = "";
    $variant_sequence = $2;
    $event_length = length($variant_sequence);
    $parsable = 1;
  } else {
    printf STDERR "WARNING: can't parse nt %s\n", $nt;
  }

#  if ($parsable and $is_substitution and my $strand = $self->auto_strand_fix) {
  if ($parsable and my $strand = $self->auto_strand_fix) {
    # apply to e.g. insertions too
    if ($strand eq "+") {
      # nothing to do
    } elsif ($strand eq "-") {
      $reference_sequence = reverse_complement($reference_sequence) unless $reference_sequence eq REF_SEQ_UNSPECIFIED;
      $variant_sequence = reverse_complement($variant_sequence);
    } else {
      die "unknown strand code $strand";
    }
  }
  
  $self->is_parsable($parsable);
  $self->is_substitution($is_substitution);
  $self->is_complex_indel($is_complex_indel);
  $self->complex_insertion_length($complex_insertion_length);
  $self->complex_deletion_length($complex_deletion_length);
  $self->is_insertion($is_insertion);
  $self->is_deletion($is_deletion);
  $self->is_coding_intronic($is_coding_intronic);
  $self->reference_sequence($reference_sequence);
  $self->variant_sequence($variant_sequence);
  $self->event_length($event_length);
  $self->start($start);
  $self->end($end);

  return $parsable;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
