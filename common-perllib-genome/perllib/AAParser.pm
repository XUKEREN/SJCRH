package AAParser;
# amino acid annotation parser, mainly intended to extract affected codon #s
# COSMIC, SJ, etc.
# MNE 9/2013
# TO DO: hgvs.org recommendations, e.g. 
# http://www.hgvs.org/mutnomen/recs-DNA.html

use strict;

use Carp qw(confess);
use Bio::Tools::CodonTable;

use Configurable;
use Exporter;

use constant CODON_VARIANT_FS => "fs";
use constant CODON_VARIANT_DELETION => "-";
# for use in indel annotations
use constant CODON_UNKNOWN_AA => "X";

@AAParser::ISA = qw(Configurable Exporter);
@AAParser::EXPORT_OK = qw();

use MethodMaker qw(
	aa_raw
codon_start
codon_end
cooked_substitution
is_silent
is_missense
is_nonsense

codon_reference
codon_number
codon_variant

parse_x_as_stop
parsable
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->parse_x_as_stop(1);
  $self->configure(%options);
  return $self;
}

sub parse {
  #
  #  rough parser for codon numbers only
  #
  my ($self, $aa, %options) = @_;

  my $codon_reference = "";
  my $codon_variant = "";

  my $aa_raw = $aa;
  $self->aa_raw($aa_raw);

  if ($aa =~ s/^p([A-Z]\d+)/p.$1/) {
    # COSMIC: repair
    # pK120_F121insLYHKGLLK
  }

  $aa =~ s/^\[//;
  # IARC: [p.V225_Q331del]

#  $aa =~ s/^p\.//;
  # COSMIC/ASU/IARC
#  $aa =~ s/^[pg][\.\,]//;
  $aa =~ s/^[pg][\.\,]\s*//;
  # allow for some flaky formatting e.g. in RB1 database
  # p. Glu748Lysfs*8

  if ($aa =~ /\(/) {
    # RB1: p.(Leu134*)
    $aa =~ s/\((.*)\)$/$1/;
  }

  $aa =~ s/\/.*//;
  # RB1: p.Arg46Lys/splice


  my ($codon_start, $codon_end);

  #
  #  generic and COSMIC-specific examples:
  #
  if ($aa =~ /^[A-Z](\d+)[A-Z]$/) {
    # standard: p.R245C
    $codon_start = $codon_end = $1;
    # should call parse_substitution() instead
#  } elsif ($aa =~ /^\w{3}(\d+)\w{3}$/) {
# FAIL: parses p.Ser2201X as codon 22!
  } elsif ($aa =~ /^[A-z]{3}(\d+)[A-z]{3}$/) {
    # standard: Trp158Ter
    $codon_start = $codon_end = $1;
#  } elsif ($aa =~ /^(\d+)\w{3}$/) {
  } elsif ($aa =~ /^([A-Z][a-z]{2})(\d+)_([A-Z][a-z]{2})(\d+)del\1([A-Z][a-z]{2})*\3$/) {
    # in-frame deletion of multiple codons
    # e.g. NHGRI p.Lys1767_Asp1769delLysLeuAsp
    # $1 = Lys
    # $2 = 1767
    # $3 = Asp
    # $4 = 1768
    # ...this regexp will work with multiple codon deletions
    # so long as 1st and last are consistent
    ($codon_start, $codon_end) = ($2, $4);
  } elsif ($aa =~ /^(\d+)[A-z]{3}$/) {
    printf STDERR "AAP: %s => %d\n", $aa, $1;
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^[A-Z](\d+)\*$/) {
    # p.R266*
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^([A-Z])(\d+)fs/) {
    # p.E409fs*17
    $codon_reference = $1;
    $codon_variant = CODON_VARIANT_FS;
    $codon_start = $codon_end = $2;
  } elsif ($aa =~ /^([A-Z])(\d+)FS\*\d+$/) {
    # COSMIC: uppercased p.N473FS*8
    $codon_reference = $1;
    $codon_start = $codon_end = $2;
    $codon_variant = CODON_VARIANT_FS;
  } elsif ($aa =~ /^[A-Z]+(\d+)_[A-Z](\d+)(ins|del)([A-Z]+)/) {
    # p.P1331_A1332insTP
    # p.L215_F219delLSRLF
    $codon_start = $1;
    $codon_end = $2;
    my $type = $3;
    my $thing = $4;
    if ($type eq "del") {
      $codon_reference = $thing;
      $codon_variant = CODON_VARIANT_DELETION;
    } elsif ($type eq "ins") {
      $codon_reference = CODON_VARIANT_DELETION;
      $codon_variant = $thing;
    } else {
      die;
    }

  } elsif ($aa =~ /^([A-Z])(\d+)del([A-Z])?/) {
    # p.F1002delF
    # p.S368del
    $codon_reference = $1;
    $codon_start = $codon_end = $2;
    my $result = $3;
    if ($result) {
      if ($result eq $codon_reference) {
	# p.F1002delF
	$codon_variant = CODON_VARIANT_DELETION;
      } else {
	die "WTF $aa";
      }
    } else {
      # V2304del
      $codon_variant = CODON_VARIANT_DELETION;
    }
  } elsif ($aa =~ /^([A-Z])(\d+)>(.*)/) {
    # p.L583>?
    # p.K294>RGG 
    $codon_reference = $1;
    $codon_start = $codon_end = $2;
    $codon_variant = $3;
#  } elsif ($aa =~ /^[A-Z](\d+)_[A-Z](\d+)>/) {
  } elsif ($aa =~ /^([A-Z])(\d+)_([A-Z])(\d+)>(.*)/) {
    # p.T751_I759>REA
    # p.K666_L667>M 
    my $p_start = $1;
    $codon_start = $2;
    my $p_end = $3;
    $codon_end = $4;
    $codon_variant = length($5) ? $5 : CODON_VARIANT_DELETION;

    $codon_reference = generate_reference_aa($p_start, $codon_start, $p_end, $codon_end);
  } elsif ($aa =~ /^\*(\d+)[A-Z]/) {
    # p.*667Y
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^\*(\d+)\*$/) {
    # *982*
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^\*(\d+)fs/) {
    # p.*417fs
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^\([A-Z]?(\d+)\)fs/) {
    # p.(Y346)fs*
    # p.(2287)fs
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^(\d+)_(\d+)del(\d+)?/) {
    # p.1670_1673del
    # p.982_1028del47
    $codon_start = $1;
    $codon_end = $2;
    my $length = $3;
    my $len1 = $codon_end - $codon_start + 1;
    if ($length) {
      die "AA deletion length sanity fail: compare $aa len $length $len1" unless $length == $len1;
    } else {
      $length = $len1;
    }
    $codon_reference = CODON_UNKNOWN_AA x $length;
    $codon_variant = CODON_VARIANT_DELETION;
  } elsif ($aa =~ /^[A-Z](\d+)_[A-Z](\d+)$/) {
    # p.N260_I262
    $codon_start = $1;
    $codon_end = $2;
  } elsif ($aa =~ /^(\d+)_(\d+)ins([A-Z]+)$/) {
    # p.120_121insLYHKGLLK
    # p.794_795insD
    $codon_start = $1;
    $codon_end = $2;
    $codon_reference = CODON_VARIANT_DELETION;
    $codon_variant = $3;
  } elsif ($aa =~ /^[A-Z](\d+)$/) {
    # p.R2468
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^[A-Z](\d+)_E\d+splice/) {
    # L314_E6splice_region (SJ)
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^[A-Z]ins[A-Z](\d+)$/) {
    # GinsP79 (SJ)
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^[A-Z](\d+)_[A-Z](\d+)fs/) {
    # Q225_L226fs (SJ)
    $codon_start = $1;
    $codon_end = $2;
  } elsif ($aa =~ /^[A-Z]{3}(\d+)del$/i) {
    # ASU: Glu441del (in-frame deletion)
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^([A-Z][a-z]{2})(\d+)del\1$/) {
    # NHGRI: p.Ala1693delAla (in-frame deletion)
    $codon_start = $codon_end = $2;
  } elsif ($aa =~ /^[A-Z]{3}(\d+)_[A-Z]{3}(\d+)del$/i) {
    # ASU: p.Leu862_Leu884del
    $codon_start = $1;
    $codon_end = $2;
  } elsif ($aa =~ /^[A-z]{3}(\d+)[A-z]{3}fs/) {
    # ASU: p.Pro112ProfsX16
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^[A-z]{3}(\d+)Term$/) {
    # HGMD: Trp158Term
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^Term(\d+)[A-z]{3}$/) {
    # HGMD: Term378Gln
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^[A-z]{3}(\d+)\*$/) {
    # Arg453*
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^[A-z]{3}(\d+)stop$/i) {
    # NHLBI: Trp736Stop
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^[A-z]{3}(\d+)X$/i) {
    # umd.be: apparently used for both stops and complex cases, e.g.
    # c.37_44del / p.Glu13X
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^[A-Z]{3}(\d+)_[A-Z]{3}(\d+)(del|ins)/i) {
    # umd.be: p.Asp1692_Ala1693insPheIle
    # NOTE: should be low-priority
    $codon_start = $1;
    $codon_end = $2;
#    printf STDERR "caught delins: %s %d-%d %s\n", $aa, $codon_start, $codon_end, $3;
  } elsif ($aa =~ /^[A-Z][a-z]{2}(\d+)$/) {
    # umd.be: p.Thr3085
    # NOTE: should be low-priority
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^[A-Z][a-z]{2}(\d+)delins/) {
    # umd.be: 
    #   p.Asp1692delinsAspIle
    #   p.Ile1855delinsMetPhe
    $codon_start = $codon_end = $1;
#    printf STDERR "caught delins %s %d\n", $aa, $codon_start;
  } elsif ($aa =~ /^(\d+)delins/ and $1) {
    # umd.be: p.1864delins
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^[A-Z][a-z]{2}(\d+)fsX\d+/) {
    # RB1: p.Lys765fsX44
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^[A-Z][a-z]{2}(\d+)[\?\=]fs$/) {
    # NHGRI: p.Val11?fs, p.Gly57=fs
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^[A-Z][a-z]{2}(\d+)_[A-Z][a-z]{2}(\d+)=fs$/) {
    # NHGRI: p.Gln895_Ser896=fs
    ($codon_start, $codon_end) = ($1, $2);
  } elsif ($aa =~ /^[A-Z][a-z]{2}(\d+)_[A-Z][a-z]{2}(\d+)\?(fs)?$/) {
    # NHGRI:
    # - p.Arg7_Glu10?
    # - p.Arg7_Asn16?fs
    ($codon_start, $codon_end) = ($1, $2);
  } elsif ($aa =~ /^[A-Z][a-z]{2}(\d+)_[A-Z][a-z]{2}(\d+)[A-Z][a-z]{2}[A-Z][a-z]{2}fs$/) {
    # NHGRI: p.Leu22_Glu23LeuValfs
    # ask Sheila: how to interpret this?
    ($codon_start, $codon_end) = ($1, $2);
  } elsif ($aa =~ /^(\d+)_(\d+)dup$/) {
    # UMD APC: p.560_559dup
    ($codon_start, $codon_end) = sort {$a <=> $b} ($1, $2);
  } elsif ($aa =~ /^(\d+)_[A-Z][a-z]{2}(\d+)dup$/) {
    # UMD APC: p.729_Asn728dup
    ($codon_start, $codon_end) = sort {$a <=> $b} ($1, $2);
  } elsif ($aa =~ /^[A-Z][a-z]{2}(\d+)\?$/) {
    # NHGRI: p.Phe15?
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^\*(\d+)\?$/) {
    # COSMIC: p.*553?
    $codon_start = $codon_end = $1;
  } elsif ($aa =~ /^[A-Z](\d+)FS\*\d+$/) {
    # COSMIC: p.N473FS*8
    $codon_start=  $codon_end = $1;
#  } elsif ($aa =~ /^(\d+)_(\d+)del(\d+)$/) {
    # COSMIC: p.982_1028del47
#    die "hey now $aa";
  } elsif ($aa =~ /^([A-Z])(\d+)_([A-Z])(\d+)del(.*)/) {
    # COSMIC: p.A344_A348del
    my $p_start = $1;
    $codon_start = $2;
    my $p_end = $3;
    $codon_end = $4;
    my $thing = $5;
    if ($thing) {
      if ($thing =~ /^\d+$/) {
	# G1804_H1815del12
	my $slen = $codon_end - $codon_start + 1;
	die "length sanity fail for $aa" unless $slen == $thing;
      } else {
	printf STDERR "ERROR: unhandled possible complex deletion %s!\n", $aa;
      }
    }

    $codon_reference = generate_reference_aa($p_start, $codon_start, $p_end, $codon_end);
    $codon_variant = CODON_VARIANT_DELETION;

#    die join ",", $aa, $codon_start, $codon_end, $codon_reference, $codon_variant;
  } elsif ($aa =~ /^([A-Z])(\d+)_([A-Z])(\d+)ins(\S+)/) {
    # COSMIC: p.N307_V308ins16
    my $p_start = $1;
    $codon_start = $2;
    my $p_end = $3;
    $codon_end = $4;
    my $thing = $5;
    die "long insertion site for $aa" unless $codon_start == $codon_end - 1;
    $codon_reference = CODON_VARIANT_DELETION;
    if ($thing =~ /^\d+$/) {
      # length specified, e.g. N307_V308ins16
      $codon_variant = CODON_UNKNOWN_AA x $thing;
    } else {
      # R1386_C1387ins*
      $codon_variant = $thing;
    }
#    die join ",", $aa, $codon_start, $codon_end, $codon_reference, $codon_variant;

  } elsif ($aa =~ /^\S+ins[A-Z]+(\d+)$/) {
    # SJ post, e.g. VAinsG3519
    $codon_start = $codon_end = $1;
  } elsif (0 and
	   # DISABLED, see below
	   $aa =~ /\-(\d+)_[A-Z][a-z]{2}(\d+)del/) {
    # umd.be:
    # nt: c.-200_441del aa: p.-66_Leu147del
    # nt: c.-200_5406del aa: p.-66_Thr1802del
    # confusing: codon sizes match, but meaning of -200?
    # seems suspicious that both starts are -200/-66
    $codon_start = $1;
    $codon_end = $2;
    printf STDERR "caught garbage %s %d-%d\n", $aa, $codon_start, $codon_end;
  }

  #
  #  SJ-specific examples:
  #
  my $parsable = 0;
  if ($codon_start and $codon_end) {
    $parsable = 1;
    confess "codon order fail $codon_start $codon_end in $aa" if $codon_start > $codon_end;
  }

  $self->parsable($parsable);
  $self->codon_start($codon_start);
  $self->codon_end($codon_end);
  $self->codon_reference($codon_reference);
  $self->codon_variant($codon_variant);

  return $parsable;
}

sub parse_substitution {
  #
  # parse substitutions which follow standard formatting:
  #
  my ($self, $thing, %options) = @_;
  my $thing_raw = $thing;
  $thing =~ s/^[pg][\.\,]//;
  # RB1: g.Lys329X (bad formatting)

  $thing =~ s/\/.*//;
  # RB1: p.Arg46Lys/splice

  if ($thing =~ /\(/) {
    # RB1: p.(Leu134*)
    $thing =~ s/\((.*)\)$/$1/;
  }
  
  $thing =~ s/%3D$/=/;
  # escaped "=" value from VEP for silent variant,
  # e.g. p.Ser286%3D => p.Ser286=

  my ($aa1_raw, $codon_number, $aa2_raw);
#  if ($thing =~ /^(\w{3})(\d+)(\w{3})$/) {
# FAIL: parses V12693D as V12/6/93D!
  if ($thing =~ /^([A-Z]{3})(\d+)([A-Z]{3})$/i) {
    # p.Pro33Ser
    ($aa1_raw, $codon_number, $aa2_raw) = ($1, $2, $3);
#  } elsif ($thing =~ /^(\w)(\d+)(\w)$/) {
  } elsif ($thing =~ /^([A-Z])(\d+)([A-Z])$/i) {
    # p.P33S
    ($aa1_raw, $codon_number, $aa2_raw) = ($1, $2, $3);
#  } elsif ($thing =~ /^([A-z]+)(\d+)(\w+)$/) {
#  } elsif ($thing =~ /^(\w+?)(\d+)(\w+)$/) {
  } elsif ($thing =~ /^([A-Z]{3})(\d+)([A-Z])$/i) {
    # mix of long and short, e.g. p.Glu30X
    # in RB1 these are sometimes actually stop codons  :/
    ($aa1_raw, $codon_number, $aa2_raw) = ($1, $2, $3);
#  } elsif ($thing =~ /^([\w\*]+?)(\d+)([\w\*]+)$/) {
  } elsif ($thing =~ /^([A-Z][a-z]{2})(\d+)\=$/) {
    # silent, e.g. NHGRI p.Val14=
    ($aa1_raw, $codon_number, $aa2_raw) = ($1, $2, $1);
  } elsif ($thing =~ /^([A-z\*]+?)(\d+)([A-z\*]+)$/) {
    # special handling for stop codons
    ($aa1_raw, $codon_number, $aa2_raw) = ($1, $2, $3);
    foreach ($aa1_raw, $aa2_raw) {
      $_ = "Ter" if lc($_) eq "term" or $_ eq "*" or lc($_) eq "stop";
      # HGMD: Trp158Term, Term378Gln
      # SJ: R453*
      # NHLBI: Trp736Stop
    }
  }

  my $parsable = 0;
  my $cooked_substitution = "";
  my $is_silent = 0;
  my $is_missense = 0;
  my $is_nonsense = 0;

  if ($aa1_raw and $codon_number and $aa2_raw) {
    my $ct = new Bio::Tools::CodonTable();

    my @cooked;
    my $idx = 0;
    foreach my $raw ($aa1_raw, $aa2_raw) {
      if (uc($raw) eq "U" or uc($raw) eq "SEC") {
	# revtranslate() won't work here
	$cooked[$idx] = "U";
      } elsif (uc($raw) eq "X") {
	# revtranslate() won't work
	if ($self->parse_x_as_stop()) {
	  # umd.be: p.Glu9X
	  # ASU TERT: c.2665C>T p.Arg889X, manually confirmed stop
	  $cooked[$idx] = "*";
	} else {
	  $cooked[$idx] = "X";
	}
      } elsif (length($raw) == 1) {
	# already uses short format
	my @codons = $ct->revtranslate($raw);
	die "can't revtranslate $raw in $thing_raw" unless @codons;
	# trust but verify
	$cooked[$idx] = $raw;
      } elsif (my @codons = $ct->revtranslate($raw)) {
	$cooked[$idx] = $ct->translate($codons[0]);
      }
      $idx++;
    }

    if ($cooked[0] and $cooked[1]) {
      $cooked_substitution = sprintf '%s%d%s',
      $cooked[0], $codon_number, $cooked[1];
      $self->codon_reference($cooked[0]);
      $self->codon_number($codon_number);
      $self->codon_variant($cooked[1]);
#      printf STDERR "cooked=%s\n", $cooked_substitution;
      if ($cooked[0] eq $cooked[1]) {
	$is_silent = 1;
      } elsif ($cooked[1] eq "*") {
	$is_nonsense = 1;
      } else {
	$is_missense = 1;
      }
      $parsable = 1;
    }
  }

  $self->is_silent($is_silent);
  $self->is_missense($is_missense);
  $self->is_nonsense($is_nonsense);

  return $self->cooked_substitution($cooked_substitution);
}

sub get_class_string {
  my ($self, $aa) = @_;
  my $class = "";
  if ($self->parse_substitution($aa)) {
    if ($self->is_silent()) {
      $class = "silent";
    } elsif ($self->is_missense) {
      $class = "missense";
    } elsif ($self->is_nonsense) {
      $class = "nonsense";
    } else {
      die;
    }
  }
  return $class;
}

sub generate_reference_aa {
  my ($p_start, $codon_start, $p_end, $codon_end) = @_;
#  $codon_reference = generate_reference_aa($p_start, $codon_start, $p_end, $codon_end);
  die unless $codon_start < $codon_end;
  my $reference_aa;
  if ($codon_start == $codon_end - 1) {
    $reference_aa = $p_start . $p_end;
  } else {
    my $len = $codon_end - $codon_start + 1;
    my @aa = split //, CODON_UNKNOWN_AA x $len;
    $aa[0] = $p_start;
    $aa[$#aa] = $p_end;
    $reference_aa = join "", @aa;
    # from annotation:
    # - create reference AA of appropriate length
    # - populate start/end codon
    # - use Xs for other codons as we don't have sufficient information
  }
  return $reference_aa;
}

sub condense_deletion {
  # Digest a variant touching 2 codons that does not alter the first
  # into a simpler format referring to only the altered codon.
  # e.g. convert K159_E160>K into E160del
  # This 
  # convert a raw substition that touches but does not change a codon.
  my ($self) = @_;
  my $result;
  my $cref = $self->codon_reference();
  my $cvar = $self->codon_variant();

  if ($self->parsable) {
    if (($self->codon_start() == $self->codon_end() - 1) and
	# 2-codon event
	length($cref) > length($cvar) and
	# deletion
	substr($cref, 0, 1) eq $cvar
	) {
      $result = sprintf "%s%ddel", substr($cref, 1), $self->codon_end();
    }
  }
  return $result;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
