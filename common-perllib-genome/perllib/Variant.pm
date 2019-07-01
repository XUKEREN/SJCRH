package Variant;
# attempt to standardize/cook variants (SNV/MNV/indel) into a standard format
# MNE 5/2015

use strict;
use Carp qw(confess);

use Configurable;
use Exporter;

use MiscUtils qw(dump_die);
use TabixFile;
use List::Util qw(min);

use constant INDEL_CHAR => "-";
use constant VCF_SYMBOLIC_DEL => "<DEL>";
# symbolic allele

@Variant::ISA = qw(Configurable Exporter);
@Variant::EXPORT_OK = qw(INDEL_CHAR new_variant_from_row);

use MethodMaker qw(
reference_name
start
end
reference_allele
variant_allele
is_substitution
is_snv
is_mnv
is_complex
is_insertion
is_deletion
event_length
additional_shift_bases
warning
raw_key
exception_warn
exception
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  return $self;
}

sub import_vcf_row {
  # import a single parsed VCF row with exactly one variant allele.
  # i.e. post VCFUtils::get_alt_rows(), or via VariantIteratorVCF.pm,
  # or "bcftools norm -m-both"
  my ($self, %options) = @_;
  my $row = $options{"-row"} || die "-row";

  my $f_chr = $options{"-f-chr"} || "CHROM";
  my $f_chr_alt = $options{"-f-chr-alt"} || "#CHROM";
  my $f_pos = $options{"-f-pos"} || "POS";
  my $f_ra = $options{"-f-ra"} || "REF";
  my $f_va = $options{"-f-va"} || "ALT";
  # default to VCF style, but can be overridden

  my $chrom_raw = $row->{$f_chr} || $row->{$f_chr_alt} || die "can't find $f_chr or $f_chr_alt";
  my $pos = $row->{$f_pos} || die $f_pos;
  my $ref_base = $row->{$f_ra} || die $f_ra;
  my $var_base = $row->{$f_va} || die $f_va;

  $self->raw_key(join ".", $chrom_raw, $pos, $ref_base, $var_base);

  my $chrom = cook_chrom($chrom_raw);

  my $is_symbolic;
  my $is_broken;
  if ($var_base =~ /\W/) {
    if ($var_base eq VCF_SYMBOLIC_DEL) {
      # ok
      $is_symbolic = 1;
    } elsif ($var_base =~ /^</) {
      # VCF tag
      $is_symbolic = 1;
    } elsif ($var_base eq ".") {
      # not provided
      # VCF spec: "If there are no alternative alleles, then the
      # missing value should be used."
      # i.e., this variant not found in this sample?
      $is_broken = "alt allele not provided";
    } else {
      $is_broken = "unhandled var base $var_base";
      # unparsed multiple alternate alleles, symbolic alleles, etc.
    }
  }


  $self->reference_name($chrom);
  my @warnings;

  if ($is_broken) {
    if ($self->exception_warn()) {
#      printf STDERR "exception: %s\n", $is_broken;
      $self->exception($is_broken);
    } else {
      confess $is_broken;
    }
  } elsif ($var_base eq VCF_SYMBOLIC_DEL) {
    # VCF 4.1 spec:
    # "If any of the ALT alleles is a symbolic allele (an angle-bracketed ID
    # String "<ID>") then the padding base is required and POS denotes the
    # coordinate of the base preceding the polymorphism."
    $self->is_deletion(1);
    die "unhandled multi base padding in symbol deletion" if length($ref_base) > 1;

    my $pad_length = 1;
    # required (spec)
    my $elen = 1;
    # not sure how this works for longer events
    $self->event_length($elen);
    $self->reference_allele("N");
    # not provided
    $self->variant_allele(INDEL_CHAR);
    $self->start($pos + $pad_length);
    $self->end($self->start() + $elen - 1);
    push @warnings, "symbolic_deletion";
  } elsif ($is_symbolic) {
    my $elen = length($ref_base);
    $self->event_length($elen);
    $self->reference_allele($ref_base);
    $self->variant_allele($var_base);
    my $pad_length = 0;
    # ?
    $self->start($pos + $pad_length);
    $self->end($self->start() + $elen - 1);
    push @warnings, "symbolic_alt_not_parsed";
  } elsif (length($ref_base) == length($var_base)) {
    # substitution (SNV, MNV, etc.)
    $self->is_substitution(1);
    my $rlen = length($ref_base);
    $self->is_mnv($rlen > 1 ? 1 : 0);
    $self->is_snv($rlen == 1 ? 1 : 0);

    my $pad_bases = 0;
    for (my $i = 0; $i < length($ref_base); $i++) {
      if (substr($ref_base, $i, 1) eq substr($var_base, $i, 1)) {
	$pad_bases++;
      } else {
	last;
      }
    }
    if ($pad_bases) {
      # e.g. chr3.10859757.TTTTTTTTTTTTTTA.TTTTTTTTTTTTTTT
      # demo_fake_mnv.vcf
      # reduces to SNV
      $ref_base = substr($ref_base, $pad_bases);
      $var_base = substr($var_base, $pad_bases);
      $self->additional_shift_bases($pad_bases);
      push @warnings, "reduced_substitution";
      $pos += $pad_bases;
    }

    my $pad_bases_r = 0;
    for (my $i = length($ref_base) - 1; $i >= 0; $i--) {
      if (substr($ref_base, $i, 1) eq substr($var_base, $i, 1)) {
	$pad_bases_r++;
      } else {
	last;
      }
    }
    if ($pad_bases_r) {
      # e.g. chr5.58638410.GG.TG, demo_broken_mnv_trailing_g.vcf
      # out of spec??
      # this is actually a SNV, trailing G is duplicated and irrelevant
      $ref_base = substr($ref_base, 0, length($ref_base) - $pad_bases_r);
      $var_base = substr($var_base, 0, length($var_base) - $pad_bases_r);
      push @warnings, "reduced_substitution_trailing";
      # out of spec???
      # trailing bases stripped so position doesn't change
    }

    push @warnings, "double_end_trim", if ($pad_bases and $pad_bases_r);
    
    $self->start($pos);
    $self->end($pos + length($var_base) - 1);
    $self->reference_allele($ref_base);
    $self->variant_allele($var_base);
    $self->event_length(length($ref_base));

    if ($ref_base eq "" and $var_base eq "") {
      $self->handle_exception("identical ref and alt alleles");
    }
  } elsif (length($ref_base) < length($var_base) and
	   index($var_base, $ref_base) == 0
#	   and length($ref_base) == length($var_base) - 1
      ) {
    # simple insertion with padding base(s) provided.
    # spec calls for exactly one base of padding, however
    # sometimes it seems there may be more, e.g.
    # 1.768116.AGTTTT.AGTTTTGTTTT which has 6 bases of padding
    # instead of 1 (this site also seems to have ambiguous mapping)
#    dump_die($row, "simple insertion", 1);
    my $pad_length = length($ref_base);
    if ($pad_length > 1) {
      push @warnings, "extra_reference_padding_insertion";
      $self->additional_shift_bases($pad_length - 1);
      # padding of 1 is expected, only report additional
    }
    $self->is_insertion(1);
    my $ins_chunk = substr($var_base, $pad_length);
    my $elen = length($ins_chunk);
    $self->event_length($elen);
    $self->reference_allele(INDEL_CHAR);
    $self->variant_allele($ins_chunk);
#    $self->start($pos);
    $self->start($pos + $pad_length - 1);
    # going back and forth on this.
    # the deletion 
    # chr1.1067576.GGTGGGGTTGTGGGGCCTCTCAG.GGTGGGGTTGTGGGGCCTCTC 
    # REQUIRES the position be adjusted past the reference sequence.
    # This logic also seems most consistent with VCF spec re: REF field.
    # However for insertions this often seems to move the target 
    # site beyond that in the BAM! (always?)
    # - should the handling be different for insertions??
    # 1.1014316.C.CG: insertion happens after 1014316 and before 1014317
    # 3.37040271.C.CCTT: after 37040271 and before 37040272
    $self->end($self->start + 1);
    # hack (both before and after bases touched)
    # interbase event; maybe end should be same as start?
    # ????
#    dump_die($self);
  } elsif (length($ref_base) > length($var_base) and
	   index($ref_base, $var_base) == 0 
#	   and length($ref_base) == length($var_base) + 1
) {
    # simple deletion with padding base provided(s)
#    dump_die($row, "simple deletion", 1);
    my $pad_length = length($var_base);
    if ($pad_length > 1) {
      push @warnings, "extra_reference_padding_deletion";
      $self->additional_shift_bases($pad_length - 1);
    }
    $self->is_deletion(1);
    my $del_chunk = substr($ref_base, $pad_length);
    my $elen = length($del_chunk);
    $self->event_length($elen);
    $self->reference_allele($del_chunk);
    $self->variant_allele(INDEL_CHAR);
    $self->start($pos + $pad_length);
    $self->end($self->start() + $elen - 1);
#    dump_die($self);
  } else {
    #
    #  complex indel (or possible a simple indel after addt'l trimming):
    #

    #
    # trim leading padding bases:
    #
    my $max = min(length($ref_base), length($var_base));
    my $leading_pad_bases = 0;
    for (my $i = 0; $i < $max; $i++) {
      if (substr($ref_base, $i, 1) eq substr($var_base, $i, 1)) {
	$leading_pad_bases++;
      } else {
	last;
      }
    }
    $ref_base = substr($ref_base, $leading_pad_bases);
    $var_base = substr($var_base, $leading_pad_bases);

    #
    # trim trailing padding bases:
    #
    $max = min(length($ref_base), length($var_base));
    my $pad_bases_r = 0;
    for (my $lookback = 1; $lookback <= $max; $lookback++) {
      if (substr($ref_base, -$lookback) eq substr($var_base, -$lookback)) {
	$pad_bases_r++;
      } else {
	last;
      }
    }
    # NHLBI:
    # 1       1323144 rs35654872      CTG     C,CG  
    # #1: 1323145.TG.-
    # #2: 1323145.T- (trailing G, so not TG->G but rather T->-)
    if ($pad_bases_r) {
#      printf STDERR "before: ra=%s va=%s\n", $ref_base, $var_base;
      $ref_base = substr($ref_base, 0, length($ref_base) - $pad_bases_r);
      $var_base = substr($var_base, 0, length($var_base) - $pad_bases_r);
    }

    if ($ref_base and $var_base) {
      # event is still complex
      $self->is_complex(1);
      $self->event_length(-1);
    } elsif ($ref_base) {
      # event trimmed to simple deletion
      $self->is_deletion(1);
      $self->event_length(length $ref_base);
      $var_base = INDEL_CHAR;
      push @warnings, "complex_trimmed_to_simple_deletion";
    } elsif ($var_base) {
      # event trimmed to simple insertion
      $ref_base = INDEL_CHAR;
      $self->is_insertion(1);
      $self->event_length(length $var_base);
      # e.g. NHLBI
      # 1       1877103 .       CG      CTG,C 
      # #1: 1.1877103.-.T: the inserted T happens after 1877103, not 04!
      # #2: 1.1877104.G.-
#      dump_die($self, "test me: complex trimmed to simple insertion and multiple leading pad bases $ref_base $var_base") if $leading_pad_bases > 1;
      $leading_pad_bases--;
      # compensate for stripping padding (below) since this
      # resolves to a simple insertion
      push @warnings, "complex_trimmed_to_simple_insertion";
#      dump_die($self, "test me, complex trimmed to simple insertion");
    } else {
      die;
    }

    $self->reference_allele($ref_base);
    $self->variant_allele($var_base);
    $self->start($pos + $leading_pad_bases);
    $self->end($pos + length($ref_base) - 1 - $leading_pad_bases);
    $self->additional_shift_bases($leading_pad_bases - 1);
    if ($leading_pad_bases == 0) {
      # ever happens?
      push @warnings, "no_complex_leading_padding";
    } elsif ($leading_pad_bases > 1) {
      # definitely happens
      push @warnings, "complex_extra_leading_padding";
    }
  }
  $self->warning(join ",", @warnings);

}

sub import_bambino_row {
  # import raw Bambino call row
  my ($self, %options) = @_;
  my $row = $options{"-row"} || die "-row";
  my $chr = cook_chrom($row->{Chr} || confess "can't find Chr");

  my ($pos, $ra, $va);
  if ($options{"-postprocessed"}) {
    $pos = $row->{WU_HG19_Pos} || $row->{WU_HG18_Pos} || die "no WU_HG19_Pos/WU_HG18_Pos";
#    $ra = $row->{ReferenceAllele} || die;
#    $va = $row->{MutantAllele} || die;
    foreach my $f (qw(ReferenceAllele MutantAllele)) {
      dump_die($row, "where is $f") unless exists $row->{$f};
    }
#    $ra = $row->{ReferenceAllele} || die;
#    $va = $row->{MutantAllele} || die;
    # while SJ-post indeed provides a "-" for indels, this is probably
    # too strict for other use cases, e.g. bambino2annovar
    $ra = $row->{ReferenceAllele} || "";
    $va = $row->{MutantAllele} || "";
  } else {
    $pos = $row->{Pos} || dump_die($row, "no Pos");
    $ra = $row->{Chr_Allele};
    $va = $row->{Alternative_Allele};
  }
  foreach ($ra, $va) {
    $_ = "" if $_ eq "-";
    # convert to raw
  }

  $self->reference_name($chr);

  if ($ra eq "") {
    # insertion
    $self->is_insertion(1);
    $pos--;
    # for insertions Bambino reports the 1-based base AFTER the insertion.
    # Most people use the base before, so standardize there.
    $self->reference_allele(INDEL_CHAR);
    $self->variant_allele($va || die);
    $self->start($pos);
    $self->end($pos);
    $self->event_length(length $va);
  } elsif ($va eq "") {
    # deletion
#    dump_die($row, "del", 1);
    $self->is_deletion(1);
    my $elen = length($ra);
    $self->event_length($elen);
    $self->reference_allele($ra);
    $self->variant_allele(INDEL_CHAR);
    $self->start($pos);
    $self->end($self->start() + $elen - 1);
  } else {
    # substitution
    $self->is_substitution(1);
    my $rl = length $ra;
    my $vl = length $va;

    # MNVs and complex alleles can come from converted
    # non-Bambino output (e.g. vcf2tab.pl)
    if ($rl == $vl) {
      $self->is_substitution(1);
      my $rl = length($ra);
      $self->is_mnv($rl > 1 ? 1 : 0);
      $self->is_snv($rl == 1 ? 1 : 0);
    } else {
      $self->is_complex(1);
    }

#    dump_die($row, "length fail $ra $va") unless length($ra) == length($va) and length($ra) == 1;
    my $event_len = length($ra);
    $self->reference_allele($ra);
    $self->variant_allele($va);
    $self->event_length($event_len);
    $self->start($pos);
    $self->end($pos + $event_len - 1);
  }


}

sub get_key {
  my ($self) = @_;
  return join ".", $self->reference_name, $self->start, $self->reference_allele, $self->variant_allele;
}

sub cook_chrom {
  # STATIC
  my ($c) = @_;
  $c =~ s/^chr//i;
  return $c;
}

sub get_type {
  my ($self) = @_;
  my $type;
  if ($self->is_complex()) {
    $type = "complex";
  } elsif ($self->is_substitution()) {
    $type = $self->is_mnv ? "mnv" : "snv";
  } elsif ($self->is_insertion()) {
    $type = "ins";
  } elsif ($self->is_deletion()) {
    $type = "del";
  } else {
    dump_die($self, "can't determine type");
  }
  return $type;
}

sub import_generic {
  # native/generic import with user-specified data.
  my ($self, %options) = @_;
  my $chr = cook_chrom($options{"-reference-name"} || confess "-reference-name");
  my $pos = $options{"-base-number"} || die;
  
  foreach my $f ("-reference-allele", "-variant-allele") {
    confess(\%options, "no $f") unless exists $options{$f};
  }
  my $reference_allele = $options{"-reference-allele"} || INDEL_CHAR;
  my $variant_allele = $options{"-variant-allele"} || INDEL_CHAR;

  foreach ($reference_allele, $variant_allele) {
    if (/^\-{2,}$/) {
      # replace e.g. -- with -
      printf STDERR "Variant: ref=%s alt=%s replacing %s with -\n", $reference_allele, $variant_allele, $_;
      $_ = INDEL_CHAR;
    }
    if (/\-/ and /[ACGT]/) {
      my $stripped = $_;
      $stripped =~ s/\-+//g;
      printf STDERR "Variant: ref=%s alt=%s replacing %s with %s\n", $reference_allele, $variant_allele, $_, $stripped;
      $_ = $stripped;
    }

  }


  $self->reference_name($chr);
  $self->reference_allele($reference_allele);
  $self->variant_allele($variant_allele);

  if ($reference_allele eq INDEL_CHAR) {
    # insertion
#    dump_die(\%options, "ins", 1);
    $self->is_insertion(1);
    my $elen = length($variant_allele);
    $self->event_length($elen);
    $self->reference_allele(INDEL_CHAR);
    $self->variant_allele($variant_allele);
    $self->start($pos);
    $self->end($pos);
    # native import, no adjustment
#    dump_die($self);
  } elsif ($variant_allele eq INDEL_CHAR) {
#    dump_die(\%options, "del", 1);
    $self->is_deletion(1);
    my $elen = length($reference_allele);
    $self->event_length($elen);
    $self->reference_allele($reference_allele);
    $self->variant_allele(INDEL_CHAR);
    $self->start($pos);
    $self->end($self->start() + $elen - 1);
#    dump_die($self);
  } elsif (length($reference_allele) == length($variant_allele)) {
    # sub
    foreach ($reference_allele, $variant_allele) {
      unless (/^[ACGT\*]+$/i) {
	# * can be present in VCF data
	$self->handle_exception("ERROR: invalid nucleotide sequence in $_");
      }
    }
    $self->is_substitution(1);
    $self->start($pos);
    $self->end($pos + length($reference_allele) - 1);
    $self->reference_allele($reference_allele);
    $self->variant_allele($variant_allele);
    my $rlen = length($reference_allele);
    $self->event_length($rlen);
    $self->is_mnv($rlen > 1 ? 1 : 0);
    $self->is_snv($rlen == 1 ? 1 : 0);
    # FIX ME: share this code?
  } else {
#    dump_die(\%options, "complex", 1);
    $self->is_complex(1);
    $self->event_length(-1);
    $self->reference_allele($reference_allele);
    $self->variant_allele($variant_allele);
    $self->start($pos);
    $self->end($pos + length($reference_allele) - 1);
  }
}

sub undo_shift {
  my ($self) = @_;
  my $shift = $self->additional_shift_bases() || die;
  $self->start($self->start - $shift);
  $self->end($self->end - $shift);
}

sub import_dbnsfp_row {
  my ($self, %options) = @_;
  my $row = $options{"-row"} || die "-row";
  my $chr = $row->{"chr"} || dump_die($row, "no chr");
  my $f_pos;
  my $f_pos_v2 = "pos(1-coor)";
  if (exists $row->{$f_pos_v2}) {
    # version 2.1: we only used main position (hg19)
    $f_pos = $f_pos_v2;
  } else {
    # version 3.0+: main field is hg38, alt is hg19
    $f_pos = $options{"-f-pos"} || confess "newer dbNSFP: specify -f-pos for position";
  }

  my $pos = $row->{$f_pos} || dump_die($row, "no pos info");
  my $ref = $row->{"ref"} || dump_die($row, "no ref");
  my $alt = $row->{"alt"} || dump_die($row, "no alt");

  $self->import_generic(
			"-reference-name" => $chr,
			"-base-number" => $pos,
			"-reference-allele" => $ref,
			"-variant-allele" => $alt,
			);
}

sub handle_exception {
  my ($self, $msg) = @_;
  if ($self->exception_warn()) {
    print STDERR "$msg\n";
    $self->exception($msg);
  } else {
    confess $msg;
  }
}

sub alleles_match {
  my ($self, $v_other) = @_;
  return (($self->reference_allele() eq $v_other->reference_allele()) and
	  ($self->variant_allele() eq $v_other->variant_allele())) ? 1 : 0;
}

sub get_interbase_range {
  # get interbase coordinates for variant start/end
  my ($self) = @_;

  my ($i_start, $i_end);
  if ($self->is_substitution() or $self->is_deletion() or $self->is_complex()) {
    $i_start = $self->start() - 1;
    $i_end = $self->end();
  } elsif ($self->is_insertion()) {
    # normalized to base BEFORE the event
    $i_start = $self->start();
    $i_end = $i_start + 1;
  } else {
    die;
  }

  return ($i_start, $i_end);
}

sub import_snv4 {
  my ($self, $snv4) = @_;
#  my @f = split /\./, $snv4;
  my @f = split /\./, $snv4, -1;
  my ($chr, $pos, $ra, $va);
  if (@f == 4) {
    ($chr, $pos, $ra, $va) = @f;
  } elsif (@f == 5) {
    # period/version # in reference sequence, e.g.
    # GL000192.1.94395.-.CCTT
    (undef, undef, $pos, $ra, $va) = @f;
    $chr = join ".", @f[0, 1];
  } else {
    die "$snv4 not in SNV4 format" unless @f == 4;
  }
  die "$snv4 has invalid position $pos" unless $pos =~ /^\d+/;

  $self->import_generic(
		       "-reference-name" => $chr,
		       "-base-number" => $pos,
		       "-reference-allele" => $ra,
		       "-variant-allele" => $va,
		      );
}

sub clone {
  my ($self) = @_;
  my $v = new Variant();

  $v->import_generic(
		     "-reference-name" => $self->reference_name(),
		     "-base-number" => $self->start(),
		     "-reference-allele" => $self->reference_allele(),
		     "-variant-allele" => $self->variant_allele()
		    );
  return $v;
}

sub get_snv4 {
  my ($self) = @_;
  return join ".", $self->reference_name, $self->start, $self->reference_allele, $self->variant_allele;
}

sub is_indel {
  my ($self) = @_;
  return ($self->is_insertion or $self->is_deletion) ? 1 : 0;
}

sub matches {
  my ($self, $other) = @_;
  return ($self->reference_name eq $other->reference_name and
	  $self->start == $other->start and
	  $self->reference_allele eq $other->reference_allele and
	  $self->variant_allele eq $other->variant_allele) ? 1 : 0;
}

sub intersects {
  # TO DO: may need work for insertions/complex alleles
  # as putative matches may need to adjust for the length of the
  # net inserted sequence
  my ($self, %options) = @_;
  my $intersects;
  my $other = $options{"-variant"} || die "-variant";
  my $buffer = $options{"-buffer"} || 0;

  my $this_start = $self->start - $buffer;
  my $this_end = $self->end + $buffer;

  if ($this_end < $other->start()) {
    $intersects = 0;
  } elsif ($this_start > $other->end()) {
    $intersects = 0;
  } else {
    $intersects = 1;
  }
  return $intersects;
}

sub new_variant_from_row {
  # STATIC, exported
  #
  # recognizes the following %FLAGS entries (parsed from @ARGV):
  #
  # "-bambino",
  # "-sj-post",
  # "-f-chr=s",
  # "-f-pos=s",
  # "-f-ra=s",
  # "-f-va=s",
  # "-insertion-position=s",
  my (%options) = @_;
  my $flags = $options{"-flags"} || die;
  my $row = $options{"-row"} || die;

  my $v = new Variant();
  my $insertion_translation_needed;

  my $f_chr;
  my $f_pos;
  my $f_ra;
  my $f_va;

  if ($flags->{bambino}) {
    $f_chr = "Chr";
    $f_pos = "Pos";
    $f_ra = "Chr_Allele";
    $f_va = "Alternative_Allele";
    $insertion_translation_needed = 1;
  } elsif ($flags->{"sj-post"}) {
    $f_chr = "Chr";
    $f_pos = "WU_HG19_Pos";
    $f_ra = "ReferenceAllele";
    $f_va = "MutantAllele";
    # this format is used sometimes for SJ-native data, sometimes for external,
    # so can't say for sure whether insertion translation is needed
  } else {
    $f_chr = $flags->{"f-chr"};
    $f_pos = $flags->{"f-pos"};
    $f_ra = $flags->{"f-ra"};
    $f_va = $flags->{"f-va"};
  }
  die "-f-chr" unless $f_chr;
  die "-f-pos" unless $f_pos;
  die "-f-ra" unless $f_ra;
  die "-f-va" unless $f_va;

  $v->import_generic(
		     "-reference-name" => $row->{$f_chr},
		     "-base-number" => $row->{$f_pos},
		     "-reference-allele" => $row->{$f_ra},
		     "-variant-allele" => $row->{$f_va},
		    );

  unless (defined($insertion_translation_needed)) {
    my $ip = $flags->{"insertion-position"} || die "-insertion-position";
    if ($ip eq "before") {
      $insertion_translation_needed = 0;
    } elsif ($ip eq "after") {
      $insertion_translation_needed = 1;
    } else {
      die "-insertion-position must be before or after";
    }
  }

  if ($insertion_translation_needed and $v->is_insertion) {
    $v->start($v->start - 1);
    $v->end($v->end - 1);
  }

  return $v;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
