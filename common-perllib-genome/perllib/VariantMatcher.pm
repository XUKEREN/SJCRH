package VariantMatcher;
# standardized variant comparisons for SNVs and indels:
# - protein AA changes (specific and codon-number based)
# - SNVs
# - etc.
#
# MNE 10/2013

use strict;
use Configurable;

use Exporter;
use Carp qw(confess cluck);

use List::Util qw(min max);

use AAParser;
use GenomeUtils qw(cook_chromosome_name);
use BucketMap;
use ReferenceSanityCheck;
use MiscUtils qw(dump_die);
use Variant qw(INDEL_CHAR);
use IndelEquivalence;
use GeneSymbolMapper qw(new_gsm_lite);

use constant BUCKET_SIZE => 100000;

use constant MATCH_TYPE_NONE => 0;
use constant MATCH_TYPE_PERFECT => 1;
use constant MATCH_TYPE_CODON_ONLY => 2;

use constant MAX_DELETION_SIZE => 100;
# maximum possible size of a deletion (better too high than too low)

@VariantMatcher::ISA = qw(Configurable Exporter);
@VariantMatcher::EXPORT_OK = qw(VM_VKEY);

my $VERBOSE = 0;

use MethodMaker qw(
		    db_aa_specific
		    db_aa_codon_only
		    db_aa_literal
                    db_snv
                    db_snv_pos
db_snv_basenum
db_variant_literal

db_deletion
db_insertion

		    aap
rsc
fasta_dir
match_type

enable_literal_aa_lookup
enable_literal_variant_lookup

mnv_warning

indel_equivalence_enable
twobit
db_ie
indel_equivalence_max_distance

gsm
gsm_enable
		  );

use constant SNV_DELIM => ".";

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->aap(new AAParser());
  $self->db_aa_specific({});
  # exact match for a specified protein AA change
  $self->db_aa_codon_only({});
  # very fuzzy lookup based on codon number ONLY; catchall bucket.
  # TO DO:
  # additional search types for more specific entries?
  # insertions/deletions, in-frame/frameshift, etc.?
  $self->db_aa_literal({});
  $self->db_variant_literal({});
  $self->db_snv({});
  $self->db_snv_pos({});
  $self->db_snv_basenum({});
  $self->db_deletion({});
  $self->db_insertion({});
  $self->db_ie({});
  $self->indel_equivalence_max_distance(100);
  $self->gsm_enable(1);
  $self->configure(%options);

  return $self;
}

sub add_aa_substitution {
  # attempt to parse an entry as a substitution
  my ($self, %options) = @_;
  my $gene = $options{"-gene"} || die "-gene";
  my $aa = $options{"-aa"} || die "-aa";
  # gene symbol to bucket along with variant.
  # could use refseq accession, etc.
  my $row = $options{"-row"} || die "-row";
  my $aap = $self->aap();

  my $parsable = 0;
  if (my $cooked = $aap->parse_substitution($aa)) {
    # formatted as simple AA substitution
    my $cnum = $aap->codon_number() || die;
    stamp_codons($row, $cnum, $cnum, $aa);
    push @{$self->db_aa_specific()->{$gene}{$cooked}}, $row;
    $parsable = 1;
  }
  return $parsable;
}

sub add_aa_codons {
  # attempt to parse codons ONLY of an AA event
  my ($self, %options) = @_;
  my $gene = $options{"-gene"} || die "-gene";
  my $aa = $options{"-aa"} || die "-aa";
  # gene symbol to bucket along with variant.
  # could use refseq accession, etc.
  my $row = $options{"-row"} || die "-row";
  my $aap = $self->aap();

  my $parsable = 0;
  if ($aap->parse($aa)) {
    # insertions, deletions, frameshifts
    # much more variety in parsing
    $parsable = 1;
    my $start = $aap->codon_start();
    my $end = $aap->codon_end();
    stamp_codons($row, $start, $end, $aa);
    for (my $i = $start; $i <= $end; $i++) {
      push @{$self->db_aa_codon_only()->{$gene}{$i}}, $row;
    }
  }
  return $parsable;
}

sub add_aa_literal {
  # literal/trusted AA annotation
  my ($self, %options) = @_;
  my $gene = $options{"-gene"} || die "-gene";
  my $aa = $options{"-aa"} || die "-aa";
  my $row = $options{"-row"} || die "-row";
  push @{$self->db_aa_literal()->{$gene}{$aa}}, $row;
  return 1;
}

sub add_aa {
  # add a reference variant AA sequence.
  # for SNVs, strictly formatted protein change.
  # for indels, looser formatting to account for real-world messiness.
  my ($self, %options) = @_;
  my $aa = $options{"-aa"} || confess "specify -aa";
  my $parsable = 0;
  $parsable = 1 if $self->add_aa_codons(%options);
  $parsable = 1 if $self->add_aa_substitution(%options);
  # for substitutions, add to both databases in case user wants
  # a fuzzy search (i.e. codon-number only alternative match for
  # a precisely-described substitution)
  $parsable = 1 if $self->enable_literal_aa_lookup() and
      $self->add_aa_literal(%options);

  if ($VERBOSE) {
    cluck sprintf "add_aa: can't parse AA entry: %s\n", $aa unless $parsable or $aa eq "p.?";
  }
  return $parsable;
}

sub find_aa_substitution {
  my ($self, %options) = @_;

  if (my $r = $options{"-sj"}) {
    # SJ format
    $options{"-gene"} = $r->{GeneName} || die;
    $options{"-aa"} = $r->{AAChange} || "";
  }

  my $gene = $self->get_lookup_gene($options{"-gene"} || die "-gene");
  my $aa = $options{"-aa"};
  die "-aa" unless defined $aa;
  my $result;
  if ($aa) {
    # might not be present, e.g. promoter
    my $aap = $self->aap();
    if (my $cooked = $aap->parse_substitution($aa)) {
      $result = $self->db_aa_specific()->{$gene}{$cooked};
    }
  }
  return $result;
}

sub get_lookup_gene {
  # translate query gene symbol into symbol used in this database
  my ($self, $gene) = @_;
  my $gene_db;
  if ($self->gsm_enable()) {
    $self->gsm_setup();
    $gene_db = $self->gsm->find($gene) || $gene;
    # if not present, just use original symbol rather than null,
    # to avoid uninitialized variable errors
  } else {
    $gene_db = $gene;
  }
  return $gene_db;
}



sub find_aa_literal {
  my ($self, %options) = @_;

  if (my $r = $options{"-sj"}) {
    # SJ format
    $options{"-gene"} = $r->{GeneName} || die;
    $options{"-aa"} = $r->{AAChange} || die;
  }

  my $gene = $self->get_lookup_gene($options{"-gene"} || die "-gene");
  my $aa = $options{"-aa"} || die "-aa";

  return $self->db_aa_literal()->{$gene}{$aa};
}

sub find_aa_single {
  my ($self, %options) = @_;
  my $gene = $self->get_lookup_gene($options{"-gene"} || die "-gene");
  my $aa = $options{"-aa"} || die "-aa";
  my $aap = $self->aap();

  my $result;
  if ($aap->parse($aa)) {
    my $start = $aap->codon_start();
    my $end = $aap->codon_end();
    if ($start == $end) {
      # single codon only
      my $hits = $self->db_aa_codon_only()->{$gene}{$start};
      if ($hits) {
	$result = [] unless $result;
	die "FIX ME: need to verify DB ENTRY HITS ONLY THIS CODON";
	push @{$result}, @{$hits} if $hits;
      }
    }
  } else {
    printf STDERR "can't parse codon info in %s\n", $aa if $VERBOSE;
  }
  return $result;
}

sub find_aa_specific_codon {
  # match to exact start/end in query event only
  # typically a substitution but could be a deletion, etc.
  my ($self, %options) = @_;
  my $gene = $self->get_lookup_gene($options{"-gene"} || die "-gene");
  my $aa = $options{"-aa"} || die "-aa";
  my $aap = $self->aap();

  my $result;
  if ($aap->parse($aa)) {
    my $start = $aap->codon_start();
    my $end = $aap->codon_end();
    for (my $i = $start; $i <= $end; $i++) {
      my $hits = $self->db_aa_codon_only()->{$gene}{$i};
      if ($hits) {
	foreach my $hit (@{$hits}) {
	  my $hit_start = $hit->{SJ_codon_start} || die;
	  my $hit_end = $hit->{SJ_codon_end} || die;
	  if ($start == $hit_start and $end == $hit_end) {
	    $result = [] unless $result;
	    push @{$result}, $hit;
	  }
	}
      }
    }
  } else {
    printf STDERR "can't parse codon info in %s\n", $aa if $VERBOSE;
  }
  return $result;
  
}

sub find_aa_codons {
  #
  # broad match: any event touching this AA's codon annotation
  #
  my ($self, %options) = @_;
  my $gene = $self->get_lookup_gene($options{"-gene"} || die "-gene");
  my $aa = $options{"-aa"} || die "-aa";
  my $aap = $self->aap();

  my $result;
  if ($aap->parse($aa)) {
    my $start = $aap->codon_start();
    my $end = $aap->codon_end();
    for (my $i = $start; $i <= $end; $i++) {
      my $hits = $self->db_aa_codon_only()->{$gene}{$i};
      if ($hits) {
	$result = [] unless $result;
	push @{$result}, @{$hits} if $hits;
      }
      # FIX ME: duplicates
    }
  } else {
    printf STDERR "can't parse codon info in %s\n", $aa if $VERBOSE;
  }
  return $result;
}

sub find_aa {
  # high level AA lookup:
  # first try specific substitution, then generic codon match
  my ($self, %options) = @_;
  my $gene = $options{"-gene"} || die "-gene";
  my $aa = $options{"-aa"} || die "-aa";
  my $result;
  my $match_type = MATCH_TYPE_NONE;
  if ($result = $self->find_aa_substitution(%options)) {
    $match_type = MATCH_TYPE_PERFECT;
  } elsif ($result = $self->find_aa_codons(%options)) {
    $match_type = MATCH_TYPE_CODON_ONLY;
  }
  $self->match_type($match_type);
  return $result;
}

sub add_snv {
  # misnomer:
  # add SNV or dinucleotide, etc. (as long as same # of bases)
  my ($self, %options) = @_;
  my $row = $options{"-row"} || die "-row";
  if ($options{"-gedi"}) {
    # row is a GeDI database record
    $options{"-reference"} = $row->{chromosome} || die;
    $options{"-base-number"} = $row->{pos} || die;
    $options{"-reference-base"} = $row->{reference_allele} || die;
    $options{"-variant-base"} = $row->{non_reference_allele} || die;
  } elsif ($options{"-sj"}) {
    # row is a SJ postprocessed variant report-style record
    $options{"-reference"} = $row->{Chr} || dump_die($row, "no Chr");
    $options{"-base-number"} = get_sj_pos($row);
    $options{"-reference-base"} = $row->{ReferenceAllele} || die;
    $options{"-variant-base"} = $row->{MutantAllele} || die;
  } elsif ($options{"-bambino"}) {
    # raw Bambino report
    $options{"-reference"} = $row->{Chr} || die;
    $options{"-base-number"} = $row->{Pos} || die;
    $options{"-reference-base"} = $row->{Chr_Allele} || die;
    $options{"-variant-base"} = $row->{Alternative_Allele} || die;
  } elsif ($options{"-variant"}) {
    # Variant.pm
    $options{"-reference"} = $row->reference_name || die;
    $options{"-base-number"} = $row->start || die;
    $options{"-reference-base"} = $row->reference_allele || die;
    $options{"-variant-base"} = $row->variant_allele || die;
  }

  my $ref_raw = $options{"-reference"} || die "-reference";
  my $ref_cooked = cook_chromosome_name($ref_raw, "-return-unknown" => 1);
  # non-canonical names, e.g. chr1_KI270706v1_random
  my $base_num = $options{"-base-number"} || die "-base-number";
  stamp_pos($row, $base_num, $base_num);

  my $ref_base = uc($options{"-reference-base"} || die "-reference-base");
  my $var_base = uc($options{"-variant-base"} || die "-variant-base");
  confess "WTF: ref=$ref_base var=$var_base" unless length($ref_base) == length($var_base);
  # SNV, dinucleotide, etc.

  foreach ($ref_base, $var_base) {
    confess "WTF: possible attempt to add indel to SNV list: ref=$ref_base var=$var_base" if /\-/;
  }
  confess "WTF: possible attempt to add indel to SNV list: ref=$ref_base var=$var_base" if length($ref_base) != length($var_base);

  if ($self->fasta_dir()) {
    my $rsc = $self->get_rsc();
    my $rsc_ok = $rsc->check("-ref-name" => $ref_cooked,
			     "-base-number" => $base_num,
#			"-ref-base" => $ref_base)) {
			     "-ref-base" => substr($ref_base,0,1),
			     # just check 1st base
			     # (might be di/tri etc.)
			    );
    unless ($rsc_ok) {
      cluck sprintf "ERROR: reference sanity check failed for %s.%s.%s->%s, NOT adding! (actual=%s)", $ref_cooked, $base_num, $ref_base, $var_base,
      $rsc->get_reference_base("-ref-name" => $ref_cooked, "-base-number" => $base_num);
      return 0;
    }
  }

  $self->db_snv_pos->{$ref_cooked}{$base_num}{$ref_base} = 1;
  # for sanity check in case query needs - complement, etc.

  my $key = join SNV_DELIM, $ref_cooked, $base_num, $ref_base, $var_base;
#  printf STDERR "add key %s\n", $key;
  push @{$self->db_snv()->{$key}}, $row;

#  dump_die($row, "save $ref_cooked $base_num", 1);
  push @{$self->db_snv_basenum()->{$ref_cooked}{$base_num}}, $row;

  return 1;
}

sub find_snv {
  my ($self, %options) = @_;
  if (my $r = $options{"-sj"}) {
    # SJ postprocessed variant report
    if (exists $r->{Chr}) {
      $options{"-reference"} = $r->{Chr};
      $options{"-base-number"} = get_sj_pos($r);
      $options{"-reference-base"} = $r->{ReferenceAllele};
      $options{"-variant-base"} = $r->{MutantAllele};
    } else {
      return undef;
      # certain test inputs (e.g. SJ gold db) don't have 
      # nucleotide-level annotations and so only AA lookups can be run.
    }
  } elsif (my $v = $options{"-variant"}) {
    $options{"-reference"} = $v->reference_name();
    $options{"-base-number"} = $v->start();
    $options{"-reference-base"} = $v->reference_allele();
    $options{"-variant-base"} = $v->variant_allele();
  }

  my $ref_raw = $options{"-reference"} || confess "-reference";
  my $ref_cooked = cook_chromosome_name($ref_raw, "-return-unknown" => 1);
  # in case chrom in query doesn't match db set
  my $base_num = $options{"-base-number"} || die "-base-number";
  my $ref_base = uc($options{"-reference-base"} || confess "-reference-base");
  my $var_base = uc($options{"-variant-base"} || confess "-variant-base");

  my $info = $self->db_snv_pos->{$ref_cooked}{$base_num};
  my $mnv_warning = 0;
  if ($info and $ref_base ne "-") {
    # variant set has a record at this position, make sure reference
    # base matches query base (unless an insertion)
    unless ($info->{$ref_base}) {
      my $lookup_len = length($ref_base);
      my %all_len;
      foreach (keys %{$info}) {
	my $l = length($_);
	$all_len{$l} = 1;
      }

      if ($all_len{$lookup_len}) {
	# reference base mismatch
	cluck sprintf "SNV lookup sanity check failed: query for %s.%s has reference %s, we have reference %s",
	$ref_cooked, $base_num, $ref_base, join ",", keys %{$info};
	# just warn rather than die, might be some known problems
	# e.g. IARC TP53 17.7579472 reference mismatch
      } else {
	# possible SNV vs. MNV mismatch
	# if e.g. matching Bambino SNV(s) vs. GATK MNV
	printf STDERR "SNV lookup warning: possible attempt to match vs. MNV, query %s.%s has reference %s, db has reference %s\n",
	$ref_cooked, $base_num, $ref_base, join ",", keys %{$info};

	$mnv_warning = max(keys %all_len);
	if (scalar keys %all_len > 1) {
	  printf STDERR "warning: possible multiple MNVs at %s.%d, lengths=%s\n",
	  $ref_cooked, $base_num, join ",", keys %all_len;
	}
	($mnv_warning) = keys %all_len;
      }
    }
  }
  $self->mnv_warning($mnv_warning);

  my $key = join SNV_DELIM, $ref_cooked, $base_num, $ref_base, $var_base;
  return $self->db_snv()->{$key};
}

sub add_insertion {
  my ($self, %options) = @_;

  if ($self->indel_equivalence_enable) {
    $self->ie_add_variant(%options);
    return;
  }

  my $row = $options{"-row"} || die "-row";

  if ($options{"-gedi"}) {
    my $ra = $row->{"reference_allele"} || die;
    unless ($ra eq "-") {
      printf STDERR "ERROR: can't handle ?complex insertion with ref allele %s\n", $ra;
      return;
    }
    my $inserted = $row->{"non_reference_allele"} || die;

    $options{"-reference"} = $row->{chromosome} || die;
    $options{"-start"} = $row->{pos} || die;
    # FIX ME: are these Bambino coordinates? adjusted??
    $options{"-end"} = $row->{pos} || die;
    # I forget what this was for; complex "delins"?
    $options{"-count"} = length $inserted;
  } elsif ($options{"-sj"}) {
    my $ra = $row->{ReferenceAllele} || die;
    die unless $ra eq "-";
    my $inserted = $row->{MutantAllele} || die;
    $options{"-reference"} = $row->{Chr} || die;
    $options{"-start"} = get_sj_pos($row);
    $options{"-end"} = get_sj_pos($row);
    $options{"-count"} = length $inserted;
  } elsif ($options{"-bambino"}) {
    my $ra = $row->{ReferenceAllele};
    die unless $ra eq "";
    my $inserted = $row->{Alternative_Allele} || die;
    $options{"-reference"} = $row->{Chr} || die;
    my $pos = $row->{Pos} || die;
    # FIX ME: adjust by -1?
    $options{"-start"} = $pos;
    $options{"-end"} = $pos;
    $options{"-count"} = length $inserted;
  } elsif ($options{"-variant"}) {
    my $ra = $row->reference_allele();
    die unless $ra eq INDEL_CHAR;
    my $inserted = $row->variant_allele() || die;
    $options{"-reference"} = $row->reference_name() || die;
    my $pos = $row->start() || die;
    # FIX ME: adjust by -1?
    $options{"-start"} = $pos;
    $options{"-end"} = $pos;
    $options{"-count"} = length $inserted;
  }

  my $ref_cooked = cook_chromosome_name($options{"-reference"} || die "-reference");
  my $start = $options{"-start"} || die "-start";
  my $end = $options{"-end"} || die "-end";
  my $count = $options{"-count"} || die "-count";

  stamp_pos($row, $start, $end);

  my $db = $self->db_insertion();
  my $bm = $db->{$ref_cooked};
  $bm = $db->{$ref_cooked} = new BucketMap("-chunk" => BUCKET_SIZE) unless $bm;

  my %record;
  $record{start} = $start;
  $record{end} = $end;
  $record{row} = $row;
  $record{size} = $count;
  # event size

  $bm->add_range(
    "-start" => $start,
    "-end" => $end,
    "-value" => \%record
      );

  $self->add_literal_variant(%options) if $self->enable_literal_variant_lookup();
}

sub add_literal_variant {
  my ($self, %options) = @_;
  my $row = $options{"-row"} || die "-row";

  if ($options{"-variant"}) {
    # Variant.pm
    $options{"-reference"} = $row->reference_name || die;
    $options{"-start"} = $row->start || die;
    $options{"-reference-base"} = $row->reference_allele || die;
    $options{"-variant-base"} = $row->variant_allele || die;
  }

  my $ref_cooked = cook_chromosome_name($options{"-reference"} || die "-reference");
  my $start = $options{"-start"} || die "-start";
  my $allele_ref = $row->{ReferenceAllele} || $options{"-reference-base"} || confess "-reference-base";
  my $allele_var = $row->{MutantAllele} || $options{"-variant-base"} || confess "-variant-base";
  my $key = join SNV_DELIM, $ref_cooked, $start, $allele_ref, $allele_var;
#  printf STDERR "adding key %s\n", $key;
  push @{$self->db_variant_literal->{$key}}, $row;
}


sub add_deletion {
  my ($self, %options) = @_;

  if ($self->indel_equivalence_enable) {
    $self->ie_add_variant(%options);
    return;
  }

  my $row = $options{"-row"} || die "-row";
  if ($options{"-gedi"}) {
    my $ra = $row->{"reference_allele"} || die;
    my $va = $row->{"non_reference_allele"} || die;
    unless ($va eq "-") {
      printf STDERR "ERROR: can't handle ?complex insertion with ref allele %s\n", $va;
      return;
    }
    $options{"-reference"} = $row->{chromosome} || die;
    $options{"-start"} = $row->{pos} || die;
    $options{"-end"} = $row->{pos} + length($ra) - 1;
  } elsif ($options{"-sj"}) {
    # row is a SJ postprocessed variant report-style record
    my $ra = $row->{"ReferenceAllele"} || die;
    my $va = $row->{"MutantAllele"} || die;
    die $va unless $va eq "-";
    $options{"-reference"} = $row->{Chr} || die;
    my $pos = get_sj_pos($row) || die;
    $options{"-start"} = $pos;
    $options{"-end"} = $pos + length($ra) - 1;
#    printf STDERR "start:%d end:%d\n", $options{"-start"}, $options{"-end"};
  } elsif ($options{"-bambino"}) {
    # raw Bambino call
    my $ra = $row->{"Chr_Allele"} || die;
    my $va = $row->{"Alternative_Allele"};
    die $va unless $va eq "";
    $options{"-reference"} = $row->{Chr} || die;
    my $pos = $row->{Pos} || die;
    $options{"-start"} = $pos;
    $options{"-end"} = $pos + length($ra) - 1;
#    dump_die($row, "Debug", 1);
#    dump_die(\%options);
  } elsif ($options{"-variant"}) {
    my $ra = $row->reference_allele || die;
    my $va = $row->variant_allele || die;
    die $va unless $va eq "-";
    $options{"-reference"} = $row->reference_name() || die;
    my $pos = $row->start || die;
    $options{"-start"} = $pos;
    $options{"-end"} = $pos + length($ra) - 1;
  }

  my $ref_cooked = cook_chromosome_name($options{"-reference"} || die "-reference");
  my $start = $options{"-start"} || die "-start";
  my $end = $options{"-end"} || confess "-end";
  stamp_pos($row, $start, $end);

  my $length = ($end - $start) + 1;

  my $db = $self->db_deletion();
  my $bm = $db->{$ref_cooked};
  $bm = $db->{$ref_cooked} = new BucketMap("-chunk" => BUCKET_SIZE) unless $bm;

  my %record;
  $record{start} = $start;
  $record{end} = $end;
  $record{row} = $row;
  $record{size} = $length;

  $bm->add_range(
    "-start" => $start,
    "-end" => $end,
    "-value" => \%record
      );

  $self->add_literal_variant(%options) if $self->enable_literal_variant_lookup();

}

sub find_indel {
  #
  # find matching indels
  #
  my ($self, %options) = @_;

  return $self->ie_find(%options) if $self->indel_equivalence_enable;

  my $fuzz = $options{"-fuzz-bases"} || 0;
  my $match_basic_type = $options{"-match-basic-type"};
  confess "specify -match-basic-type" unless defined $match_basic_type;
  # if true, reference matches must be of same basic type (insertion/deletion)
  # counterargument: both insertions and deletions can cause a frameshift
  # at the same site
  my $match_size = $options{"-match-size"};
  confess "specify -match-size" unless defined $match_size;
  # if set, requires that matching indel have the same # of bases
  # inserted or deleted

  my ($is_insertion, $is_deletion);

  my $query_event_size;

  if (my $r = $options{"-sj"}) {
    # query is from SJ postprocessed variant report
#    dump_die($r, "SJ post debug", 1);
    $options{"-reference"} = $r->{Chr};

#    foreach (sort keys %{$r}) {
#      printf STDERR "%s: %s\n", $_, $r->{$_};
#    }

    my $start = get_sj_pos($r);
    my $ref_a = $r->{ReferenceAllele} || die;
    my $var_a = $r->{MutantAllele} || die;

    my $end;

    if ($ref_a eq "-") {
      # insertion
      # locations may need some tweaking:
      # - in Bambino, reported base # is first mapped base AFTER inserted bases
      # - postprocessed reports may shift the locations; not sure whether
      #   this also affects the base #
      # - so, lookups may require some tolerance, probably a good idea
      #   to use at least 1 nt of fuzz
      $is_insertion = 1;
      $end = $start;
      $query_event_size = length($var_a);
    } elsif ($var_a eq "-") {
      # deletion
      $is_deletion = 1;
      $end = $start + length($ref_a) - 1;
      $query_event_size = length($ref_a);
    } else {
      # this variant is not an indel!
      # don't search.
      die "ERROR: unparsed indel $ref_a $var_a" if $ref_a =~ /\-/ or $var_a =~ /\-/;
      return undef;
    }

    dump_die($r, "no end") unless $end;

    $options{"-start"} = $start;
    $options{"-end"} = $end;
  } elsif (my $v = $options{"-variant"}) {
    # Variant.pm
    if ($v->is_insertion()) {
      $is_insertion = 1;
    } elsif ($v->is_deletion()) {
      $is_deletion = 1;
    } else {
      die "unhandled variant type";
    }
    $options{"-reference"} = $v->reference_name();
    $query_event_size = $v->event_length() || die;
    $options{"-start"} = $v->start();
    $options{"-end"} = $v->end();
  }

  die "fix me: don't know indel type" unless $is_insertion or $is_deletion;
  die "no query event size" unless $query_event_size;

  my $ref_cooked = cook_chromosome_name($options{"-reference"} || die "-reference");
  my $start = $options{"-start"} || die "-start";
  my $end = $options{"-end"} || confess "-end";

#  printf STDERR "query: %d-%d\n", $start, $end;

  my $fuzzy_start = $start - $fuzz;
  my $fuzzy_end = $end + $fuzz;

  my @hits;
  my %saw;

  my @search_db;
  if ($match_basic_type) {
    if ($is_insertion) {
      @search_db = $self->db_insertion();
    } elsif ($is_deletion) {
      @search_db = $self->db_deletion();
    } else {
      die;
    }
  } else {
    @search_db = ($self->db_insertion(), $self->db_deletion());
  }

  foreach my $db (@search_db) {
    my $bm = $db->{$ref_cooked} || next;
    my $hits = $bm->fuzzy_find("-start" => $fuzzy_start, "-end" => $fuzzy_end);
    foreach my $hit (@{$hits}) {
      next if $hit->{end} < $fuzzy_start;
      next if $hit->{start} > $fuzzy_end;

      if ($match_size) {
	unless ($hit->{size}) {
	  dump_die($hit, "database variant doesn't have size, ???", 1);
	  dump_die($hit->{row}, "database row dump", 1);
	}
	next unless $query_event_size == $hit->{size};
      }

      if (0) {
	printf STDERR "WARNING: DISTANCE FIX DISABLED!\n";
      } else {
	my $distance = abs($start - $hit->{start});
	next if $distance > $fuzz;
      }

      next if $saw{$hit};
      $saw{$hit} = 1;
      push @hits, $hit;
    }
  }

  if ($options{"-any-overlap"}) {
    my $qs = $start - MAX_DELETION_SIZE;
    my $qe = $end + MAX_DELETION_SIZE;
#    print STDERR "query $start $end $qs $qe\n";

    foreach my $db (@search_db) {
      my $bm = $db->{$ref_cooked} || next;
      my $hits = $bm->fuzzy_find("-start" => $qs, "-end" => $qe);
      foreach my $hit (@{$hits}) {
	next if $hit->{end} < $start;
	next if $hit->{start} > $end;
	next if $saw{$hit};
	$saw{$hit} = 1;
	push @hits, $hit;
      }
    }
  }

  if (@hits) {
    return [ map {$_->{row}} @hits ];
  } else {
    return undef;
  }
}

sub get_rsc {
  my ($self) = @_;
  my $rsc = $self->rsc();
  unless ($rsc) {
    # init
    $rsc = new ReferenceSanityCheck("-fasta_dir" => 
				    $self->fasta_dir || die "need -fasta-dir");
    $self->rsc($rsc);
  }
  return $rsc;
}

sub get_aa_genes {
  # find all genes (AA additions only)
  my ($self) = @_;
  my %genes;
  foreach my $db (
    $self->db_aa_specific(),
    $self->db_aa_codon_only(),
      ) {
    foreach (keys %{$db}) {
      $genes{$_} = 1;
    }
  }
  return \%genes;
}

sub find_snv_site {
  my ($self, %options) = @_;
  if (my $r = $options{"-sj"}) {
    # SJ postprocessed variant report
    if (exists $r->{Chr}) {
      $options{"-reference"} = $r->{Chr};
      $options{"-base-number"} = get_sj_pos($r);
    } else {
      return undef;
      # certain test inputs (e.g. SJ gold db) don't have 
      # nucleotide-level annotations and so only AA lookups can be run.
    }
  } elsif (my $v = $options{"-variant"}) {
    $options{"-reference"} = $v->reference_name();
    $options{"-base-number"} = $v->start();
  }

  my $ref_cooked = cook_chromosome_name(($options{"-reference"} || confess "-reference"), "-return-unknown" => 1);
  my $base_num = $options{"-base-number"} || die "-base-number";

#  printf STDERR "check %s %s\n", $ref_cooked, $base_num;

  my $hits = $self->db_snv_basenum()->{$ref_cooked}{$base_num};
  return $hits;
}

sub stamp_codons {
  my ($row, $start, $end, $aa) = @_;
#  printf STDERR "DEBUG codon stamp for %s: %d-%d\n", $aa, $start, $end;
  my $contention;
  $row->{SJ_AA_history}{$aa} = 1;
  # overkill?

  if (my $old_start = $row->{SJ_codon_start}) {
    my @v = ($start, $old_start);
    my %v = map {$_, 1} @v;
    if (scalar keys %v > 1) {
      # legit for dinucleotides, e.g. TP53 [p.K132N;M133L]
      my $msg = sprintf "conflicting AA in same record for %s start, v=%s\n", $aa, join ",", @v;
      confess "epic fail $msg" if scalar keys %{$row->{SJ_AA_history}} == 1;
      # possible parsing inconsistency between substitition and general parser?
    }
    $row->{SJ_codon_start} = min(@v);
  } else {
    $row->{SJ_codon_start} = $start;
  }

  if (my $old_end = $row->{SJ_codon_end}) {
    my @v = ($end, $old_end);
    my %v = map {$_, 1} @v;
    if (scalar keys %v > 1) {
      my $msg = sprintf "conflicting AA in same record for %s end, v=%s\n", $aa, join ",", @v ;
      die "epic fail $msg" if scalar keys %{$row->{SJ_AA_history}} == 1;
    }
    $row->{SJ_codon_end} = max(@v);
  } else {
    $row->{SJ_codon_end} = $end;
  }
}

sub stamp_pos {
  my ($row, $start, $end) = @_;
  if (my $old_start = $row->{SJ_genomic_start}) {
    my @v = ($start, $old_start);
    $row->{SJ_genomic_start} = min(@v);
  } else {
    $row->{SJ_genomic_start} = $start;
  }

  if (my $old_end = $row->{SJ_genomic_end}) {
    my @v = ($end, $old_end);
    $row->{SJ_genomic_end} = max(@v);
  } else {
    $row->{SJ_genomic_end} = $end;
  }
}

sub find_literal_variant {
  my ($self, %options) = @_;
  my $key;
  if (my $r = $options{"-sj"}) {
    $key = join SNV_DELIM,
    (cook_chromosome_name($r->{Chr} || die)),
    (get_sj_pos($r) || die),
    ($r->{ReferenceAllele} || die),
    ($r->{MutantAllele} || die);
  } elsif (my $r2 = $options{"-variant"}) {
    $key = join SNV_DELIM,
    ($r2->reference_name || die),
    ($r2->start || die),
    ($r2->reference_allele || die),
    ($r2->variant_allele || die);
  }

  return $self->db_variant_literal->{$key};
}

sub get_sj_pos {
  my ($row) = @_;
  die unless $row;
  my $pos = $row->{WU_HG19_Pos} || $row->{WU_HG18_Pos};
  dump_die($row, "can't find base number") unless $pos;
  return $pos;
}

sub add_gedi_row {
  my ($self, %options) = @_;
  my $row = $options{"-row"} || die;

  my $ra = $row->{reference_allele};
  unless ($ra) {
    printf STDERR "ERROR: blank reference allele\n";
    return;
  }

  my $va = $row->{non_reference_allele} || die;
  die "say what?" if $ra =~ /\-/ and $va =~ /\-/;

  #
  #  add nucleotide-based entry:
  #
  if ($ra =~ /\-/) {
    $self->add_insertion(
			 "-row" => $row,
			 "-gedi" => 1
			);
  } elsif ($va =~ /\-/) {
    $self->add_deletion(
			"-row" => $row,
			"-gedi" => 1
		       );
  } else {
    $self->add_snv(
		   "-row" => $row,
		   "-gedi" => 1,
		  );
  }
}

sub add_complex {
  #
  # complex variants are only partially/clunkily supported.
  #
  my ($self, %options) = @_;
  my $row = $options{"-row"} || die;
  die unless $options{"-variant"};
  # Variant.pm only for now

  $self->add_literal_variant("-row" => $row, "-variant" => 1);

#  dump_die($row, "Debug", 1);
  
  if ($self->indel_equivalence_enable()) {
    # before equivalence, kludge below created 2 subvariants,
    # one deletion and one insertion.  This interferes with
    # indel equivalence mode where input data has custom user-provided
    # hash keys in Variant.pm references.
    # workaround: store in SNV lookup table.
    my $key = join SNV_DELIM, $row->reference_name, $row->start, $row->reference_allele, $row->variant_allele;
    push @{$self->db_snv()->{$key}}, $row;
    # hack
#    die "hey now";
  } else {
    my $v_del = new Variant();
    $v_del->import_generic(
      "-reference-name" => $row->reference_name,
      "-base-number" => $row->start,
      "-reference-allele" => $row->reference_allele,
      "-variant-allele" => INDEL_CHAR
	);
    $self->add_deletion("-row" => $v_del, "-variant" => 1);

    my $v_ins = new Variant();
    $v_ins->import_generic(
      "-reference-name" => $row->reference_name,
      "-base-number" => $row->start,
      "-reference-allele" => INDEL_CHAR,
      "-variant-allele" => $row->variant_allele,
	);
    $self->add_insertion("-row" => $v_ins, "-variant" => 1);
  }
}

sub find_indel_window {
  #
  # find any indels within the specified start/end.
  # used for matches to complex variants.
  # not sure if fuzz bases should apply as the complex variant
  # should have resolved the full range of the variant.
  #
  my ($self, %options) = @_;

  if (my $v = $options{"-variant"}) {
    # Variant.pm
    $options{"-reference"} = $v->reference_name();
    $options{"-start"} = $v->start();
    $options{"-end"} = $v->end();
  }

  my $ref_cooked = cook_chromosome_name($options{"-reference"} || die "-reference");
  my $start = $options{"-start"} || die "-start";
  my $end = $options{"-end"} || die "-end";

  my $fuzzy_start = $start;
  my $fuzzy_end = $start;

  my @hits;
  my %saw;

  my @search_db = ($self->db_insertion(), $self->db_deletion());

  foreach my $db (@search_db) {
    my $bm = $db->{$ref_cooked} || next;
    my $hits = $bm->fuzzy_find("-start" => $fuzzy_start, "-end" => $fuzzy_end);
    foreach my $hit (@{$hits}) {
      next if $hit->{end} < $fuzzy_start;
      next if $hit->{start} > $fuzzy_end;
      next if $saw{$hit};
      $saw{$hit} = 1;
      push @hits, $hit;
    }
  }

  if (@hits) {
    return [ map {$_->{row}} @hits ];
  } else {
    return undef;
  }
}

sub add_variant {
  my ($self, $v, %options) = @_;
  if (my $row = $options{"-row"}) {
    # clone additional data into variant reference.
    # awkward: we typically want to track the raw db row data,
    # but this add method uses a Variant.pm reference.
    foreach my $k (keys %{$row}) {
      if (exists $v->{$k}) {
#	printf STDERR "skipping %s\n", $k;
      } else {
#	printf STDERR "clone %s\n", $k;
	$v->{$k} = $row->{$k};
      }
    }
  }

  if ($v->is_complex) {
    $self->add_complex(
		       "-row" => $v,
		       "-variant" => 1,
		      );
  } elsif ($v->is_substitution) {
    $self->add_snv(
		   "-row" => $v,
		   "-variant" => 1,
		  );
  } elsif ($v->is_deletion) {
    $self->add_deletion(
			"-row" => $v,
			"-variant" => 1,
		       );
  } elsif ($v->is_insertion) {
    $self->add_insertion(
			"-row" => $v,
			"-variant" => 1,
		       );
  } else {
    dump_die($v, "unhandled variant type");
  }
}

sub ie_add_variant {
  #
  #  add an indel to equivalence database, bucket by chrom/type/size
  #
  my ($self, %options) = @_;
  my $v = $options{"-row"} || die;
  if ($v->is_insertion() or $v->is_deletion()) {
    my $key = $self->get_ie_hash_key($v);
    my $db = $self->db_ie() || die;
#    printf STDERR "IE db add %s under %s\n", $v->get_snv4(), $key;
    push @{$db->{$key}}, $v;
    # bucket indels by chromosome, type, and size to minimize comparisons
  } else {
    dump_die($v, "not an indel");
  }
}

sub get_ie_hash_key {
  my ($self, $v) = @_;
  return join "_", $v->reference_name, $v->get_type, $v->event_length();
}

sub ie_find {
  #
  # find an equivalent indel
  #
  my ($self, %options) = @_;
  my $v_query = $options{"-variant"} || die "-variant";
  my $key = $self->get_ie_hash_key($v_query);
  my $results;
  if (my $raw_set = $self->db_ie->{$key}) {
    my $ie = new IndelEquivalence("-twobit" => $self->twobit);

    if (0) {
      dump_die($v_query, "query", 1);
      foreach (@{$raw_set}) {
	dump_die($_, "db", 1);
      }
    }

    my @query;
    my $max_distance = $self->indel_equivalence_max_distance();
    foreach my $v (@{$raw_set}) {
      my $distance = abs($v_query->start - $v->start);

      push @query, $v if $distance <= $max_distance;
    }

#    printf STDERR "IE: raw:%d filtered:%d\n", scalar @{$raw_set}, scalar @query;

    if (@query) {
      $results = $ie->compare_one(
				  "-query" => $v_query,
				  "-database" => \@query,
				 );
    }

  }
  return $results;
}

sub gsm_setup {
  my ($self) = @_;
  unless ($self->gsm()) {
    # initialize gene symbol disambiguation; delay until first query
    my $gsm = new_gsm_lite();
    $self->gsm($gsm);
    foreach my $gene (keys %{$self->db_aa_specific()}) {
      $gsm->add_gene("-gene" => $gene);
    }
  }
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
