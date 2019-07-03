package dbSNP_tabix;
# find dbSNP entries via tabix flatfile
# MNE 6/2015

use strict;

use Configurable;
use Exporter;
use TabixFile;
use MiscUtils qw(dump_die);
use GenomeUtils qw(reverse_complement);
use VariantMatcher;

use constant CLASS_SNV => "single";
use constant CLASS_MNV => "mnp";

@dbSNP_tabix::ISA = qw(Configurable Exporter);
@dbSNP_tabix::EXPORT_OK = qw();

use MethodMaker qw(
	file
substitutions_only
tf
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->substitutions_only(1);
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self) = @_;
  my $tf = new TabixFile("-file" => $self->file || die "-file");
  $tf->indel_wiggle_bases(3);
  # hack, not even implemented yet
  $self->tf($tf);
}

sub find {
  my ($self, %options) = @_;
  my $variant = $options{"-variant"} || die "-variant";

  my $results;

  my $subs_only = $self->substitutions_only();
  my $query_usable = $subs_only ? $variant->is_substitution() : 1;

  my $ra = $variant->reference_allele();
  my $va = $variant->variant_allele();
  die "MNVs not implemented, fix me" if $variant->is_mnv;
  die "indel lookup not implemented" unless $subs_only;

  if ($query_usable) {
    my $user_pos = $variant->start();

    my $rows = $self->tf->query(
				"-reference" => $variant->reference_name,
				"-pos" => $variant->start,
				"-end" => $variant->end,
				"-hash" => 1);
    if ($rows) {
#      dump_die($variant, "query", 1);

      foreach my $hit (@{$rows}) {
	my $db_ref = $hit->{refUCSC};
	my $class = $hit->{class};

	my $db_pos = ($hit->{chromStart} || die) + 1;

	my $usable = $db_pos == $user_pos;

	if ($subs_only) {
	  if ($ra) {
	    if ($variant->is_mnv) {
	      # user query type must match database entry type
	      $usable = 0 unless $class eq CLASS_MNV;
	    } else {
	      $usable = 0 unless $class eq CLASS_SNV;
	    }
	  }
	  $usable = 0 if $hit->{chromStart} == $hit->{chromEnd};
	  # rs41293467: zero span
	}

	if ($usable and $ra) {
	  if ($class eq CLASS_SNV) {
	    # for SNVs, expect the reference batch to match perfectly
	    if (length($db_ref) != length($ra)) {
	      # e.g. chr2.233243796.T.C hits rs1130338, annotated
	      # as an SNV but with multi-nucleotide reference and variant bases
	      # (feh)
	      $usable = 0;
	    } else {
	      dump_die($hit, "reference base mismatch, user=$ra db=$db_ref") if $db_ref ne uc($ra);
	    }
	  } elsif ($class eq CLASS_MNV) {
	    if (uc($ra) eq $db_ref) {
	      # good: db entry has reference sequence entry that's the
	      # same length as the MNV sequence, and it matches.
	      # (e.g. rs34714427)
	    } else {
	      #	  dump_die($hit, "MNV reference mismatch, user=$ra db=$db_ref");
	      $usable = 0;
	      # some of these might salvagable if the desired allele
	      # appears in the allele set even if the given reference
	      # base is only 1 nt.
	    }
	  } else {
	    dump_die($hit, "unhandled class $class, user=$ra db=$db_ref");
	  }
	}

	if ($usable and $va) {
	  # user specified variant allele
	  my $vb = $self->get_variant_bases($hit);

	  if (grep {$va eq $_} @{$vb}) {
	    # ok
	  } else {
	    #	dump_die($options{"-sj"}, "SJ query:", 1);
	    #	dump_die($hit, "dbSNP mismatch!");
	    # different variant at site
	    $usable = 0;
	  }
	}

	#    dump_die($hit, "debug $usable", 1);

	if ($usable) {
	  $results = [] unless $results;
	  push @{$results}, $hit;
	}
      }
    }
  }

  return $results;
}

sub get_variant_bases {
  # return genomic variant bases for a SNV
  # (correcting for strand)
  my ($self, $hit) = @_;
  my $observed = $hit->{observed} || die;
  my $strand = $hit->{strand} || die;
  my @entries = split /\//, $observed;
  my $reference_base = uc($hit->{refUCSC} || die);
  my $saw_reference;
  my @variant_bases;
  foreach my $entry (@entries) {
    if ($strand eq "-") {
      $entry = reverse_complement($entry);
    } elsif ($strand ne "+") {
      die;
    }

    $entry = uc($entry);
    if ($entry eq $reference_base) {
      $saw_reference = 1;
    } else {
      push @variant_bases, $entry;
    }
  }

#  rs144571919

#  dump_die($hit, sprintf("WARNING: didn't find reference %s in %s (%s)!",
#			 $reference_base, $observed, $hit->{name}), 1)
#      unless $saw_reference;
  # turn off: happens frequently for ?complex? indels

  # rs144571919:
  # reference is a C, observed bases are A/C, strand is +...WTF?
  return \@variant_bases;
}

sub get_vm_for_variants {
  #
  # query intervals for a set of variants, parsing results into
  # a VariantMatcher instance.
  #
  my ($self, %options) = @_;
  my $variants = $options{"-variants"} || die "-variants";
  my $subs_only = $self->substitutions_only();
  die "indels not supported" unless $subs_only;

#  printf STDERR "intervals: %s\n", join " ", @intervals;

  my $hits = $self->tf->query(
			      "-variants" => $variants,
			      "-hash" => 1,
			     );
  my $vm = new VariantMatcher();

  if ($hits) {
    foreach my $r (@{$hits}) {
      if ($subs_only) {
	my $class = $r->{class};
	next unless $class eq CLASS_SNV or $class eq CLASS_MNV;
      }
      my $ra = uc($r->{refUCSC} || die);
      next if $ra eq "-";
      # zero span, etc.

      my $pos = $r->{chromStart} + 1;
      # convert from interbase to in-base
      my $vbs = $self->get_variant_bases($r);
      foreach my $va (@{$vbs}) {
	next unless length($ra) == length($va);
	# some mnps have ref length of 1 and variant length > 1

	my $v = new Variant();
#      printf STDERR "add variant %s\n", join ".", $r{chrom}, $pos, $ra, $va;
	$v->exception_warn(1);

	$v->import_generic(
	  "-reference-name" => ($r->{chrom} || die),
	  "-base-number" => $pos,
	  "-reference-allele" => $ra,
	  "-variant-allele" => $va
	    );
	if ($v->exception) {
#	dump_die($v, "failed variant parsing", 1);
	  dump_die($r, "failed variant parsing", 1);
	} else {
	  $v->{dbsnp} = $r;
	  $vm->add_snv(
	    "-row" => $v,
	    "-variant" => 1,
	      );
	}
      }
    }
  }

  return $vm;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
