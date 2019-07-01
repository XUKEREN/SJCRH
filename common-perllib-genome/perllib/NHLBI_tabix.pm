package NHLBI_tabix;
# tabix search of NHLBI records

use strict;
use Configurable;
use Exporter;

use NHLBIParser;
use TabixFile;
use GenomeUtils qw(cook_chromosome_name);
use MiscUtils qw(dump_die);

@NHLBI_tabix::ISA = qw(Configurable Exporter);
@NHLBI_tabix::EXPORT_OK = qw();

use MethodMaker qw(
file
tf
np
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self) = @_;
  my $file = $self->file || die;
  my $tf = new TabixFile("-file" => $file);
  $tf->indel_wiggle_bases(0);
  # indels not implemented for now
  $self->tf($tf);
  $self->np(new NHLBIParser());
}

sub find_one {
  my ($self, %options) = @_;
  my $v;
  if (my $r = $options{"-sj-post"}) {
    $v = new Variant();
    $v->import_bambino_row("-row" => $r,
			   "-postprocessed" => 1);
  }
  die unless $v;

  my $base_reference = $v->reference_allele;
  my $base_variant = $v->variant_allele;

  foreach ($base_reference, $base_variant) {
    my $bad;
    $bad = 1 if $_ and $_ eq "-";
    # query is an indel: won't match
    $bad = 1 if length($_) > 1;
    # MNV: won't match
    return [] if $bad;
  }

  my $tf = $self->tf();

  my $hits_tabix = $tf->query(
			    "-variants" => [ $v ],
			    "-hash" => 1
			   );
  my $np = $self->np;

  my $hits_raw = [];
  if (@{$hits_tabix}) {
    my @hits;
    foreach my $h (@{$hits_tabix}) {
      my $parsed = $np->parse_line($h);
      push @{$hits_raw}, $parsed;
    }
  } else {
    return [];
  }

  my @results;

  foreach my $hit (@{$hits_raw}) {
    my $usable = 1;

    if ($base_reference) {
      # reference base given: ensure consistent
      #	dump_die($hit, "debug vs user $base_reference", 1);
      if ($hit->{reference_base} ne $base_reference) {
	dump_die($hit, "reference sanity fail vs. user $base_reference");
	die "reference sanity fail!!";
	$usable = 0;
      }
    }

    if ($base_variant) {
      # variant base given: report only records containing this base
      my $found;
      foreach my $vb (@{$hit->{variant_bases}}) {
	$found = 1 if uc($vb) eq uc($base_variant);
      }
      #	dump_die($hit);
      #	die join ",", @{$hit->{variant_bases}};
      $usable = 0 unless $found;
    }

    push @results, $hit if $usable;
    #      dump_die($hit, "hit debug", 1) if $usable;
  }

  return \@results;
}

sub get_parser {
  my ($self) = @_;
  return $self->np;
}

sub get_vm_for_variants {
  #
  # query intervals for a set of variants, parsing results into
  # a VariantMatcher instance.
  #
  my ($self, %options) = @_;
  my $variants = $options{"-variants"} || die "-variants";

  my $hits = $self->tf->query(
			      "-variants" => $variants,
			      "-hash" => 1,
			     );

  my $vm = new VariantMatcher();
  my $np = $self->np() || die;

  if ($hits) {
    foreach my $r (@{$hits}) {
      my $parsed = $np->parse_line($r);
      my $chr = $parsed->{Chromosome} || die;
      my $pos = $parsed->{Position} || die;
      my $ref_base = $parsed->{reference_base} || die;
      foreach my $var_base (@{$parsed->{variant_bases}}) {
#      dump_die($parsed);
	my $v = new Variant();
	$v->import_generic(
	  "-reference-name" => $chr,
	  "-base-number" => $pos,
	  "-reference-allele" => $ref_base,
	  "-variant-allele" => $var_base
	    );
	if ($v->exception) {
#	dump_die($v, "failed variant parsing", 1);
	  dump_die($r, "failed variant parsing", 1);
	} else {
	  $v->{nhlbi} = $r;
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
