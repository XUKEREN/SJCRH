package TabixBatchAnnotation;
# batch-annotate Variant.pm references from a tabix file
# MNE 4/2016

use strict;

use MiscUtils qw(split_list dump_die get_hash_option);
use FileUtils qw(find_binary);
use Configurable;
use Exporter;
use Counter;
use VariantMatcher;

use constant TABIX_KEY => "_tabix";
use constant TABIX_VARIANT => "_tabix_variant";

@TabixBatchAnnotation::ISA = qw(Configurable Exporter);
@TabixBatchAnnotation::EXPORT_OK = qw(TABIX_VARIANT);

use MethodMaker qw(
		    tabix
		    split_count
		    twobit
		    indel_equivalence_enable
		    indel_fuzzy_nt

		    f_tabix_chr
		    f_tabix_pos
		    f_tabix_ref_allele
		    f_tabix_var_allele

		    user_row_key
		    annotation_map

site_only
store_hits
store_site
verbose

vcf2tab
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->indel_equivalence_enable(1);
  $self->indel_fuzzy_nt(3);
  # for non-equivalence matches
#  $self->verbose(1);
  $self->configure(%options);
  return $self;
}

sub query {
  my ($self, %options) = @_;
  my $query_raw = $options{"-query"} || die "-query";
  my $site_only = $self->site_only();
  my $tag = $options{"-tag"} || "tabix batch";
  my $callback = $options{"-callback"};

  #
  #  break into chunks suitable for a single tabix batch call:
  #
  my @query_sets;
  if (my $split_count = $self->split_count) {
    @query_sets = split_list($query_raw, $split_count);
  } else {
    @query_sets = $query_raw;
  }

  my $f_tabix_chr;
  my $f_tabix_pos;
  my $f_tabix_ref_allele;
  my $f_tabix_var_allele;

  if ($self->vcf2tab()) {
    $self->f_tabix_chr("Chr");
    $self->f_tabix_pos("WU_HG19_Pos");
    $self->f_tabix_ref_allele("ReferenceAllele");
    $self->f_tabix_var_allele("MutantAllele");
    $self->tabix->vcf2tab_mode($self->vcf2tab);
  }
  
  $f_tabix_chr = $self->f_tabix_chr() || die "f_tabix_chr";
  $f_tabix_pos = $self->f_tabix_pos() || die "f_tabix_pos";
  unless ($site_only) {
    $f_tabix_ref_allele = $self->f_tabix_ref_allele() || die "f_tabix_ref_allele";
    $f_tabix_var_allele = $self->f_tabix_var_allele() || die "f_tabix_var_allele";
  }
  # TO DO:
  # - option to have user provide a callback to return variant(s)
  #   in tabix row instead

  my $field_map = $self->annotation_map();
  # mapping of tabix fields to local column names
  # TO DO: callback option
  my $store_hits = $self->store_hits();

  die "need either -annotation_map or -store_hits" unless $field_map or $store_hits;
  
  my $store_site = $self->store_site();
  # in addition to main query, also store site-only matches in specified key

  my $user_row_key = $self->user_row_key() || die "user_row_key";

  my $verbose = $self->verbose();
  my $c;
  $c = new Counter(\@query_sets) if $verbose;

  foreach my $qs (@query_sets) {
    #
    # query intervals from tabix:
    #
    my $t_rows = $self->tabix->query(
				     "-variants" => $qs,
				     "-hash" => 1
				    );
    my %by_pos;

    my $vm = new VariantMatcher(
      "-indel_equivalence_enable" => $self->indel_equivalence_enable,
      "-twobit" => $self->twobit
	);
    # always init in case tabix has no hits, e.g.
    # gnomad 16:46387006-46431456?

    if ($t_rows) {
      #
      #  parse tabix rows into Variant.pm references and load into
      #  VariantMatcher:
      #
      foreach my $r (@{$t_rows}) {
	my $v = new Variant();
	my $chr = get_hash_option($r, $f_tabix_chr);
	my $pos = get_hash_option($r, $f_tabix_pos);
	my ($ra, $va);
	if ($site_only) {
	  $ra = "A";
	  $va = "C";
	  # site-only search: bogus/placeholder alleles.
	} else {
	  $ra = get_hash_option($r, $f_tabix_ref_allele);
	  $va = get_hash_option($r, $f_tabix_var_allele);
	}
	$v->exception_warn(1);
	$v->import_generic(
			   "-reference-name" => $chr,
			   "-base-number" => $pos,
			   "-reference-allele" => $ra,
			   "-variant-allele" => $va,
			  );
	unless ($v->exception) {
	  $v->{TABIX_KEY()} = $r;
	  $vm->add_variant($v);
#	  printf STDERR "load tabix hit %s\n", $v->get_snv4();

	  my $pkey = join "_", $v->reference_name, $v->start;
	  push @{$by_pos{$pkey}}, $v;
	  # find_snv_site() will only work if primary search type is
	  # site-only, i.e. alleles are munged to SNVs.
	  # this should work for all.
	}
      }
    } elsif (0) {
      printf STDERR "no tabix hits!: dumping variants:\n";
      foreach my $v (@{$qs}) {
	printf STDERR "  %s\n", $v->get_snv4();
      }
      die "no tabix hits!";
    }

    #
    #  match query variants to tabix variants:
    #
    foreach my $v (@{$qs}) {
      my $hits;
      if ($site_only) {
	# can only query by position, e.g. ExAC coverage where
	# there is no allele information
	$hits = $vm->find_snv_site("-variant" => $v);
      } elsif ($v->is_substitution()) {
	$hits = $vm->find_snv("-variant" => $v);
      } elsif ($v->is_indel()) {
	$hits = $vm->find_indel(
				"-variant" => $v,
				"-match-basic-type" => 1,
				"-match-size" => 1,
				"-fuzz-bases" => $self->indel_fuzzy_nt,
			       );
      } elsif ($v->is_complex()) {
	$hits = $vm->find_snv("-variant" => $v);
	# if equivalence enabled, stored in SNV db
      } else {
	dump_die($v, "ERROR: unhandled query variant type", 1);
      }

      my $user_row = $v->{$user_row_key} || dump_die($v, "no $user_row_key");
      delete $user_row->{TABIX_VARIANT()};
      # field will be re-used with each batch query, so remove stale info

      foreach my $user_key (keys %{$field_map}) {
	$user_row->{$user_key} = "";
	# init to blank so column will still be present when variant
	# not found
      }

      if ($hits) {
	my $hit;
	if (@{$hits} > 1) {
	  # multiple hits in the database match (typically due to
	  # indel equivalence).  If this happens, use the match
	  # closest to the query variant.
	  # TO DO: leftalign policy?
	  my %dups = map {$_->get_snv4(), 1} @{$hits};
	  printf STDERR "WARNING: query %s matches %d database records, matches %s\n", $v->get_snv4(), scalar (@{$hits}), join ", ", sort keys %dups if $verbose;
	  my $qpos = $v->start();
	  my $dkey = "__distance__";

	  foreach (@{$hits}) {
	    $_->{$dkey} = abs($qpos - $_->start());
#	    dump_die($_, sprintf("ambiguous tabix hit for %s", $v->get_snv4), 1);
	  }
	  my @sorted = sort {$a->{$dkey} <=> $b->{$dkey}} @{$hits};
	  $hit = $sorted[0];
	} else {
	  $hit = $hits->[0];
	}
	$user_row->{TABIX_VARIANT()} = $hit;
	# Variant.pm reference to closest match

	if ($store_hits) {
	  # store tabix hits
	  my @tabix = map {$_->{TABIX_KEY()}} @{$hits};
	  $user_row->{$store_hits} = \@tabix;
#	  $user_row->{TABIX_VARIANT()} = $hits->[0];
	  # hack: store single hit for medal equivalency tag purposes
	  # might need work!
	} else {
	  # map tabix fields to user fields
	  my $tabix = $hit->{TABIX_KEY()} || dump_die($hit, "no tabix ref");
	  # we now have the tabix row corresponding to our query variant

#	printf STDERR "user:%s tabix:%s\n", $v->get_snv4, $hit->get_snv4();

	  foreach my $user_key (keys %{$field_map}) {
	    my $tabix_key = $field_map->{$user_key};
	    my $tabix_value = $tabix->{$tabix_key};
	    dump_die($tabix, "$tabix_key not defined") unless defined $tabix_value;
	    $user_row->{$user_key} = $tabix_value;
	  }
	}
      }

      if ($store_site) {
	# secondary site-only match type (e.g. ClinVar)
	my $pkey = join "_", $v->reference_name, $v->start;
	if (my $hits = $by_pos{$pkey}) {
	  my @tabix = map {$_->{TABIX_KEY()}} @{$hits};
	  $user_row->{$store_site} = \@tabix;
	}
      }
    }

    if ($callback) {
      my @co;
      @co = ("-rows" => [ map {$_->{$user_row_key} || die} @{$qs} ]);
      if (ref $callback eq "ARRAY") {
	my ($code, @args) = @{$callback};
	push @co, @args;
	&$code(@co);
      } else {
	&$callback(@co);
      }
    }

    $c->next($tag) if $c;
  }  # $qs
}

sub get_snv4 {
  # SNV4 of a tabix hit
  my ($self, $row) = @_;
  return join ".", @{$row}{
    $self->f_tabix_chr,
    $self->f_tabix_pos,
    $self->f_tabix_ref_allele,
    $self->f_tabix_var_allele,
    
  };
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
