package IndelEquivalence;
# wrapper to Stephen Rice's "indelmatch"
# http://hc-wiki.stjude.org/display/compbio/How+to+find+equivalent+indels
# MNE 12/2015
#
# for insertions, all input Variant.pm references are expected to be
# standardized to the base number BEFORE the indel, as used by most
# databases.  This intermediate code will translate to the base
# number AFTER the insertion, as used by indelmatch.  This process
# should be transparent to the user.
#
# TO DO:
# - ability to stream in possible candidate variants, filtering
#   to only those matching the same type/size as the query variant?

use strict;
use Exporter;

use MiscUtils qw(get_hash_option dump_die);
use Configurable;
use FileUtils qw(find_binary write_simple_file);
use FileHandle;
use TemporaryFileWrangler;

use constant TAG_QUERY => "query";
use constant TAG_DATABASE => "db";

use constant INDELMATCH_INSERTIONS_BASE_AFTER => 1;
my $WHERE_DELIMITER = ";";
# e.g. if the query and database variant match exactly,
# both locations will be reported in one column

#use constant MUNGE_MT => 1;

@IndelEquivalence::ISA = qw(Configurable Exporter);
@IndelEquivalence::EXPORT_OK = qw();

use MethodMaker qw(
	input_lines
snv4db
twobit
tfw
outfile
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  $self->tfw(new TemporaryFileWrangler());
  return $self;
}

sub reset {
  my ($self) = @_;
  $self->input_lines([]);
  $self->snv4db({});
  unlink $self->outfile if $self->outfile;
}

#sub query {
sub compare_one {
  #
  # see if a single query variant matches one or more database variants
  #
  my ($self, %options) = @_;
  find_binary("indelmatch", "-die" => 1);
  my $v_query = get_hash_option(\%options, "-query");
  # the query variant: a single Variant.pm reference 
  my $database =  get_hash_option(\%options, "-database");
  # an arrayref of Variant.pm references to search
  die "no db variants" unless @{$database};

  $self->reset();

  my @results;
  if (@{$database} == 1 and $v_query->matches($database->[0])) {
    # if database contains only one variant and it's formatted the same way,
    # no need to run equivalence code
    push @results, $database->[0];
  } else {
    $self->add_input_line($v_query, 0);
    foreach my $v (@{$database}) {
      $self->add_input_line($v, 1);
    }

    my $fh = $self->run_equiv();

    my %hits;
    while (<$fh>) {
      chomp;
      my @f = split /\t/, $_;
      die unless @f >= 3;
      splice(@f, 2, 1);
      # remove result count field
      my $saw_query;
      my %equiv;
      while (@f) {
	my ($variant, $where_raw) = splice(@f, 0, 2);
	my @where = split /$WHERE_DELIMITER/, $where_raw;

	foreach my $where (@where) {
	  if ($where eq TAG_QUERY) {
	    $saw_query = 1;
	  } elsif ($where eq TAG_DATABASE) {
	    $equiv{$variant} = 1;
	  } else {
	    die "unhandled $where";
	  }
	}
      }

      if ($saw_query) {
	die "WTF" if %hits;
	# should only happen once
	%hits = %equiv;
      }
    }

    my $db = $self->snv4db || die;
    foreach my $hit (keys %hits) {
      my $key = $self->get_db_key("-snv4" => $hit);
      my $set = $db->{$key} || die "can't find db ref for $key";
      push @results, @{$set};
    }
  }

  return @results ? \@results : undef;
}

sub add_input_line {
  my ($self, $v, $is_db) = @_;

  if ($v->is_indel) {
    my $ref = $v->reference_name;
    my $pos = $v->start;
    $pos++ if $v->is_insertion and INDELMATCH_INSERTIONS_BASE_AFTER;
    my $ra = $v->reference_allele;
    my $va = $v->variant_allele;
    my $tag = $is_db ? TAG_DATABASE : TAG_QUERY;
    my $snv4 = join ".", $ref, $pos, $ra, $va;

    push @{$self->snv4db()->{$snv4}}, $v if $is_db;

    push @{$self->input_lines}, join "\t", $snv4, $tag;
  } else {
    dump_die($v, "variant is not an indel");
  }
}

sub run_equiv {
  my ($self) = @_;
  my $tfw = $self->tfw;
  my $outfile = $tfw->get_tempfile("-append" => ".indelequiv");
  $self->outfile($outfile);
#  printf STDERR "outfile: $outfile\n";

  my $twobit = $self->twobit || die "-twobit";

  my $cmd = sprintf '|indelmatch %s > %s', $self->twobit, $outfile;
  my $fh = new FileHandle();
  $fh->open($cmd) || die;
  foreach my $l (@{$self->input_lines}) {
    printf $fh "%s\n", $l;
  }
  $fh->close() || die;
  die "$cmd exited with $?" if $?;
  # TO DO: capture STDERR to look for warnings about e.g. reference mismatches

  if (0) {
    # debug
    write_simple_file($self->input_lines, "ie.in");
    printf STDERR "cat ie.in | indelmatch %s\n", $self->twobit;
    die;
  }


  $fh->open($outfile);
  return $fh;
}

sub find_equivalences {
  #
  # find equivalences in a set of variants
  #
  my ($self, %options) = @_;
  find_binary("indelmatch", "-die" => 1);
  my $database =  get_hash_option(\%options, "-database");
  # an arrayref of Variant.pm references to search
  die "no db variants" unless @{$database};

  $self->reset();
  foreach my $v (@{$database}) {
    $self->add_input_line($v, 1);
  }

  my $fh = $self->run_equiv();
  my $db = $self->snv4db || die;
  my @equivs;
  while (<$fh>) {
    chomp;
    my @f = split /\t/, $_;
    die unless @f >= 3;
    splice(@f, 2, 1);
    # remove result count field
    my @equiv;
    push @equivs, \@equiv;
    while (@f) {
      my ($variant, $where_raw) = splice(@f, 0, 2);
      my $key = $self->get_db_key("-snv4" => $variant);
      my $set = $db->{$key} || die "can't find db ref for $key";
      push @equiv, @{$set};
    }
  }
  return \@equivs;
}

sub get_db_key {
  my ($self, %options) = @_;
  my ($chr, $pos, $ra, $va);
  if (my $snv4 = $options{"-snv4"}) {
    ($chr, $pos, $ra, $va) = split /\./, $snv4;
  } else {
    die;
  }
  $chr =~ s/^chr//;
  # cache contains cooked Variant.pm reference names,
  # indel equivalence code always uses "chr" prefix even if input doesn't

#  $chr = "MT" if MUNGE_MT and $chr eq "M";

  return join ".", $chr, $pos, $ra, $va;
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
