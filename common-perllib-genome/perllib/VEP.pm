package VEP;
# wrapper/prep for VEP (Variant Effect Predictor) 
# http://www.ensembl.org/info/docs/tools/vep/script/index.html

use strict;
use Exporter;
use FileHandle;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use TemporaryFileWrangler;
use FileUtils qw(find_binary);
use WorkingFile;

my $VEP_KEY_DELIM = "_";

@VEP::ISA = qw(Configurable Exporter);
@VEP::EXPORT_OK = qw();

use MethodMaker qw(
cache_dir
fasta
v2vep
tfw
vep_in
vep_out
vep_command
prep_only
vep_headers

buffer_size
fork_count
vep_only

prebuilt

streaming_mode
fh_results
key2v
variant_number
fh_input
wf_input
result_queue

held_result

vep_extra_params
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  $self->reset();
  return $self;
}

sub reset {
  my ($self) = @_;
  $self->tfw(new TemporaryFileWrangler());

  my $f_prebuilt = $self->prebuilt();
  my $streaming_mode = $self->streaming_mode();

  my $tfw = $self->tfw;
  my $vep_in = $tfw->get_tempfile("-append" => ".vep");
  my $vep_out;
  my $wf = 0;
  my $fh = 0;
  # 0 rather than undef so accessor methods will write peroperly
      
  if ($f_prebuilt) {
    $vep_out = $f_prebuilt;
  } else {
    $vep_out = $vep_in . ".out";
    $tfw->add_tempfile($vep_out);
    $wf = new WorkingFile($vep_in);
    $fh = $wf->output_filehandle;
  }
  $self->vep_in($vep_in);
  $self->vep_out($vep_out);
  $self->fh_input($fh);
  $self->wf_input($wf);

  my %key2v;
  $self->key2v(\%key2v);
  # non-streaming only

  $self->variant_number(0);

}

sub add_variant {
  my ($self, $v) = @_;

  my $streaming_mode = $self->streaming_mode;
  my $key = $v->get_snv4();
  my $fh = $self->fh_input();

  my $variant_number = $self->variant_number() + 1;
  $self->variant_number($variant_number);
  # for synchronization purposes in streaming mode, helpful to 
  # have the key include the variant/row number.  Otherwise
  # Bad Things might happen if e.g. input file had 2 of the same
  # variants in a row.

  #
  #  prep input file:
  #
  # http://www.ensembl.org/info/docs/tools/vep/vep_formats.html#input
  if ($streaming_mode) {
#    die;
  } else { 
    my $key2v = $self->key2v();
    die "duplicate input $key" if $key2v->{$key};
    # sanity check: shouldn't happen
    $key2v->{$key} = $v;
  }
  my $chr = $v->reference_name;
  my $ref_base = $v->reference_allele;
  my $var_base = $v->variant_allele;
  my $strand = "+";
  my $id = join $VEP_KEY_DELIM, $variant_number, $key;
  my ($start, $end);
  if ($v->is_substitution() or $v->is_deletion() or $v->is_complex()) {
    $start = $v->start;
    $end = $v->end;
  } elsif ($v->is_insertion()) {
    # "An insertion (of any size) is indicated by start coordinate = end coordinate + 1."
    $start = $v->start + 1;
    $end = $v->end;
  } else {
    die "unhandled variant type";
  }

  printf $fh "%s\t%d\t%d\t%s/%s\t%s\t%s\n",
  $chr,
  $start,
  $end,
  $ref_base,
  $var_base,
  $strand,
  $id if $fh;

}

sub run_vep {
  my ($self, %options) = @_;

  my $vep_binary;
  if (find_binary("vep")) {
    # starting with vep 88, variant_effect_predictor.pl is now called vep.
    # however research cluster has an ancient version of vep on the
    # path with variant_effect_predictor.pl which won't work with 
    # newer caches
    $vep_binary = "vep";
  } else {
    $vep_binary = "variant_effect_predictor.pl";
  }
  find_binary($vep_binary, "-die" => 1);

  if (my $wf = $self->wf_input) {
    $wf->finish();
    # close input file if results aren't prebuilt
  }
  return if $self->prep_only();

#  printf "%s\n", $vep_in; sleep 60;

  my $fasta = $self->fasta() || die "-fasta";
  my $cache_dir = $self->cache_dir() || die "-cache_dir";

  my $vep_out = $self->vep_out() || die;

  my $cmd = sprintf '%s --quiet --no_stats --refseq --hgvs --force_overwrite --offline --fasta %s -i %s -o %s --dir %s --dir_cache %s --dir_plugins %s --pubmed',
  $vep_binary,
  $fasta,
  $self->vep_in,
  $vep_out,
  $cache_dir,
  $cache_dir,
  $cache_dir;

  $cmd .= sprintf ' --buffer_size %d', $self->buffer_size if $self->buffer_size;
  $cmd .= sprintf ' --fork %d', $self->fork_count if $self->fork_count;

  $cmd .= " " . $self->vep_extra_params() if $self->vep_extra_params();

  $self->vep_command($cmd);

  unless ($self->prebuilt) {
    my $start_time = time;
    system($cmd);
    printf STDERR "VEP run time: %d\n", time - $start_time;
    die "$cmd exited with $?" if $?;
  }
  return if $self->vep_only();

  #
  #  parse output:
  #
  my %v2vep;
  
  die "no outfile $vep_out" unless -s $vep_out;

  my @headers;
  $self->vep_headers(\@headers);
  $self->v2vep({});
  my $fh_out = new FileHandle();
  $fh_out->open($vep_out) || die;
  $self->fh_results($fh_out);
  $self->get_result() unless $self->streaming_mode;
}

sub get_result {
  # in streaming mode, return next set of results.
  # in non-streaming mode, read and cache all results.
  my ($self) = @_;
  my $streaming_mode = $self->streaming_mode;
  my $fh = $self->fh_results || die;
  my $vep_headers = $self->vep_headers || die;
  my $v2vep = $self->v2vep || die;
  my $key2v = $self->key2v || die;

  my @set;
  my $last_vnum;

  if (my $r = $self->held_result()) {
    $self->held_result(0);
    return $r;
  }

  if ($streaming_mode) {
    if (my $q = $self->result_queue) {
      my $keys = $q->{Uploaded_variation} || die;
      my ($vnum, $key) = split /_/, $keys;
      $last_vnum = $vnum;

      push @set, $q;
      $self->result_queue(0);
    }
  }

  while (<$fh>) {
    chomp;
    if (/^##/) {
      next;
    } elsif (/^#/) {
      die if @{$vep_headers};
      s/^#//;
      @{$vep_headers} = split /\s+/, $_;
    } else {
#	my @f = split /\s+/, $_, -1;
      my @f = split /\t/, $_, -1;
      die sprintf "row/header mismatch: h=%d r=%d line=%s", scalar @{$vep_headers}, scalar @f, $_ unless @f == @{$vep_headers};
      my %r;
      @r{@{$vep_headers}} = @f;

      my $keys = $r{Uploaded_variation} || die;
      my ($vnum, $key) = split /_/, $keys;
      if ($streaming_mode) {
	if (not(defined($last_vnum)) or $vnum == $last_vnum) {
	  # first variant, or result for same input variant
	  push @set, \%r;
	  $last_vnum = $vnum;
	} else {
	  # entry for the next result, save for next time
	  $self->result_queue(\%r);
	  last;
	}
      } else {
	if (my $v = $key2v->{$key}) {
	  push @{$v2vep->{$v}}, \%r;
	} else {
	  dump_die(\%r, "no source Variant for $key") unless $self->prebuilt;
	  # if prebuilt, may be used with mux.pl in which case we
	  # may only be interested in a subset of results
	}
      }
    }
  }

  return @set ? \@set : undef;
}

sub get_results {
  my ($self, $v) = @_;
  # $v = same Variant.pm reference as added
  return $self->v2vep->{$v};
}

sub parse_tracking_key {
  my ($self, $key) = @_;
#  my ($vnum, $snv4) = split /_/, $key;
  # don't use split in case reference name also contains the
  # same delimiter, e.g. chr16_KI270728v1_random
  my $idx_delim = index $key, $VEP_KEY_DELIM || die "can't find key delimiter for $key";
  my $vnum = substr($key, 0, $idx_delim);
  my $snv4 = substr($key, $idx_delim + 1);

  return ($vnum, $snv4);
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
