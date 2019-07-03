package DNANexusWrapper;
# helper for dnanexus runs, parsing app's JSON
# MNE 2/2017

use strict;
use Exporter;
use FileHandle;
use POSIX qw(ceil);

use JSON;

use Configurable;
use MiscUtils qw(get_hash_option dump_die log_message shell_cmd);
use FileUtils qw(write_simple_file read_simple_file count_file_lines);

@DNANexusWrapper::ISA = qw(Configurable Exporter);
@DNANexusWrapper::EXPORT_OK = qw();

use MethodMaker qw(
	cache_job
        cache_outputs

json
input_data
output_policy
clean_policy
app

wait_for_job
clean_io
clean_job

auto_instance
job_chunk_size
min_cores
max_cores
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->cache_job(1);
  $self->wait_for_job(1);
  $self->cache_outputs(1);
  $self->input_data({});
  $self->output_policy({});
  $self->clean_policy({});
  $self->clean_io(1);
  $self->clean_job(1);
  $self->configure(%options);
  return $self;
}

sub json_setup {
  my ($self, %options) = @_;
  $self->json($self->get_json(%options));
}

sub get_json {
  my ($self, %options) = @_;
  my $fh = $options{"-fh"} || die "-fh";
  local $/ = undef;
  my $json_text = <$fh>;
  my $json = JSON->new->allow_nonref;
  return $json->decode($json_text);
}

sub set_input {
  my ($self, %options) = @_;
  my $parameter = get_hash_option(\%options, "-parameter");
  my $value = get_hash_option(\%options, "-value");
  my $data = $self->input_data();
  die "duplicate param $parameter" if exists $data->{$parameter};
  $data->{$parameter} = $value;
  if (defined $options{"-clean"}) {
    $self->clean_policy()->{$parameter} = $options{"-clean"};
  }

}

sub set_output_policy {
  my ($self, %options) = @_;
  my $parameter = get_hash_option(\%options, "-parameter");
  my $value = get_hash_option(\%options, "-policy");
  my $data = $self->output_policy();
  die "duplicate param $parameter" if exists $data->{$parameter};
  $data->{$parameter} = $value;
}

sub run {
  my ($self) = @_;
  # TO DO:
  # - instance type
  my $json = $self->json();
  my $input_data = $self->input_data();
  my $output_policy = $self->output_policy();
  my $cache_job = $self->cache_job();

  my $instance_type;
  if (my $job_file = $self->auto_instance()) {
    # try to choose an instance appropriate for the amount of work
    my $line_count = count_file_lines($job_file);
    if (0) {
      printf STDERR "fix me\n";
      $line_count = 300000;
    }

    my $job_chunk_size = $self->job_chunk_size() || die;
    my $job_count = ceil($line_count / $job_chunk_size);

    my $min_cores = $self->min_cores() || 0;
    my $max_cores = $self->max_cores();
    my @all_core_choices = (2, 4, 8, 16, 32);
    # dx run --instance-type-help

    my @core_choices = grep {$_ >= $min_cores and $_ <= $max_cores} @all_core_choices;

    my $minimum_busy_ratio = 0.80;
    # require that this fraction of cores be utilized when moving
    # to a larger instance.  Don't want to jump to a much higher
    # instance (e.g. from 8 to 16 cores) if the additional ones
    # will be mostly idle.

    my $cores_wanted = $core_choices[0];

    printf STDERR "auto-instance: work:%d\n", $line_count;
    foreach my $c (@core_choices) {
      last if $max_cores and $c > $max_cores;

      my $jobs_handled = $c * $job_chunk_size;
      my $ratio = $line_count / $jobs_handled;

      printf STDERR "  cores:%d handled:%d ratio:%f\n", $c, $jobs_handled, $ratio;
      last if $ratio < $minimum_busy_ratio;
      # don't move to a higher instance count: not busy enough

      $cores_wanted = $c;
      # this instance can be kept busy enough

      last if $job_count <= $c;
    }
    printf STDERR "  requesting %d cores\n", $cores_wanted;

    $instance_type = sprintf 'mem1_ssd1_x%d', $cores_wanted;
  }


  # identify tracking file for job ID cache; default to first input file.
  my $inputSpec = $json->{inputSpec} || die;
  my $p_tracking;
  foreach my $input (@{$inputSpec}) {
    my $p = $input->{name};
    if ($input->{class} eq "file") {
      $p_tracking = $p;
      last;
    }
  }

  # verify handling of all output files is specified:
  my $outputSpec = $json->{outputSpec} || die;
  my $f_job_cache;

  foreach my $output (@{$outputSpec}) {
    if ($output->{class} eq "file") {
      my $p = $output->{name};
      die sprintf "output %s has no policy specified\n", $p unless $output_policy->{$p};
      unless ($f_job_cache) {
	# 
	my $outfile = $self->get_outfile("-base" => $p_tracking,
					 "-name" => $p);
	$f_job_cache = $outfile . ".dx.job";
      }
    } else {
      die "non-file output, ???";
    }
  }

  my $job_id;
  if ($cache_job) {
    die unless $f_job_cache;
    if (-s $f_job_cache) {
      my $lines = read_simple_file($f_job_cache);
      $job_id = $lines->[0];
    }
  }

  my @p_input;
  my @p_output;

  my %param2id;
  my @p;
  # verify all inputs are present and upload local files if necessary:
  foreach my $input (@{$inputSpec}) {
    my $p = $input->{name};
    # currently all parameters treated as manadatory;
    # needs improvement for optional params
    my $v = $input_data->{$p};
    die "input for parameter $p not specified" unless defined $v;
    push @p, $p;

    if ($input->{class} eq "file") {
      my $clean_ok = $self->clean_policy->{$p};
      die "no clean policy specified for input file $p" unless defined $clean_ok;

      push @p_input, $p if $clean_ok;

      unless ($job_id) {
	printf STDERR "uploading %s...\n", $v;
	$param2id{$p} = run_dx_id(sprintf 'dx upload "%s" --brief', $v);
      }
    }
  }

  if ($job_id) {
    log_message("job $job_id already submitted");
  } else {
    my $run_name = sprintf '%s%s', ($self->app ? "app-" : ""), $json->{name};
    my $cmd = sprintf 'dx run %s', $run_name;
    foreach my $p (@p) {
      my $v = $param2id{$p} || $input_data->{$p};
      $cmd .= sprintf ' -i%s=%s', $p, $v;
    }
    $cmd .= " -y --brief";
    $cmd .= sprintf " --instance-type %s", $instance_type if $instance_type;

    $job_id = run_dx_id($cmd);
    write_simple_file([ $job_id ], $f_job_cache);
  }

  # TO DO:
  # maybe skip this check if outfiles exist??
  if ($self->wait_for_job()) {
    my $cmd = sprintf 'dx wait %s', $job_id;
    shell_cmd("-cmd" => $cmd);
  }

  #
  # job is finished:
  #

  foreach my $output (@{$outputSpec}) {
    if ($output->{class} eq "file") {
      my $p = $output->{name};
      my $f_out_local = $self->get_outfile("-base" => $p_tracking,
					   "-name" => $p);
      my $need_download = not(-f $f_out_local);
      push @p_output, $p;

      if ($need_download) {
	my $cmd = sprintf "dx download %s:%s -o %s", $job_id, $p, $f_out_local;
	shell_cmd("-cmd" => $cmd);
      }
    } else {
      die "non-file output, ???";
    }
  }

  if ($self->clean_io()) {
    # remove cloud inputs and outputs
    my $json;
    my $fh;
    $fh = new FileHandle();
    if (0) {
      print STDERR "DEBUG: hardcoded json\n";
      $fh->open("debug_json.txt");
    } else {
      $fh->open(sprintf 'dx describe %s --json|', $job_id);
    }
    my $json = $self->get_json("-fh" => $fh);

    my @rm;
    foreach my $p (@p_input) {
      # delete input files only if explicitly OK,
      # (in case of eventual shared/immutable inputs)
      # may be too paranoid
      if ($self->clean_policy->{$p}) {
	push @rm, $json->{input}{$p}{'$dnanexus_link'} || die;
      }
    }

    foreach my $p (@p_output) {
      # assume output files are all safe to delete once downloaded
      push @rm, $json->{output}{$p}{'$dnanexus_link'} || die;
    }

    foreach my $id (@rm) {
      shell_cmd("-cmd" => "dx rm $id");
    }
  }

  if ($self->clean_job() and $cache_job) {
    unlink $f_job_cache;
  }

}

sub run_dx_id {
  # STATIC
  my ($cmd) = @_;
  my $id = shell_cmd("-cmd" => $cmd, "-backtick" => 1) || die;
  return $id;
}

sub get_outfile {
  my ($self, %options) = @_;
  my $p_base = $options{"-base"} || die "-base";
  my $p = $options{"-name"} || die "-name";
  my $infile = $self->input_data->{$p_base} || die;
  my $policy = $self->output_policy->{$p} || die "policy for $p";  
  my $outfile;
  if (my $suffix = $policy->{suffix}) {
    $suffix = "." . $suffix unless $suffix =~ /^\./;
    $outfile = $infile . $suffix;
  } else {
    die "suffix is only implemented policy";
  }
  return $outfile;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
