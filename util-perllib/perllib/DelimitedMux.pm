package DelimitedMux;
# split and rejoin delimited files (helpful for parallelization)
# MNE 12/2014

use strict;
use Carp qw(confess);

use Configurable;
use Exporter;
use File::Basename;
use File::Path;
use Cwd;
use FileHandle;

use DelimitedFile;
use Reporter;
use MiscUtils qw(dump_die shell_cmd);
use Cluster;
use FileUtils qw(write_simple_file read_simple_file newer_than);

@DelimitedMux::ISA = qw(Configurable Exporter);
@DelimitedMux::EXPORT_OK = qw();

use MethodMaker qw(
	file
        mux_tag
split_file_cache
	split_files
	subdir_mode

	job_out_files
cache
cleanup_type
time_start
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->mux_tag("mux");
  $self->cache(1);
  $self->configure(%options);
  return $self;
}

sub split_file {
  # demux
  my ($self, %options) = @_;
  my $file = $options{"-file"} || $self->file() || die;
  my $rows_per_file  = $options{"-lines"} || die "-lines";
  # TO DO: alternative method breaking into user-specified number of
  # fragments
  my $cache = $self->cache();
  my $cf = sprintf "%s.mux_files_%d", basename($file), $rows_per_file;
  $self->split_file_cache($cf);
  my $needed = 1;
  my @files;
  if ($cache and -s $cf and not(newer_than($file, $cf))) {
    # usable cache
    my $files = read_simple_file($cf);
    $needed = 0;
    foreach (@{$files}) {
      $needed = 1 unless -f $_;
    }
    @files = @{$files} unless $needed;
  }
  if ($needed) {
    my $fh_in = new FileHandle();
    $fh_in->open($file) || die;
    my $header_line = <$fh_in>;
    my $processed = 0;
    my $file_counter = 0;
    my $subdir_mode = $self->subdir_mode;
    my $wf;
    my $fh_out;
    while (my $line = <$fh_in>) {
      if ($processed++ % $rows_per_file == 0) {
	my $file_num = $file_counter++;
	$wf->finish() if $wf;

	my $fmt = join ".", basename($file), $self->mux_tag(), $rows_per_file, $file_num;

	my $fn;
	if ($subdir_mode) {
	  $fn = join "/", $fmt, $fmt;
	  my $dir = dirname($fn);
	  unless (-d $dir) {
	    mkpath($dir) || die;
	  }
	  die unless -d $dir;
	} else {
	  $fn = $fmt;
	}

	push @files, $fn;

	$wf = new WorkingFile($fn);
	$fh_out = $wf->output_filehandle();
	print $fh_out $header_line;
      }
      print $fh_out $line;
    }
    $wf->finish() if $wf;
    write_simple_file(\@files, $cf) if $cache;
  }

  return $self->split_files(\@files);
}

sub run_jobs {
  my ($self, %options) = @_;

  $self->time_start(time);
  if ($options{"-parallel"}) {
    # use GNU parallel instead of submitting jobs to cluster
    return $self->run_jobs_parallel(%options);
  } elsif ($options{"-pool"}) {
    # use a pool of jobs: punt to different implementation
    return $self->run_jobs_pooled(%options);
  }

  my $serial_mode = $options{"-serial"};
  my $bsub_cores = $options{"-bsub-cores"};

  my $template = $options{"-template"} || die "-template";
  my $ram = $options{"-ram"} || die "-ram";
  die "ram must be an integer" unless $ram =~ /^\d+$/;
  my $out_suffix = $options{"-out-suffix"} || die "-out-suffix";
  $out_suffix = "." . $out_suffix unless $out_suffix =~ /^\./;

  my @outfiles;

  my $start_dir = getcwd();
  my $subdir_mode = $self->subdir_mode();

  foreach my $fn (@{$self->split_files()}) {
    my ($cmd, $outfile);

    if ($subdir_mode) {
      my $sd = dirname($fn);
      chdir($sd) || die "can't cd to $sd";
      my $bn = basename($fn);
      $cmd = sprintf $template, $bn;
      $outfile = $bn . $out_suffix;
      # use basename for Cluster.pm since we are now in a subdir
      push @outfiles, $fn . $out_suffix;
      # use fully qualified name for final output list
    } else {
      $cmd = sprintf $template, $fn;
      $outfile = $fn . $out_suffix;
      push @outfiles, $outfile;
    }
    die "failed to place file in template $template" if $template eq $cmd;

    if ($serial_mode) {
      shell_cmd("-cmd" => $cmd);
    } else {
      my $c = new Cluster(
	"-outfile" => $outfile,
	"-project" => "PCGP",
	  );

#      $c->debug(1);

      $c->processors($bsub_cores) if $bsub_cores;
      $c->node_class("");
      $c->memory_reserve_mb($ram);
      $c->memory_limit_mb($ram);
      $c->command($cmd);
      $c->run();
    }

   (chdir($start_dir) || die) if $subdir_mode;
  }

  if (my $sleep_time = $options{"-wait"}) {
    $sleep_time = 60 if $sleep_time == 1;
    my $total = scalar @outfiles;
    while (1) {
      my $have = 0;
      foreach my $fn (@outfiles) {
	$have++ if -f $fn;
      }

      if ($have == $total) {
	# done
	last;
      } else {
	printf STDERR "have %d/%d (%.1f%%); elapsed:%d, waiting for %d sec...\n",
	  $have, $total, ($have * 100 / $total),
	    (time - $self->time_start), $sleep_time;
	sleep $sleep_time;
      }
    }
  }

  return $self->job_out_files(\@outfiles);
}

sub join_files {
  my ($self, %options) = @_;
  my $stdout = $options{"-stdout"};
  my $outfile = $options{"-out"};
  die unless $outfile or $stdout;
  my $compress = $options{"-compress"};
  $outfile = "-" if $stdout;
  my $job_files = $self->job_out_files() || die;

  my $rpt;
  my @clean_files;
  my $cleanup_type = $self->cleanup_type();

  my $have_all = 1;
  my $count_missing = 0;
  foreach my $jf (@{$job_files}) {
    unless (-f $jf) {
      $have_all = 0;
      $count_missing++;
    }
  }

  my $result;

  if ($have_all) {
    printf STDERR "joining results to %s...\n", $outfile;

    my @o;
    push @o, ("-compress" => $compress) if $compress;
    my $wf;
    my $fh_out;
    if ($stdout) {
      $fh_out = *main::STDOUT;
    } else {
      $wf = new WorkingFile($outfile, @o);
      $fh_out = $wf->output_filehandle();
    }
    my $wrote_header;

    foreach my $jf (@{$job_files}) {
      printf STDERR "  %s...\n", $jf;

      my $fh_in = new FileHandle();
      $fh_in->open($jf);

      my $header = <$fh_in>;
      if ($wrote_header) {
	die "header mismatch" unless $header eq $wrote_header;
      } else {
	print $fh_out $header;
	$wrote_header = $header;
      }

      while (my $line = <$fh_in>) {
	print $fh_out $line;
      }

      if ($cleanup_type and $cleanup_type eq "lite") {
	# clean cluster log files only, leave split files
	# so we can rerun if desired
	push @clean_files, $jf;
	foreach my $suffix (qw(out err bsub job)) {
	  push @clean_files, sprintf "%s.cluster.%s", $jf, $suffix;
	}
      }
    }
    $wf->finish();

    if ($cleanup_type) {
      # remove all intermediate files based on glob of split file.
      if ($cleanup_type eq "glob") {
	push @clean_files, $self->split_file_cache();
	foreach my $sf (@{$self->split_files}) {
	  push @clean_files, $sf;
	  push @clean_files, glob($sf . ".*");
	}
      } elsif ($cleanup_type eq "lite") {
	# already handled
      } else {
	die "unimplemented cleanup mode: $cleanup_type";
      }
    }

    unlink(@clean_files) if @clean_files;
    $result = 1;
  } else {
    printf STDERR "can't join: %d files not available\n", $count_missing;
    $result = 0;
  }

  return $result;

}

sub run_jobs_pooled {
  my ($self, %options) = @_;
  my $slot_count = $options{"-pool"} || die "-pool";
  # maximum number of jobs to run 
  my $template = $options{"-template"} || die "-template";
  my $ram = $options{"-ram"} || die "-ram";
  die "ram must be an integer" unless $ram =~ /^\d+$/;
  my $sleep_time = $options{"-wait"};
  $sleep_time = 60 if ($sleep_time || 1) == 1;
  my $bsub_cores = $options{"-bsub-cores"};

  my $out_suffix = $options{"-out-suffix"} || die "-out-suffix";
  $out_suffix = "." . $out_suffix unless $out_suffix =~ /^\./;

  my @outfiles;

  die "subdir more not implemented for -pool" if $self->subdir_mode();

  my %out2cmd;
  foreach my $fn (@{$self->split_files()}) {
    my ($cmd, $outfile);

    $cmd = sprintf $template, $fn;
    die "failed to place file in template $template" if $template eq $cmd;
    $outfile = $fn . $out_suffix;
    push @outfiles, $outfile;

    die "duplicate" if $out2cmd{$outfile};
    $out2cmd{$outfile} = $cmd unless -e $outfile;
  }

  my %needed = %out2cmd;
  my @slots;

  my $total = scalar @outfiles;

  while (1) {
    my $start_time = time;

    my $have = 0;
    foreach my $fn (@outfiles) {
      $have++ if -e $fn;
    }

    if ($have == $total) {
      # done
      last;
    } else {
      if (%needed) {
	# jobs still need to be submitted
	for (my $i = 0; $i < $slot_count; $i++) {
	  if (not($slots[$i]) or -f $slots[$i]) {
	    # empty slot or job for this slot is finished
	    my ($outfile) = sort keys %needed;
	    my $cmd = $out2cmd{$outfile} || die;

#	    printf STDERR "submit for slot %d: %s\n", $i, $cmd;
	    printf STDERR "submit for slot %d\n", $i;

	    my $c = new Cluster(
	    			"-outfile" => $outfile,
	    			"-project" => "PCGP",
	    		       );
	    $c->processors($bsub_cores) if $bsub_cores;
	    $c->node_class("");
	    $c->memory_reserve_mb($ram);
	    $c->memory_limit_mb($ram);
	    $c->command($cmd);
	    $c->run();

	    $slots[$i] = $outfile;
	    delete $needed{$outfile};
	    last unless %needed;
	  }
	}
      }

      my $elapsed = time - $start_time;
      my $wait = $sleep_time - $elapsed;
      $wait = 0 if $wait < 0;

      printf STDERR "have %d/%d (%.1f%%), waiting for %d sec...\n",
	$have, $total, ($have * 100 / $total),
	  $wait;
      sleep $wait if $wait;
    }
  }

  return $self->job_out_files(\@outfiles);
}

sub run_jobs_parallel {
  my ($self, %options) = @_;
  my $template = $options{"-template"} || die "-template";
  my $out_suffix = $options{"-out-suffix"} || die "-out-suffix";
  $out_suffix = "." . $out_suffix unless $out_suffix =~ /^\./;

  my (@cmds, @outfiles);

  my $gp_template = $template;
  $gp_template =~ s/\%s/{}/ || die "can't find parameter in template";
  # parallel will perform filename substitution for us

  if (0) {
    $gp_template .= ' >{}.out 2>{}.err';
    # also redirect stdout/stderr to logfiles
  }

  my $split_files = $self->split_files();

  foreach my $fn (@{$split_files}) {
    my $outfile = $fn . $out_suffix;
    push @outfiles, $outfile;
  }

  my $gp_cmd = sprintf "|parallel -v --eta";
#  $gp_cmd .= " -j+0";
#  $gp_cmd .= " -j-6";

  if (my $job_count = $options{"-pool"}) {
    $gp_cmd .= sprintf " -j %s", $job_count;
    # works with either exact number of jobs, or modifier like +1/-1
  } else {
    $gp_cmd .= " -j+0";
    # run one job per core
  }
  $gp_cmd .= sprintf " '%s'", $gp_template;

  if (1) {
    printf STDERR "parallel debug:\n";
    printf STDERR "  cmd: %s\n", $gp_cmd;
    foreach my $fn (@{$split_files}) {
      printf STDERR "  file:%s size:%d\n", $fn, -s $fn;
    }
    
  }

#  die $gp_cmd;

  my $fh = new FileHandle();
  $fh->open($gp_cmd) || die;

  foreach my $fn (@{$split_files}) {
    printf $fh "%s\n", $fn;
  }
  $fh->close;

  return $self->job_out_files(\@outfiles);
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
