package LiftOverBatch;
# batch liftOver w/recursion on failures
# MNE 10/2015

use strict;
use Configurable;
use Exporter;

use LiftOver;
use LiftOverStrand;
use TemporaryFileWrangler;
use FileUtils qw(read_simple_file);
use MiscUtils qw(split_list dump_die);
use File::Copy;
use File::Basename;

use constant F_LO_IN_CHR => "_FIC";
use constant F_LO_IN_START => "_FIS";
use constant F_LO_IN_END => "_FIE";
# input chromosome and interbase start/end
use constant F_LO_OUT_CHR => "_FOC";
use constant F_LO_OUT_START => "_FOS";
use constant F_LO_OUT_END => "_FOE";
# output chromosome and interbase start/end
use constant F_LO_OUT_MINMATCH => "_FOMM";
use constant F_LO_OUT_ERROR => "_FOERR";
use constant F_LO_OUT_SAME_STRAND => "FOSTR";

@LiftOverBatch::ISA = qw(Configurable Exporter);
@LiftOverBatch::EXPORT_OK = qw(
				F_LO_IN_CHR
				F_LO_IN_START
				F_LO_IN_END
				F_LO_OUT_CHR
				F_LO_OUT_START
				F_LO_OUT_END
				F_LO_OUT_ERROR
				F_LO_OUT_MINMATCH
				F_LO_OUT_SAME_STRAND
			     );

use MethodMaker qw(
	genome_from
        genome_to
	chain_file
	tfw

	enable_retry
	enable_strand_check
	los
enable_chromosome_cache
		 );
# enable_chromosome_cache override must be in constructor

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->enable_retry(1);
  $self->enable_strand_check(1);
  $self->enable_chromosome_cache(1);
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub genome_to_chain_string {
  # STATIC
  my ($genome) = @_;
  my %map = qw(
		GRCh37-lite   hg19
		GRCh38        hg38
		hg18          hg18
      );
  return $map{$genome};
}

sub setup {
  my ($self) = @_;
  my $genome_from = $self->genome_from() || die;
  my $genome_to = $self->genome_to() || die;

  my $str_from = genome_to_chain_string($genome_from) || die "can't get chain for $genome_from";
  my $str_to = genome_to_chain_string($genome_to) || die "can't get chain for $genome_to";
  my $lo = new LiftOver();
  my $chain = $lo->get_chain_file("-from" => $str_from, "-to" => $str_to);
  $self->chain_file($chain);
  $self->tfw(new TemporaryFileWrangler());

  if ($self->enable_strand_check) {
    $self->los(new LiftOverStrand(
				  "-genome_from" => $genome_from,
				  "-genome_to" => $genome_to,
				  "-chunk_cache_chromosomes" => $self->enable_chromosome_cache
				 ));
  }

}

sub query_batch {
  my ($self, %options) = @_;
  my $rows = $options{"-rows"} || die "-rows";
  my $min_match = $options{"-min-match"};
  my $min_match_floor = $options{"-min-match-floor"} || 0.01;
#  printf STDERR "starting set of %d, mm=%s\n", scalar(@{$rows}), $min_match || "n/a";

  # - prep input .bed file
  # - run liftOver
  # - check if completed OK for all rows

#  my $bn = tmpnam() . ".liftover";
  my $fn_chain = $self->chain_file || die;

  my $tfw = $self->tfw();
  my $bn = $tfw->get_tempfile("-append" => ".liftover");
  my $fn_in = $bn . ".in";
  my $fn_out = $bn . ".out";
  my $fn_unmapped = $bn . ".unmapped";

  open(LOTMP, ">" . $fn_in) || die;
  my %queue;
  foreach my $row (@{$rows}) {
    my $chr = $row->{F_LO_IN_CHR()} || die;
    my $start = $row->{F_LO_IN_START()} || die;
    my $end = $row->{F_LO_IN_END()} || die;
    dump_die($row, "no chr prefix") unless $chr =~ /^chr/;
    my $in_key = get_in_key($row);
    push @{$queue{$in_key}}, $row;
    printf LOTMP "%s\n", join "\t", $chr, $start, $end, $in_key;
  }
  close LOTMP;

  my $los = $self->los();
  my $strand_check = $self->enable_strand_check();
  my $enable_retry = $self->enable_retry();

  while (1) {
    unlink($fn_out, $fn_unmapped);
    # paranoia
    die if -s $fn_out or -s $fn_unmapped;

    my $cmd = sprintf 'liftOver %s %s %s %s',
    $fn_in, $fn_chain, $fn_out, $fn_unmapped;
    $cmd .= sprintf ' -minMatch=%.2f', $min_match if $min_match;
    $cmd .= " >/dev/null 2>&1";
#    printf STDERR "cmd: %s\n", $cmd;
    system($cmd);
    if ($?) {
      my $error = $?;
      my $f_backup = $ENV{HOME} . "/" . basename($fn_in);
      copy($fn_in, $f_backup);
      my $cmd2 = sprintf 'df -k %s', dirname($fn_in);
      print STDERR "running $cmd2\n";
      system $cmd2;
      die "ERROR: exit $error from command:$cmd input copied to $f_backup";
    }


    #
    #  parse successful rows:
    #
    my $results = read_simple_file($fn_out);
    foreach my $rr (@{$results}) {
      my @f = split /\t/, $rr;
      die scalar @f unless @f == 4;
      my ($out_chr, $out_start, $out_end, $in_key) = @f;
#      print STDERR "processing $in_key\n";
      my $src_rows = $queue{$in_key} || next;
      # won't be present if already handled

      my $strand_ok;
      if ($strand_check) {
	$strand_ok = $los->strand_check(
	  "-from-chr" => $src_rows->[0]->{F_LO_IN_CHR()},
	  "-from-start" => $src_rows->[0]->{F_LO_IN_START()},
	  "-from-end" => $src_rows->[0]->{F_LO_IN_END()},
	  "-to-chr" => $out_chr,
	  "-to-start" => $out_start,
	  "-to-end" => $out_end);
      }

      foreach my $row (@{$src_rows}) {
	$row->{F_LO_OUT_CHR()} = $out_chr;
	$row->{F_LO_OUT_START()} = $out_start;
	$row->{F_LO_OUT_END()} = $out_end;
	$row->{F_LO_OUT_ERROR()} = "";
	$row->{F_LO_OUT_SAME_STRAND()} = $strand_ok;
	$row->{F_LO_OUT_MINMATCH()} = $min_match if defined $min_match;
      }
      delete $queue{$in_key};
    }

    if (-s $fn_unmapped) {
      my $um = read_simple_file($fn_unmapped);
      while (@{$um}) {
	my $reason = shift @{$um};
	my $line = shift @{$um};
	my @f = split /\t/, $line;
	die unless @f == 4;
	my $in_key = $f[3];
#	print STDERR "$reason $line\n";
	my $problem;
	my $retryable = 0;
	if ($reason =~ /Deleted in new/) {
	  $problem = "deleted";
	} elsif ($reason =~ /Partially deleted in new/) {
	  $problem = "partially_deleted";
	  $retryable = 1;
	  $retryable = 0 if $min_match and $min_match <= $min_match_floor;
	} elsif ($reason =~ /Split in new/) {
	  $problem = "split_in_new";
	  $retryable = 1;
	  $retryable = 0 if $min_match and $min_match <= $min_match_floor;
	} elsif ($reason =~ /Duplicated in new/) {
	  $problem = "duplicated_in_new";
	} else {
	  die sprintf "ERROR: unhandled liftOver result line $reason";
	}

	#
	#  record errors for each unmappable variant:
	#
	if (my $src_rows = $queue{$in_key}) {
	  # unless already processed
	  foreach my $row (@{$src_rows}) {
	    foreach my $f (F_LO_OUT_CHR(),
			   F_LO_OUT_START(),
			   F_LO_OUT_END(),
			   F_LO_OUT_SAME_STRAND(),
			   F_LO_OUT_MINMATCH()
			  ) {
	      $row->{$f} = "";
	    }
	    $row->{F_LO_OUT_ERROR()} = $problem;
	  }
	}
	delete $queue{$in_key} unless $retryable;
      }
    }

    if (%queue) {
      if ($enable_retry) {
	$min_match = 0.95 unless $min_match;
	# 0.95 = liftOver default
	$min_match -= 0.01;
#	printf STDERR "retrying: queue=%d min_match=%.2f\n", scalar(keys %queue), $min_match;
      } else {
	# not retrying, leave failures alone
	last;
      }
    } else {
      # all complete or unretryable
      last;
    }
  }

  unlink($fn_in, $fn_out, $fn_unmapped);
  return;

  die;

  my $usable = 1;
  my $problem;
  if (-s $fn_unmapped) {
    $usable = 0;
    $problem = "has_unmapped";
  }
  my $results = read_simple_file($fn_out);
  unless (@{$results} == @{$rows}) {
    $usable = 0;
    $problem = "count_mismatch" unless $problem;
  }

  if ($usable) {
    #
    # OK
    #
    for (my $i = 0; $i < @{$rows}; $i++) {
      my @f = split /\t/, $results->[$i];
      die scalar @f unless @f == 3;
      my ($out_chr, $out_start, $out_end) = @f;
      my $row = $rows->[$i];
      $row->{F_LO_OUT_CHR()} = $out_chr;
      $row->{F_LO_OUT_START()} = $out_start;
      $row->{F_LO_OUT_END()} = $out_end;
      $row->{F_LO_OUT_ERROR()} = "";

      if ($strand_check) {
#	printf STDERR "LOB debug: %s\n", $row->{accession};
	my $strand_ok = $los->strand_check(
					   "-from-chr" => $row->{F_LO_IN_CHR()},
					   "-from-start" => $row->{F_LO_IN_START()},
					   "-from-end" => $row->{F_LO_IN_END()},

					   "-to-chr" => $out_chr,
					   "-to-start" => $out_start,
					   "-to-end" => $out_end
					  );
	$row->{F_LO_OUT_SAME_STRAND()} = $strand_ok;
      }

    }

    if ($min_match and @{$rows} == 1) {
      # success retrying single row with lower -minMatch
      $rows->[0]->{F_LO_OUT_MINMATCH()} = $min_match;
    }

  } else {
    #
    # problem somewhere in list:
    #
    if (@{$rows} == 1) {
      #
      # problematic row isolated
      #
      my $problem = "";
      my $retryable = 0;
      if (-s $fn_unmapped) {
	open(UMCHK, $fn_unmapped) || die;
	while (<UMCHK>) {
	  if (/^#/) {
	    if (/Deleted in new/) {
	      $problem = "deleted";
	    } elsif (/Partially deleted in new/) {
	      $problem = "partially_deleted";
	      $retryable = 1;
	    } elsif (/Split in new/) {
	      $problem = "split_in_new";
	    } elsif (/Duplicated in new/) {
	      $problem = "duplicated_in_new";
	    } else {
	      die sprintf "ERROR: unhandled liftOver result line $_";
	    }
	  }
	}
      }
#      dump_die($rows->[0], "problem=$problem retryable=$retryable", 1);

      my $retry = $self->enable_retry;
      $retry = 0 unless $retryable;
      $min_match = 0.95 unless $min_match;
      # 0.95 = liftOver default
      $min_match -= 0.01;
      $retry = 0 if $min_match < 0.01;

      if ($retry) {
#	dump_die($rows->[0], "retry $min_match", 1);
	$self->query_batch(
			   "-rows" => $rows,
			   "-min-match" => $min_match
			  );
      } else {
	# can't retry: stop
	my $row = $rows->[0];
	$row->{F_LO_OUT_CHR()} = "";
	$row->{F_LO_OUT_START()} = "";
	$row->{F_LO_OUT_END()} = "";
	$row->{F_LO_OUT_SAME_STRAND()} = "";
	$row->{F_LO_OUT_ERROR()} = $problem;
      }
    } else {
      # recurse
      my $next_count = int(scalar @{$rows} / 2);
      $next_count = 1 if $next_count < 1;

      my $lists = split_list($rows, $next_count);
      foreach my $list (@{$lists}) {
	$self->query_batch("-rows" => $list);
      }
    }
  }

  unlink($fn_in, $fn_out, $fn_unmapped);
}

sub get_in_key {
  # STATIC
  my ($row) = @_;
  return join ".", @{$row}{F_LO_IN_CHR(),
			   F_LO_IN_START(),
			   F_LO_IN_END()};
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
