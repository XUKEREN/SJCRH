package LiftOver;

use strict;
use File::Copy;

use Configurable;
use TemporaryFileWrangler;
use TdtConfig;
use GenomeUtils qw(cook_chromosome_name);

@LiftOver::ISA = qw(Configurable Exporter);
@LiftOver::EXPORT_OK = qw();

use MethodMaker qw(
		   liftover_binary
		   liftover_chain_dir
                   single_chain_dir

		   translated_ok
		   translated_chr
		   translated_base

                   translated_start_interbase
                   translated_end_interbase
                   translated_min_match

                   enable_retry
error
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
#  $self->liftover_binary("/app1/dnload/newersrc/bin/x86_64/liftOver");
  $self->liftover_binary("liftOver");
  # expected on PATH
  $self->enable_retry(1);
  # enable retrying with looser criteria on failure
  $self->configure(%options);

  # $self->liftover_chain_dir("/h1/edmonsom/liftover/");
  # $self->liftover_chain_dir("/nfs_exports/apps/gnu-apps/NextGen/lwbin/database/CrossAssembly/LiftOver/chains/");
  # $self->single_chain_dir(1);

  # $self->liftover_chain_dir("/nfs_exports/genomes/1/LIFTOVER_CHAIN/");
  unless ($self->liftover_chain_dir()) {
    my $config = TdtConfig::readConfig('app', 'liftover') || die "can't find liftover app config";
    $self->liftover_chain_dir($config->{LIFTOVER_CHAIN_DIR} || die "no chain dir in config");
  }
  $self->single_chain_dir(0);

  return $self;
}

sub translate_base {
  my ($self, %options) = @_;
  #
  # translate base number from one genome build to another.
  # base number is 1-based, not 0-based.
  #
  my $from = $options{"-from"} || die "-from";
  my $to = $options{"-to"} || die "-to";
  my $chr = $options{"-chr"} || die "-chr";

  my ($start_interbase, $end_interbase);

  my $interval_mode;
  if (my $base = $options{"-base"}) {
    # 1-based
    $start_interbase = $base - 1;
    $end_interbase = $base;
  } elsif ($start_interbase = $options{"-start-interbase"}) {
    $end_interbase = $options{"-end-interbase"} || die;
    $interval_mode = 1;
  } else {
    die;
  }

  my $tfw = new TemporaryFileWrangler();
#  my $bn = tmpnam() . ".liftover";
  my $bn = $tfw->get_tempfile() . ".liftover";

  my $fn_in = $bn . ".in";
  my $fn_out = $bn . ".out";
  my $fn_unmapped = $bn . ".unmapped";

  my $binary = $self->liftover_binary() || die "no liftover binary";
#  die "$binary not executable" unless -x $binary;

  my $fn_chain;
  if ($self->single_chain_dir) {
    $fn_chain = sprintf '%s/%sTo%s.over.chain',
    $self->liftover_chain_dir(),
    $from,
    ucfirst($to);
  } else {
    # subdirectories based on "from" genome
    $fn_chain = sprintf '%s/%s/%sTo%s.over.chain',
    $self->liftover_chain_dir(),
    ucfirst($from),
    $from,
    ucfirst($to);
  }
#  die "can't find or read $fn_chain" unless -s $fn_chain and -r $fn_chain;
  die "can't find or read $fn_chain" unless -s $fn_chain;
  # -r might not work on sonas even though file actually readable

  #  /app1/dnload/newersrc/bin/x86_64/liftOver hg18_chr17_base_7519167.bed hg18ToHg19.over.chain out unmapped

  foreach ($from, $to) {
#    die "must be hg18 or hg19" unless $_ eq "hg18" or $_ eq "hg19";
    die "genome code \"$_\" formatting problem: must be lc" unless /^[a-z]+\d+$/;
  }

  my $ok = 0;
  my $trans_chr = "";
  my $trans_base = 0;
  my $trans_interbase_start = -1;
  my $trans_interbase_end = -1;
  my $min_match;
  my $problem = "";

  my $cc = cook_chromosome_name($chr, "-undef-unknown" => 1);
  if ($cc) {
    $chr = "chr" . $chr unless $chr =~ /^chr/;
    $chr =~ /^chr(\w+)$/ || die;
    die "formatting problem for $chr" if length($1) > 2;
    open(LOTMP, ">" . $fn_in) || die "can't write to $fn_in";
    printf LOTMP "%s\n", join "\t", $chr, $start_interbase, $end_interbase;
    close LOTMP;
  } else {
    $problem = "incompatible chromosome $chr";
  }

  while (1) {
    foreach my $f ($fn_out, $fn_unmapped) {
      unlink $f;
      die if -s $f;
    }

    last unless defined $cc;

    my $cmd = join " ", $binary, $fn_in, $fn_chain, $fn_out, $fn_unmapped;
    $cmd .= sprintf " -minMatch=%.2f", $min_match if $min_match;
    $cmd .= " >/dev/null 2>&1";

    if (0) {
      printf STDERR "cmd: %s\n", $cmd;
      copy($fn_in, "liftover_in.txt") || die;
    }
    system $cmd;
    if ($?) {
      print STDERR "ERROR: exit code $?\n";
      $ok = 0;
      die "error; is $binary on PATH?";
      # likely can't find binary / PATH problem
    } elsif (-s $fn_out) {
      open(LOTMP, $fn_out) || die;
      my $line = <LOTMP>;
      chomp $line;
      close LOTMP;
      my ($chr, $start, $end) = split /\t/, $line;
      $trans_chr = $chr;
      $trans_base = $start + 1;
      $trans_interbase_start = $start;
      $trans_interbase_end = $end;
      $ok = 1;
    } else {
      # no liftOver output file: base may have been deleted
      if (-s $fn_unmapped) {
	# verify
	open(UMCHK, $fn_unmapped) || die;
	while (<UMCHK>) {
	  if (/^#/) {
	    if (/Deleted in new/) {
	      $problem = "deleted";
	    } elsif (/Partially deleted in new/) {
	      $problem = "partially_deleted";
	    } elsif (/Split in new/) {
	      $problem = "split_in_new";
	    } elsif (/Duplicated in new/) {
	      $problem = "duplicated_in_new";
	    } else {
	      die sprintf "ERROR: unhandled liftOver result line for %s/%s/%s: %s", $chr, $start_interbase, $end_interbase, $_;
	    }
	  }
	}
	if ($problem) {
	  $ok = 0;
	  $trans_chr = $trans_base = $trans_interbase_start = $trans_interbase_end = $problem;
	} else {
	  printf STDERR "ERROR: unhandled liftOver unmapped case, command=%s\n", $cmd;
	  die "DEBUG ME";
	}
      } else {
	printf STDERR "ERROR: liftOver failed to generate output file, command=%s\n", $cmd;
	$ok = 0;
	die "TEST ME";
	# ???
      }
    }
    $self->error($problem);

    if ($ok) {
      # success: done
      $ok = 2 if $min_match;
      # required looser match
      last;
    } elsif ($interval_mode and $self->enable_retry()) {
      # retry intervals w/more relaxed params on failure
      $min_match = 1 unless defined $min_match;
#      $min_match -= 0.01;
      # odd rounding issues
      $min_match = sprintf '%.2f', $min_match - 0.01;
#      printf STDERR "retry %s\n", $min_match;
      last unless $min_match >= 0.01;
      # stop eventually
    } else {
      # SNV: only try once
      last;
    }
  }


  $self->translated_ok($ok);
  $self->translated_chr($trans_chr);
  $self->translated_base($trans_base);
  $self->translated_start_interbase($trans_interbase_start);
  $self->translated_end_interbase($trans_interbase_end);
  $self->translated_min_match($min_match || "");

  unlink($fn_in, $fn_out, $fn_unmapped);
  return $ok;
}

sub get_chain_file {
  my ($self, %options) = @_;
  my $from = $options{"-from"} || die "-from";
  my $to = $options{"-to"} || die "-to";
  foreach ($from, $to) {
#    die "must be hg18 or hg19" unless $_ eq "hg18" or $_ eq "hg19";
    die "genome code \"$_\" formatting problem: must be lc" unless /^[a-z]+\d+$/;
  }
  my $fn_chain;
  if ($self->single_chain_dir) {
    $fn_chain = sprintf '%s/%sTo%s.over.chain',
    $self->liftover_chain_dir(),
    $from,
    ucfirst($to);
  } else {
    # subdirectories based on "from" genome
    $fn_chain = sprintf '%s/%s/%sTo%s.over.chain',
    $self->liftover_chain_dir(),
    ucfirst($from),
    $from,
    ucfirst($to);
  }
#  die "can't find or read $fn_chain" unless -s $fn_chain and -r $fn_chain;
  die "can't find or read $fn_chain" unless -s $fn_chain;
  return $fn_chain;
}

1;
