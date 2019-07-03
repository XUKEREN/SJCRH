package FASTQ;

use strict;
use Configurable;
use FileHandle;
use Carp qw(cluck);

#use BinaryIndex;
# probably obsolete

#use misc;
use MiscUtils qw(average median);
# might be forgetting something

use POSIX qw(SEEK_SET);

@FASTQ::ISA = qw(Configurable Exporter);
@FASTQ::EXPORT_OK = qw();

use constant INDEX_TYPE_ALPHANUM => 1;
# standard: random string IDs
use constant INDEX_TYPE_SUBACCESSION => 2;
# all IDs are subaccessions of a single accession
# (e.g. NCBI SRA FASTQ files)

use MethodMaker qw(
	file
        fh
        index_type

        minimum_quality
        minimum_run_length
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->index_type(INDEX_TYPE_ALPHANUM);
  $self->configure(%options);
  return $self;
}

sub parse_file {
  # @SRR002786.1 3939:2:1:756:890 length=32
  # GGAGAACTTCTTTTAAGGTATGATATCACTGA
  # +SRR002786.1 3939:2:1:756:890 length=32
  # &'&+(<1%$$$$$$&$$$L$$$R$A$$$$$GI
  #
  # we can't just search for @ or + at start of line because
  # these might appear as quality values!
  #
  my ($self, %options) = @_;
  my $fh;
  if ($fh = $options{"-fh"}) {
  } elsif (my $file = $options{"-file"}) {
    $fh = new FileHandle();
    $fh->open($file) || die "can't open $file";
  }
  my $save_header = $options{"-save-header"};
  my $callback = $options{"-callback"};
  my $single = $options{"-single"};

  die "no filehandle" unless $fh;

  my $object;
  my @entries;
  my $data_ref;
  my $tell = $fh->tell();

  # parse
  # http://maq.sourceforge.net/fastq.shtml
  while (<$fh>) {
    chomp;
    next unless /\w/;

    if (/^\@(.*)/) {
      # ID line
      $object = {};
      $object->{id} = $1;
      $object->{tell} = $tell;
      $object->{header_line} = $_ if $save_header;
      push @entries, $object;

      if ($_ = <$fh>) {
	chomp;
	$object->{sequence} = $_;
	if ($_ = <$fh>) {
	  chomp;
	  if (/^\+(.*)/) {
	    my $qid = $1 || "";
	    if ($qid eq "" or $qid eq $object->{id}) {
	      # quality line ID matches
	      if ($_ = <$fh>) {
		chomp;
		$object->{quality} = $_;

		die "seq/quality length mismatch" unless length($object->{sequence}) == length $object->{quality};
	      } else {
		die "expected quality block";
	      }
	    } else {
	      die "ID mismatch, $_ vs " . $object->{id};
	    }
	  } else {
	    die "expected qual line, got $_";
	  }
	} else {
	  die "expected quality ID line";
	}
      } else {
	die "expected sequence line";
      }
    } else {
      die "expected ID line, got $_";
    }

    if ($callback) {
      my $result = &$callback($object);
      last if $result and $result eq "-1";
      # I forget what this is for
      @entries = ();
    } elsif ($single) {
      # return single record
      return $object;
    }

    $tell = $fh->tell();
  }

  $fh->close if (ref($fh) || "") eq "FileHandle";
  return \@entries;
}

sub encode_quality {
  # encode arrayref of phred-quality scores to FASTQ 
  my ($self, $qlist) = @_;
  return join "", map {chr($_ + 33)} @{$qlist};
}

sub decode_quality {
  my ($self, $string, $solexa) = @_;

  my @q;
  if ($solexa) {
    # http://maq.sourceforge.net/fastq.shtml
    # this code is given, but has too many parens!:
    # $Q = 10 * log(1 + 10 ** (ord($sq) - 64) / 10.0)) / log(10);
    # ?????
    
    # ...I think this is meant to convert the standard Solexa quality
    # range of -40 to +40.  Not sure WHAT some of the data I've seen is...

    foreach (split //, $string) { 
#      my $q = 10 * log(1 + 10 ** (ord($_) - 64) / 10.0) / log(10);
      # I don't think this is correct; returns very large values
      my $q = 10 * log(1 + 10 ** (ord($_) - 64) / 10.0) / log(10);
      push @q, int($q);
    }
  } else {
#    my $template = sprintf 'c%d', length $string;
    # 8/11/2009: hmm, for some reason unpacking with "c*" sometimes
    # leads to extra values!  Dunno why offhand, so just specify
#    @q = unpack $template, $string;
    # this also hosed if control characters are present,
    # seems to spill over into next base

    @q = map {unpack 'c', $_} split //, $string;

    if (@q != length $string) {
      die join ",", scalar(@q), length($string);
    }
    my $possibly_solexa = 0;
    foreach (@q) {
#      if ($_ > 90) {
#      if ($_ > 92) {
      if ($_ > 101) {
	# phred can return quality 68 (encoded as 101)
	$possibly_solexa = $_ if !$_ or $_ > $possibly_solexa;
      }

      # http://maq.sourceforge.net/fastq.shtml
      $_ -= 33;
    }
    printf STDERR "WARNING: FASTQ might be Solexa format! ($possibly_solexa)\n" if $possibly_solexa;
#    cluck sprintf "WARNING: FASTQ might be Solexa format! ($possibly_solexa)\n" if $possibly_solexa;
  }

  foreach (@q) {
    if ($_ < 0) {
      print STDERR "horrible error: negative quality score $_!\n";
      last;
    }
  }

  return \@q;
}

# sub find {
#   #
#   # find an identifier using binary index
#   #
#   require BinaryIndex;
#   my ($self, %options) = @_;
#   my $id = $options{"-id"} || die "-id";

#   $self->check_index();
#   # build index if necessary
#   my $bi = new BinaryIndex("-file" => $self->file);

#   my $index_id = $id;
#   if ($self->index_type() == INDEX_TYPE_SUBACCESSION) {
#     $id =~ /\.(\d+)$/ || die;
#     $index_id = $1;
#   }

#   my $ref;
#   if (my $pos = $bi->find("-id" => $index_id)) {
#     my $fh = $self->fh;
#     unless ($fh) {
#       $fh = new FileHandle();
#       $fh->open($self->file) || die $!;
#       $self->fh($fh);
#     }

#     $fh->seek($pos, SEEK_SET) || die "seek failed";

#     $ref = $self->parse_file(
# 			     "-fh" => $fh,
# 			     "-single" => 1,
# 			    );
#     die "search error for $id" unless $ref->{id} eq $id;
#   }
#   return $ref;
# }

# sub check_index {
#   my ($self) = @_;
#   my $index_type = $self->index_type();
#   my $bi = new BinaryIndex("-file" => $self->file);

#   if ($bi->index_needed()) {
#     $bi->set_index_type($index_type == INDEX_TYPE_SUBACCESSION ?
# 			BinaryIndex::INDEX_TYPE_INTEGER :
# 			BinaryIndex::INDEX_TYPE_ALPHANUM
# 		       );

#     my $last_acc;

#     my $callback = sub {
#       my ($ref) = @_;

#       my $id = $ref->{id};
#       if ($index_type == INDEX_TYPE_SUBACCESSION) {
# 	$id =~ /^(.*)\.(\d+)$/ || die "can't detect accession/subnumber";
# 	my ($main_acc, $sub) = ($1, $2);
# 	if (defined $last_acc) {
# 	  die "error: different accs $last_acc $main_acc, illegal w/this index type" unless $last_acc eq $main_acc;
# 	} else {
# 	  $last_acc = $main_acc;
# 	}
# 	$id = $sub;
#       }

#       $bi->add(
# 		"-id" => $id,
# 		"-position" => $ref->{tell}
# 	       );
#     };

#     $self->parse_file(
# 		      "-file" => $self->file,
# 		      "-callback" => $callback
# 		     );

#     $bi->build_index();
#   }
# }

sub get_header {
  my ($self, $ref) = @_;
  my $header = $ref->{header_line};
  if ($header) {
    $header =~ s/^[\@\+]//;
    $header =~ s/ length=\d+$//;
  } else {
    $header = $ref->{id};
  }
  return $header;
}

sub write_sequence {
  #
  # write a single sequence in FASTQ format
  #
  my ($self, %options) = @_;
  my $ref = $options{"-ref"} || die "-ref";
  my $fh = $options{"-fh"} || select();
  $fh = \*main::STDOUT unless ref $fh;

  my $len = length $ref->{sequence};
  die unless $len == length $ref->{quality};

  my $header = $self->get_header($ref);
  
  printf $fh '@%s' . "\n", $ref->{id};
  printf $fh "%s\n", $ref->{sequence};
  print $fh "+\n";
  printf $fh "%s\n", $ref->{quality};
}

sub write_fasta {
  #
  # write a single sequence in FASTA format
  #
  my ($self, %options) = @_;
  my $ref = $options{"-ref"} || die "-ref";
  my $fh = $options{"-fh"} || select();
  $fh = \*main::STDOUT unless ref $fh;
  my $header = $self->get_header($ref);
  printf $fh ">%s\n", $header;
  printf $fh "%s\n", $ref->{sequence};
}

sub quality_trim_run {
  #
  #  mask low-quality bases and return fragment of minimum
  #  continuous run length
  #
  my ($self, %options) = @_;
  my $ref = $options{"-ref"} || die;
  
  my $min_quality = $self->minimum_quality() || die "minimum_quality not set";
  # minimum continuous quality in run to be viable
  my $min_mer_length = $self->minimum_run_length() || die "minimum_run_length not set";

  my $seq = $ref->{sequence};
  my $qual_packed = $ref->{quality};
  my $qual = $self->decode_quality($qual_packed, 0);

  my $break_char = "\001";
  my $len = length $seq;
  my $ok;
  my $i;
  for ($i=0; $i < $len; $i++) {
    $ok = 1;
    if (substr($seq, $i, 1) !~ /[acgt]/i) {
      # missing basecall
#      printf STDERR "missing base %s\n", substr($seq,$i,1);
      $ok = 0;
    } elsif ($qual->[$i] < $min_quality) {
#      printf STDERR "poor quality %s\n", $qual->[$i];
      $ok = 0;
    }

    if (!$ok) {
      substr($seq, $i, 1) = $break_char;
      substr($qual_packed, $i, 1) = $break_char;
    }
  }

  my @seq_chunks = split /$break_char/, $seq;
  my @qual_chunks = split /$break_char/, $qual_packed;
  die unless @seq_chunks == @qual_chunks;
    
  my $longest=0;
  my $longest_i;
  for ($i=0; $i < @seq_chunks; $i++) {
    if (!defined($longest_i) or
	length($seq_chunks[$i]) > $longest) {
      $longest_i=$i;
      $longest=length $seq_chunks[$i];
    }
  }

  my $result;
  if ($longest >= $min_mer_length) {
    # usable
    my %r = %{$ref};
    $r{sequence} = $seq_chunks[$longest_i];
    $r{quality} = $qual_chunks[$longest_i];
    $result = \%r;
  } else {
    my $masked = $seq;
    $masked =~ s/$break_char/N/g;
    printf STDERR "rejecting sequence %s, max viable length=%d, masked=%s\n", $ref->{id}, $longest, $masked;
  }
  return $result;
}

sub quality_trim_end {
  #
  #  quality-trim read from end (where lowest-quality nucleotides
  #  are most likely to appear)
  #
  my ($self, %options) = @_;
  my $ref = $options{"-ref"} || die;
  my $verbose = $options{"-verbose"};
#  printf STDERR "trim: id=%s\n", $ref->{id};

#  my $WINDOW_SIZE = 5;
  my $WINDOW_SIZE = 6;
  # 6+ required to continue trim of SRR002786.15
#  my $WINDOW_SIZE = 7;
  my $WINDOW_CHECK_FRACTION = 0.20;
  # 0.1 = lowest decile 
  # 0.2 = lowest quintile
  # 0.5 = median

  my $HARD_WINDOW_SIZE = 3;
  my $HARD_WINDOW_CHECK_FRACTION = 0.5;
  my $hard_check_index = int($HARD_WINDOW_SIZE * $HARD_WINDOW_CHECK_FRACTION);

  my $WINDOW_CHECK_QUAL_THRESHOLD = $options{"-quality-threshold"} || 1;

#  my $WINDOW_CHECK_QUAL_THRESHOLD = 2;
#  my $WINDOW_CHECK_QUAL_THRESHOLD = 3;

  my $check_index = int($WINDOW_SIZE * $WINDOW_CHECK_FRACTION);

  my $ABS_MASK_QUALITY = 0;
  # threshold quality value where we trim regardless of window size
  
  my $seq = $ref->{sequence};
  my $qual_packed = $ref->{quality};
  my $qual = $self->decode_quality($qual_packed, 0);

  my $len = length $seq;

  my $i;

  my (@q, @qs, $q);

  my $trim_start = 0;
  my $trim_end = $len - 1;

  foreach my $direction (1, -1) {
    #
    # two passes: from start (forward) and end (reverse)
    #
    my $ok;
    my $stop = $direction == 1 ? $len - $WINDOW_SIZE : $WINDOW_SIZE - 1;
    my %taint;

    for ($i = $direction == 1 ? 0 : $len - 1;
	 $direction == 1 ? $i <= $stop : $i >= $stop;
	 $i += $direction) {
      #
      # this section:
      #  - first base in window flunking minimum quality threshold
      #  - main lookahead window test
      #
      $ok = 1;
      if ($direction == 1) {
	@q = @{$qual}[$i .. $i + ($WINDOW_SIZE - 1)];
      } else {
	@q = @{$qual}[$i - ($WINDOW_SIZE - 1) .. $i];
      }
      @qs = sort {$a <=> $b} @q;
      $q = $qs[$check_index];
      if (substr($seq, $i, 1) !~ /[acgt]/i) {
	# invalid basecall at window start position
	#      printf STDERR "missing base %s\n", substr($seq,$i,1);
	$ok = 0;
      } elsif ($direction == 1 and $q[0] <= $ABS_MASK_QUALITY) {
	# first base in window is below absolute trim level (forward)
	$ok = 0;
      } elsif ($direction == -1 and $q[$#q] <= $ABS_MASK_QUALITY) {
	# first base in window is below absolute trim level (reverse)
	$ok = 0;
      } elsif ($q < $WINDOW_CHECK_QUAL_THRESHOLD) {
	$ok = 0;
      }
      $taint{$i} = 1 unless $ok;
      printf STDERR "dir=%s window:%s  check:%s  median:%s  mean:%s  ok:%s\n", $direction, join(",", @q), $q, median(\@q), average(\@q), $ok if $verbose;

      #
      #  this section: hard window check
      #
      my $hard_ok = 1;
      my ($wstart, $wend);
      if ($direction == 1) {
	$wstart = $i;
	$wend = $i + ($HARD_WINDOW_SIZE - 1);
      } else {
	$wstart = $i - ($HARD_WINDOW_SIZE - 1);
	$wend = $i;
      }
      @q = @{$qual}[$wstart .. $wend];
      @qs = sort {$a <=> $b} @q;
      $q = $qs[$hard_check_index];
      printf STDERR "dir=%s window:%s q_fixed:%s\n", $direction, join(",", @q), $q if $verbose;
      if ($q < $WINDOW_CHECK_QUAL_THRESHOLD) {
	#
	# taint entire window UNLESS inner edge exceeds threshold
	#
	# examples using window size 3 with check at base 2/3:
	#
	#  1. pattern is "0 3 0"; trim entire window.
	#     Since we know the 3 is flanked by LQ sequence we're
	#     comfortable trimming it.
	#
	#  2. trimming in reverse and pattern is "3 0 0"
	#     Here we DON'T know the 3 is flanked by LQ sequence;
	#     it should be left alone for now.
	#
	$hard_ok = 0;
	my $j;
	my $trim;
	if ($direction == 1) {
	  # start trimming from FIRST instance of LQ sequence
	  # at end of current window, to beginning of window
	  for ($j = $wend; $j >= $wstart; $j--) {
	    $trim = 1 if $qual->[$j] < $WINDOW_CHECK_QUAL_THRESHOLD;
	    $taint{$j} = 1 if $trim;
	  }
	} else {
	  for ($j = $wstart; $j <= $wend; $j++) {
	    $trim = 1 if $qual->[$j] < $WINDOW_CHECK_QUAL_THRESHOLD;
	    $taint{$j} = 1 if $trim;
	  }
	}
      }

      last if $ok and $hard_ok;
    }

    if (%taint) {
      my @s = sort {$a <=> $b} keys %taint;
      if ($direction == 1) {
	# trimming forward: next base after final tainted base
	$trim_start = $s[$#s] + 1;
      } else {
	# trimming reverse: first base before final tainted base
	$trim_end = $s[0] - 1;
      }
    }

  }

  my %r = %{$ref};

  if ($trim_start >= $trim_end) {
    # entire sequence trimmed
    $r{sequence} = $r{quality} = "";
  } else {
    my $tlen = ($trim_end - $trim_start) + 1;
    $r{sequence} = substr($ref->{sequence}, $trim_start, $tlen);
    $r{quality} = substr($ref->{quality}, $trim_start, $tlen);
  }
  return \%r;
}


1;

