package BAMSet2;

use strict;
use Carp qw(confess cluck);

use File::Basename;
use File::Spec;

use Configurable;

use SampleName;

@BAMSet2::ISA = qw(Configurable Exporter);
@BAMSet2::EXPORT_OK = qw();

use MethodMaker qw(
		    bucket
		    disease
		    auto_dash_munge
		    ignore_symlinks
resolve_symlinks
unique

warn
ignore_bam_string
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->auto_dash_munge(1);
#  $self->ignore_symlinks(1);
  # usually used to link old/brokenly-named BAMs to fixed filenames
  $self->resolve_symlinks(1);
  # for symlinks, parse the target rather than the link name
  # (hopefully the the target conforms to sample naming convention)

  $self->bucket({});
  $self->unique({});
  $self->configure(%options);
  return $self;
}

sub add_bams {
#  cluck "add_bams()\n";
  my ($self, %options) = @_;
  my $bams = $options{"-bams"} || [];
  my $dir = $options{"-dir"};
  my $regexp_subject = $options{"-regexp-subject"};
  if ($options{"-standard"}) {
    my $disease = $self->get_disease(%options) || die "-disease";
    $dir = sprintf '/nfs_exports/genomes/1/PCGP/BucketRaw/SJ%s/', $disease;
  } elsif ($options{"-excap"}) {
    my $disease = $self->get_disease(%options) || die "-disease";
#    $dir = sprintf '/nfs_exports/genomes/1/PCGP/BucketRaw/SJ%s/', $disease;
    $dir = sprintf '/nfs_exports/genomes/1/projects/FREQEXCAP/PCGP/BucketRaw/SJ%s/', $disease;
  } elsif ($options{"-wgs"}) {
    my $disease = $self->get_disease(%options) || die "-disease";
#    $dir = sprintf '/nfs_exports/genomes/1/PCGP/BucketRaw/SJ%s/', $disease;
    $dir = sprintf '/nfs_exports/genomes/1/projects/WHOLEGENOME/PCGP/BucketRaw/SJ%s/', $disease;
  }

  if ($dir) {
    die "where is $dir" unless -d $dir;
    my @bams = glob($dir . "/*.bam");
#    printf STDERR "DEBUG: raw BAM list:\n%s\n", join "\n", @bams;
    die "no BAMs" unless @bams;
    my $ignore_string = $self->ignore_bam_string;
    foreach my $b (@bams) {
      next if $ignore_string and index($b, $ignore_string) > -1;
      push @{$bams}, $b;
    }
  }
  die "no BAMs" unless @{$bams};

  my $bucket = $self->bucket();
  my $DEBUG = 0;
  foreach my $bam (@{$bams}) {
    print STDERR "debug bam = $bam\n" if $DEBUG;
    my $unique = $self->unique();
    if (-l $bam) {
      next if $self->ignore_symlinks();
      if ($self->resolve_symlinks) {
	# use the filename referred to by the symlink instead
	my $orig = $bam;
	$bam = readlink($bam);
	if (dirname($bam) eq ".") {
	  $bam = dirname($orig) . "/" . $bam;
	}
#	die $bam;
#	printf STDERR "raw:%s readlink:%s dir:%s cooked:%s\n", $orig, $bam, dirname($bam), File::Spec->canonpath($bam);
      }
    }

    unless (-e $bam) {
      printf STDERR "ignoring %s, doesn't exist (broken link?)\n", $bam;
      next;
    }
    # $ ls -l  /nfs_exports/genomes/1/PCGP/BucketRaw/SJINF/SJINF059*bam
    # lrwxrwxrwx 1 mparker zhanggrp            20 Jul 20  2012 /nfs_exports/genomes/1/PCGP/BucketRaw/SJINF/SJINF059_G-SJ-13.bam -> SJINF060_G-SJ-14.bam
    # broken link, target doesn't exist

    next if $unique->{basename($bam)};
    # process each target file only once,
    # otherwise we create duplicate BAM match problems later

#    print STDERR "processing $bam\n";
    $unique->{basename($bam)} = 1;

    my $bn = basename($bam);

    if ($self->auto_dash_munge and $bn =~ /SJ[A-Z]+\d+\-/) {
      $bn =~ s/(SJ[A-Z]+\d+)\-/$1_/ || die;
    }

    my %info = SampleName::parse($bn, "WARN");
    if ($DEBUG) {
      printf STDERR "debug for %s:\n", $bn;
      foreach (sort keys %info) {
	printf STDERR "  %s: %s\n", $_, $info{$_};
      }
    }
    unless (%info) {
      printf STDERR "skipping bad name $bn\n";
      next;
    }
    if ($options{"-regexp-subject"}) {
      $bn =~ /^(SJ[A-Z]+\d\d\d)/ || die;
      my $hack_subject = $1;
      printf STDERR "map %s => %s (was %s)\n", $bn, $hack_subject, $info{subject};
      $info{subject_orig} = $info{subject};
      $info{subject} = $hack_subject;
    }

    $info{filename} = $bam;
    $info{broad_type} = $self->get_broad_type(\%info);
    $info{strict_type} = $self->get_strict_type(\%info);
    if (0) {
      foreach (sort keys %info) {
	printf STDERR "%s: %s\n", $_, $info{$_};
      }
      print STDERR "\n";
    }
    push @{$bucket->{$info{disease}}{$info{subject}}}, \%info;
  }
}

sub get_disease {
  my ($self, %options) = @_;
  my $disease = $self->disease || $options{"-disease"} || die "-disease";
  $disease =~ s/^SJ//;
  # formal code doesn't include SJ
  return $disease;
}

sub get_subjects {
  my ($self, %options) = @_;
  my $disease = $self->get_disease(%options);
  my $bucket = $self->bucket;
  return [ sort keys %{$bucket->{$disease}} ];
}

sub get_bam {
  my ($self, %options) = @_;
  my $disease = $self->get_disease(%options);
  my $bucket = $self->bucket();
  my $subject = $options{"-subject"} || die "-subject";
  my $index = $options{"-index"};
  my $disambiguate_require = $options{"-disambiguate-require"};
  my $disambiguate_exclude = $options{"-disambiguate-exclude"};
  # attempt to resolve ambiguities w/string match

  my @hits;
  my $error;
  if (my $type = $options{"-type"}) {
    # note this will be a little fuzzy! e.g.
    #   SJRHB011_D-TB-10-2213.bam
    #   SJRHB011_E-TB-09-0827.bam
    # will both be parsed as type D.
    @hits = grep {$_->{type} eq $type} @{$bucket->{$disease}{$subject}};
    $error = "no matches for type $type" unless @hits;
    if (@hits > 1) {
      die "FIX ME: disambig" if $disambiguate_require or $disambiguate_exclude;
      foreach my $ref (@hits) {
	foreach (sort keys %{$ref}) {
	  printf "%s: %s\n", $_, $ref->{$_};
	}
	print "\n";
      }
      $error = sprintf "ambiguous BAMs for type %s: %s", $type, join ",", map {$_->{filename}} @hits;
    }
  } elsif (my $broad_type = $options{"-broad-type"}) {
    @hits = grep {$_->{broad_type} eq $broad_type} @{$bucket->{$disease}{$subject}};
    $error = "no matches for broad type $type" unless @hits;
    $error = "ambiguous BAMs for broad type $broad_type" if @hits > 1;
    die "FIX ME: disambig" if $disambiguate_require or $disambiguate_exclude;
  } elsif (my $strict_type = $options{"-strict-type"}) {
    @hits = grep {$_->{strict_type} eq $strict_type} @{$bucket->{$disease}{$subject}};
    if ($index) {
      @hits = grep {$_->{index} == $index} @hits;
    }
    $error = "no matches for $subject and strict type $strict_type" unless @hits;
    if (0 and not(@hits)) {
      foreach my $subject (keys %{$bucket->{$disease}}) {
	foreach my $hash (@{$bucket->{$disease}{$subject}}) {
	  foreach my $k (sort keys %{$hash}) {
	    printf "%s: %s\n", $k, $hash->{$k};
	  }
	  print "\n";
	}
	die $subject;
      }
      die join "\n", $error, join ",", keys %{$bucket->{$disease}};
      die $error;
    }

    if (@hits > 1 and $disambiguate_require or $disambiguate_exclude) {
      my @filtered;
      foreach my $ref (@hits) {
	next if $disambiguate_require and index(basename($ref->{filename}), $disambiguate_require) == -1;
	next if $disambiguate_exclude and index(basename($ref->{filename}), $disambiguate_exclude) > -1;
	push @filtered, $ref;
      }
      @hits = @filtered;
    }

    if ($options{"-list"}) {
      # just return raw set, don't worry about duplicates
      return \@hits;
    }

    if (@hits > 1) {
      foreach my $ref (@hits) {
	foreach (sort keys %{$ref}) {
	  printf "%s: %s\n", $_, $ref->{$_};
	}
	print "\n";
      }
      $error = sprintf "ambiguous BAMs for strict type %s: %s", $strict_type, join ",", map {$_->{filename}} @hits;
    }

    $error .= sprintf " avail=%s", join ",", map {$_->{filename}} @{$bucket->{$disease}{$subject}} if $error;
  } elsif ($options{"-all"}) {
    @hits = @{$bucket->{$disease}{$subject}};
  } else {
    die "specify -type / -broad-type / -strict-type / -all";
  }
  my $result;
  if ($error) {
    if ($self->warn() or $options{"-warn"}) {
      print STDERR "$error\n";
    } else {
      confess $error;
    }
  } elsif ($options{"-all"}) {
    $result = [ map {$_->{filename}} @hits ];
  } else {
    $result = $hits[0]->{filename};
  }
  return $result;
}

sub get_broad_type {
  # generate broad sample type from annotations
  # e.g.
  #  - D and M => D
  #  - v1 D/E/F => D
  # http://hc-wiki.stjude.org:8080/display/DTAWHS/Sample+Name+Library
  my ($self, $info) = @_;
  my $version = $info->{conventionVersion} || die;
  my $result;
  my $type = $info->{type} || die;
  if ($version == 1) {
    if ($type eq "G") {
      $result = "G";
    } elsif (
	     $type eq "D" or
	     $type eq "M" or
	     $type eq "R" or
	     $type eq "X" or
	     # xenograft: I guess they wouldn't transplant normal cells...
	     $type eq "A"
	     # autopsy
	    ) {
      $result = "D";
    } else {
      die "unhandled v1 type code $type";
    }
  } elsif ($version == 2) {
    if ($type eq "G") {
      $result = "G";
    } elsif (
	     $type eq "D" or
	     $type eq "R" or
	     $type eq "M" or
	     $type eq "C" or
	     $type eq "X" or
	     # xenograft: I guess they wouldn't transplant normal cells...
	     $type eq "A"
	    ) {
      $result = "D";
    } else {
      die "unhandled v2 type code $type";
    }
  } else {
    die "BAM version $version, FIX ME";
  }
  return $result;
}

sub get_strict_type {
  # generate strict sample type from annotations
  my ($self, $info) = @_;
  my @stuff = split /_/, $info->{sample};
  die unless @stuff == 2;
#  die join ",", @stuff, substr($stuff[1], 0, 1) unless length $stuff[1] == 1;
  return substr($stuff[1], 0, 1);
  # i.e. for SJCBF001018_D2-TB-00-5223.bam return D, not D2
}

sub get_alternate_disease_codes {
  # return all non-D disease codes in BAMs
  my ($self, %options) = @_;
  my $disease = $self->get_disease(%options);
  my $bucket = $self->bucket();
  my $subject = $options{"-subject"} || die "-subject";

  my %codes;

  foreach my $ref (@{$bucket->{$disease}{$subject}}) {
    if ($ref->{type} eq "G") {
      printf STDERR "skipping %s: type=%s\n", $ref->{filename}, $ref->{type};
      next;
    }
    next if $ref->{strict_type} eq "D";
    $codes{$ref->{strict_type}} = 1;
#    printf "%s: %s %s\n", @{$ref}{qw(filename type strict_type)};
  }
  return [sort keys %codes];
}

sub get_diseases {
  my ($self) = @_;
  return [ sort keys %{$self->bucket()} ];
}


1;
