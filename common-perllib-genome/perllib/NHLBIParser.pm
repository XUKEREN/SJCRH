package NHLBIParser;
# parse NHLBI flatfiles and perform basic munging
# e.g.
# /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/EVS_SNPs/ESP6500/SNP_txt/
# MNE 11/2013
#
# adds helper fields:
# - reference_sequence
# - reference_position
# - reference_base
# - variant_bases (arrayref!)
# - aa_list (arrayref!)
# - genes (arrayref!)
# - raw_line (for rewriting)

use strict;
use FileHandle;

use Configurable;
use MiscUtils qw(dump_die);

@NHLBIParser::ISA = qw(Configurable);

use MethodMaker qw(
	file
fh
headers
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
  if ($self->file) {
    # might not be specified, e.g. parse via tabix hit
    my $fh = new FileHandle();
    $fh->open($self->file || die "-file") || die;
    $self->fh($fh);
  }
}

sub next {
  my ($self, %options) = @_;
  my $fh = $self->fh() || die;
  my $headers = $self->headers();
  unless ($headers) {
    while (<$fh>) {
      if (/^\#\#/) {
	# secondary? info/headers
	next;
      } elsif (/^\#/) {
	# header line
	s/^\#//;
	$self->headers($headers = [ split /\s+/, $_ ]);
	last;
      } else {
	die;
      }
    }
  }

  my $line = <$fh>;
  return defined($line) ? $self->parse_line($line) : undef;
}

sub get_genotype_total {
  my ($self, $row, $ref_base, $var_base) = @_;
#  printf STDERR "ref:%s var:%s raw:%s\n", $ref_base, $var_base, $row->{AllGenotypeCount};
  my @all_g = split /\//, $row->{AllGenotypeCount};
  die unless @all_g;
  my $g_homozygous_ref = $ref_base x 2;
  die unless $g_homozygous_ref =~ /^[A-Z][A-Z]$/;
  my $total_patients = 0;
  my $hits = 0;
  my %wanted = map {$_, 1} ($ref_base, $var_base);
  my $population_total = 0;
  foreach my $gs (@all_g) {
    my ($genotype, $count) = split /\=/, $gs;
    $population_total += $count;

    my $found_variant;
    my $found_other;
    my @gb = (split //, $genotype);
    if (@gb == 2) {
      # expected 2-base genotype
      foreach my $b (@gb) {
	$found_variant = 1 if $b eq $var_base;
	$found_other = 1 unless $wanted{$b};
      }
    } elsif (@gb == 1) {
      printf STDERR "ignoring wack single-base NHLBI genotype %s, raw=%s\n", $genotype, $row->{AllGenotypeCount};
      # FIX ME!
    } else {
      die "unknown genotype $genotype";
    }
    if ($found_variant and not($found_other)) {
      # - genotype contains desired variant (homozygous or heterozygous)
      # - genotype does not include a base other than the variant
#      print STDERR "ok $genotype\n";
      $total_patients += $count;
      $hits++;
    }
  }
  printf STDERR "AllGenotypeCount:%s total_variant_patients:%d EA:%s AA:%s\n",
  $row->{AllGenotypeCount}, $population_total,
  $row->{EuropeanAmericanGenotypeCount},
  $row->{AfricanAmericanGenotypeCount};

  die "WTF $hits" unless $hits == 2;
  return $total_patients;
}

sub get_genotype_frequency {
  my ($self, $row, $var_base, %options) = @_;
#  printf STDERR "ref:%s var:%s raw:%s\n", $ref_base, $var_base, $row->{AllGenotypeCount};
  die "var base $var_base" unless $var_base =~ /^[A-Z]$/;
  my $pkey = sprintf '%sGenotypeCount', $options{"-population"} || "All";
  dump_die($row, "no key $pkey") unless exists $row->{$pkey};
  my @all_g = split /\//, $row->{$pkey};
  die unless @all_g;
  my $maf_mode = $options{"-genotype"};

  my $freq;
  my $verbose;

  if ($maf_mode) {
    my $homozygous_alt_double = $var_base x 2;
    my $variant_alleles = 0;
    my $total_alleles = 0;
    my $check;
    printf STDERR "genotype start: %s\n", join " ", @all_g if $verbose;
    foreach my $gs (@all_g) {
      my ($genotype, $subject_count) = split /\=/, $gs;
      printf STDERR "genotype:%s count:%d\n", $genotype, $subject_count if $verbose;

      if (length($genotype) == 2) {
	# diploid
	if ($genotype eq $homozygous_alt_double) {
	  # homozygous for variant
	  $variant_alleles += $subject_count * 2;
	  printf STDERR "  adding %d\n", $subject_count * 2 if $verbose;
	} elsif (index($genotype, $var_base) != -1) {
	  # heterozygous for variant
	  $variant_alleles += $subject_count;
	  printf STDERR "  adding %d\n", $subject_count if $verbose;
	}
	$total_alleles += $subject_count * 2;
      } elsif (length($genotype) == 1) {
	# chrX.7812017.C.T
	# TT=1006/TC=1047/T=1059/CC=1794/C=1072
	if ($genotype eq $var_base) {
	  $variant_alleles += $subject_count;
	  printf STDERR "  adding %d\n", $subject_count if $verbose;
	  printf STDERR "WARNING: single allele genotype\n";
	}
	$total_alleles += $subject_count;
	$check = 1;
      } else {
	die length $genotype;
      }
    }

    $freq = $variant_alleles / $total_alleles;
    printf STDERR "var_alleles:%d  total_alleles:%d  MAF:%.2f\n", $variant_alleles, $total_alleles, $freq if $verbose;
#    die if $check;
  } else {
    # original version: FAIL?  i.e. population frequency rather than MAF?
    my $variant_genotypes = 0;
    my $total_genotypes = 0;
    foreach my $gs (@all_g) {
      my ($genotype, $count) = split /\=/, $gs;
      $total_genotypes += $count;
      $variant_genotypes += $count if index($genotype, $var_base) != -1;
    }
    $freq = $variant_genotypes / $total_genotypes;
  }

#  printf STDERR "AllGenotypeCount:%s variant:%d total:%d freq:%f\n", $row->{AllGenotypeCount}, $variant_genotypes, $total_genotypes, $freq;

  return $freq;
}

sub parse_line {
  my ($self, $line) = @_;
  my $row;
  if (ref $line) {
    # already-parsed row (e.g. tabix hit)
    $row = $line;
  } else {
    chomp $line;
    my $headers = $self->headers() || die;
    $row = {};
    @{$row}{@{$headers}} = split /\s+/, $line;
    $row->{raw_line} = $line;
  }

  my $alleles_raw = $row->{Alleles};
  my @all_alleles;
  my $original_format;
  if ($alleles_raw =~ /\//) {
    # older format
#    @all_alleles = split /\//, $row->{Alleles} || die sprintf "can't parse Alleles field: alleles=%s raw=%s", $row->{Alleles}, $line;
    @all_alleles = split /\//, $alleles_raw;
    $original_format = 1;
  } elsif ($alleles_raw =~ />/) {
    # might be multiple entries, e.g.
    # TA>TAAAAA;TA>TAA;TA>TAAA;TA>T
    my @things = split /;/, $alleles_raw;
    my %saw;
    foreach my $thing (@things) {
      my @f = split />/, $thing;
      die unless @f == 2;
      foreach my $allele (@f) {
	next if $saw{$allele};
	$saw{$allele} = 1;
	push @all_alleles, $allele;
      }
    }
  } else {
    die $alleles_raw;
  }

  my $ref_base = $row->{RefBaseNCBI37} || die;

  my ($reference_sequence, $reference_position);
  if (exists $row->{Position}) {
    # tabix version with standardized position field for target genome
    $reference_sequence = $row->{Chromosome} || die;
    $reference_position = $row->{Position} || die;
  } else {
    # old format
    my $position = $row->{"base(NCBI.37)"} || die;
    my @t = split /\:/, $position;
    die $line unless @t == 2;
    ($reference_sequence, $reference_position) = @t;
  }

  $row->{reference_sequence} = $reference_sequence;
  $row->{reference_position} = $reference_position;
  $row->{reference_base} = $ref_base;

  foreach my $nt (@all_alleles) {
    dump_die($row, "unknown allele $nt") unless length($nt) == 1 and $nt =~ /^[acgt]$/i;
    # newer version contains indels, dinucs: how to handle?
  }

  my $found_ref;
  my @variant_alleles;
  foreach my $allele (@all_alleles) {
    if ($allele eq $ref_base) {
      $found_ref = 1;
    } else {
      push @variant_alleles, $allele;
    }
  }
  die "can't find ref base in allele list" unless $found_ref;
  # sanity check for RC, etc.
  if ($original_format) {
    # in earlier format, variant alleles appear first, then reference allele
    die "ref base is not last allele" unless $ref_base eq $all_alleles[$#all_alleles];
  }
  $row->{variant_bases} = [ @variant_alleles ];

  # if possible, track duplicates by AA change rather than
  # nucleotide position, since we're most interested in the net
  # effect which may be caused by different variants.
  $row->{genes} = [ split /,/, $row->{Genes} ];

  #
  # get total count of patients showing variant allele
  #  - any population
  #  - either homozygous or heterozygous
  #
  my @all_g = split /\//, $row->{AllGenotypeCount};
  die unless @all_g;
  my $g_homozygous_ref = $ref_base x 2;
  die $g_homozygous_ref unless $g_homozygous_ref =~ /^[A-Z][A-Z]$/;
  my $total_patients = 0;
  foreach my $gs (@all_g) {
    my ($genotype, $count) = split /\=/, $gs;
    next if $genotype eq $g_homozygous_ref;
#    print "$ref_base $genotype $count\n";
    $total_patients += $count;
  }
  $row->{total_variant_genotype_patient_count} = $total_patients;
  # count of patients showing some version of the variant genotype

  my $aa;
  if ($original_format) {
    $aa = $row->{AminoAcidChange} || die;
  } else {
    # not present in later version
    $aa = "none";
  }
  if ($aa ne "none") {
    # usable AA
    my @things = split /,/, $aa;
    foreach (@things) {
      if (length($_) == 3) {
	# ok
      } elsif ($_ eq "stop") {
	# reformat to more standard notation
	$_ = "*";
      } else {
	die;
      }
    }
    my $ref_aa = pop @things;
    # reference appears last,
    # e.g. TP53 silver P72R appears as "ARG,PRO" rather
    # than "PRO,ARG"
    # there may be triplets, e.g. LYS,GLN,GLU
    my @var_aa = @things;

    @things = split /\//, $row->{ProteinPos} || die;
    die unless @things == 2;
    my $aa_num = $things[0];
    # 2nd entry is total number of codons

    my @formatted;
    foreach my $var_aa (@var_aa) {
      push @formatted, sprintf '%s%d%s',
      ucfirst(lc($ref_aa)), $aa_num, ucfirst(lc($var_aa));
    }
    $row->{aa_list} = [ @formatted ];
  }
  return $row;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
