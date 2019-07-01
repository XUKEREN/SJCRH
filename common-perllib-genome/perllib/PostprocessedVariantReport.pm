package PostprocessedVariantReport;
# standard SJ variant postprocessing report: deal with header variants
# MNE 11/2012

use strict;
use Carp qw(confess);

use Exporter;

use constant F_GENE => "GeneName";
use constant F_REFSEQ => "mRNA_acc";
use constant F_AACHANGE => "AAChange";

use Configurable;

use DelimitedFile;

@PostprocessedVariantReport::ISA = qw(Configurable Exporter);
@PostprocessedVariantReport::EXPORT_OK = qw(
F_GENE
F_REFSEQ
F_AACHANGE
);

use MethodMaker qw(
	rows
	ignore_blank_rows
standardize_chr
skip_indels
skip_format_check
headers
auto_repair_swapped_alleles

headers_raw
		  );

my %SJ_CALL_CODES = map {$_, 1} qw(
Good
2SNP
SNPCluster
GoodClone
2SNPClone
GoodSNP
SNPClusterClone
);

use constant HEADERS_16 => (
			    # standard headers
			    F_GENE,
			    "SJQuality",
			    "Sample",
			    "Chr",
			    "WU_HG19_Pos",
			    "Class",
			    F_AACHANGE,
			    "ProteinGI",
			    F_REFSEQ,
			    "#Mutant_In_Tumor",
			    "#Total_In_Tumor",
			    "#Mutant_In_Normal",
			    "#Total_In_Normal",
			    "ReferenceAllele",
			    "MutantAllele",
			    "Flanking"
			   );

use constant HEADERS_17_MAF => (
			    F_GENE,
			    "SJQuality",
			    "Sample",
			    "Chr",
			    "WU_HG19_Pos",
			    "Class",
			    F_AACHANGE,
			    "ProteinGI",
			    F_REFSEQ,
			    "#Mutant_In_Tumor",
			    "#Total_In_Tumor",
			    "#Mutant_In_Normal",
			    "#Total_In_Normal",

			    "MAF",
			    # some downstream reports insert
			    # mutant allele frequency

			    "ReferenceAllele",
			    "MutantAllele",
			    "Flanking"
			   );


use constant HEADERS_17_WITH_CALL => (
			    F_GENE,
			    "SJQuality",
			    "Sample",
			    "Chr",
			    "WU_HG19_Pos",
			    "Class",
			    F_AACHANGE,
			    "ProteinGI",
			    F_REFSEQ,
			    "#Mutant_In_Tumor",
			    "#Total_In_Tumor",
			    "#Mutant_In_Normal",
			    "#Total_In_Normal",
			    "ReferenceAllele",
			    "MutantAllele",
			    "Flanking",
                            "SJCall",
                            # not sure of column name: a judgment call
			   );

use constant HEADERS_13 => (
			    F_GENE,
			    "SJQuality",
			    "Sample",
			    "Chr",
			    "WU_HG19_Pos",
			    "Class",
			    F_AACHANGE,
#			    "ProteinGI",
			    F_REFSEQ,
			    "#Mutant_In_Tumor",
			    "#Total_In_Tumor",
#			    "#Mutant_In_Normal",
#			    "#Total_In_Normal",
			    "ReferenceAllele",
			    "MutantAllele",
			    "Flanking"  # only a portion
			   );
# 3/2013: yet another flavor for 454-style

use constant HEADERS_9 => (
			    F_GENE,
			    "SJQuality",
			    "Sample",
			    "Chr",
			    "WU_HG19_Pos",
			    "Class",
#			    F_AACHANGE,
#			    "ProteinGI",
			    F_REFSEQ,
#			    "#Mutant_In_Tumor",
#			    "#Total_In_Tumor",
#			    "#Mutant_In_Normal",
#			    "#Total_In_Normal",
			    "ReferenceAllele",
			    "MutantAllele",
#			    "Flanking"  # only a portion
			   );
# 12/2013: yet another flavor for 454-style
# /nfs_exports/genomes/1/PCGP/BucketIntermediate/bvadodaria/MiSeq/COGALL_validation/*targets
# NRAS    SJHQ    SJCOGGALL010223 1       115258744       missense        NM_002524       C       T

use constant HEADERS_9_MK2 => (
			    F_GENE,
			    "SJQuality",
			    "Sample",
			    "Chr",
			    "WU_HG19_Pos",
			    "Class",
			    "ReferenceAllele",
			    "MutantAllele",
			    F_REFSEQ,
			   );
# 2/2014:
# KILL ME NOW: A DIFFERENT 9-column format
# /nfs_exports/genomes/1/PCGP/BucketIntermediate/bvadodaria/MiSeq/Clingen_uncaptured_targets.txt/Clingen_uncapped_bin1.targets
# NBPF10  SJHQ    SJE2A007D       chr1    145322843       MISSENSE        C       T       NM_001039703


sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->rows([]);
  $self->ignore_blank_rows(1);
  $self->standardize_chr(1);
  $self->auto_repair_swapped_alleles(1);
  $self->skip_format_check({});
  $self->configure(%options);
  return $self;
}

sub parse {
  my ($self, %options) = @_;
  my $file = $options{"-file"} || die "-file";
  my $header = $options{"-header"};
#  $header = 1 unless defined $header;
  unless (defined $header) {
    open(HCHK, $file) || confess "can't open $file: $!";
    my $first = <HCHK>;
    close HCHK;
    $header = $first =~ /SJQuality/ ? 1 : 0;
  }
  my $rows = $self->rows();
  if ($header) {
    # file already has headers
    my $df = new DelimitedFile(
			       "-file" => $file,
			       "-headers" => 1,
			     );
    $self->headers_raw($df->headers_raw);
    # existing user-supplied headers, which may add additional columns
    # that we don't know about or handle

    while (my $row = $df->get_hash()) {
#      die join "\n", sort keys %{$row};
      foreach my $h (HEADERS_16) {
	unless (exists $row->{$h}) {
	  if ($h eq "WU_HG19_Pos") {
	    die "need either WU_HG19_Pos / WU_HG18_Pos / Pos" unless exists $row->{WU_HG18_Pos} or $row->{Pos};
	  } elsif ($h eq "#Mutant_In_Tumor") {
	    tolerate_check($row, "#Mutant_In_Tumor", "Mutant_In_Tumor");
	  } elsif ($h eq "#Total_In_Tumor") {
	    tolerate_check($row, "#Total_In_Tumor", "Total_In_Tumor");
	  } elsif ($h eq "#Mutant_In_Normal") {
	    tolerate_check($row, "#Mutant_In_Normal", "Mutant_In_Normal");
	  } elsif ($h eq "#Total_In_Normal") {
	    tolerate_check($row, "#Total_In_Normal", "Total_In_Normal");
	  } elsif (
	    $h eq F_GENE or 
	    $h eq F_AACHANGE or
	    $h eq "ProteinGI" or
	    $h eq F_REFSEQ 
	      ) {
	    # tier 1 has these fields, tier 2+ does not
	    tolerate_check($row, $h, "Class");
	  } else {
	    confess "missing header $h in $file";
	  }
	}
      }
#      $self->format_check($row);

      if ($self->ignore_blank_rows()) {
	my @usable = grep {$_} values %{$row};
	next unless @usable;
      }

      $self->repair_swapped_alleles($row) if $self->auto_repair_swapped_alleles();

      push @{$rows}, $row;
    }
    $self->headers([ HEADERS_16 ]);
  } else {
    my $df = new DelimitedFile(
			       "-file" => $file,
			       "-headers" => 0
			     );
    my @headers;
    while (my $row = $df->next("-ref" => 1)) {
      if ($self->ignore_blank_rows) {
	my @usable = grep {$_} @{$row};
	next unless @usable;
      }

      my %r;

      unless (@headers) {
	# detect column formatting
	if (@{$row} == 16 or $options{"-force-16"}) {
	  @headers = HEADERS_16;
	} elsif ($options{"-force-17"}) {
	  @headers = HEADERS_17_WITH_CALL;
	} elsif (@{$row} == 17) {
	  my %test;
	  @test{HEADERS_17_MAF()} = @{$row};
#	  my $has_call = ($test{SJCall} and $SJ_CALL_CODES{$test{SJCall}});
	  # 9/2013: more variety in call codes, and this seems to be the
	  # more common format, so reverse logic
	  my $has_call = 1;
	  if ($test{MAF} =~ /\d/) {
	    $has_call = 0;
	    die "TEST ME: possible MAF format";
	  }
	  @headers = $has_call ? HEADERS_17_WITH_CALL : HEADERS_17_MAF;
	} elsif (@{$row} == 13) {
	  @headers = HEADERS_13;
	  $self->skip_format_check()->{"ProteinGI"} = 1;
	} elsif (@{$row} == 9) {
	  $self->skip_format_check()->{"Flanking"} = 1;
	  $self->skip_format_check()->{"ProteinGI"} = 1;

	  my %test;
	  @test{HEADERS_9()} = @{$row};
	  my $mrna_acc = $test{mRNA_acc};
	  if ($mrna_acc =~ /_/) {
	    # version 1
	    @headers = HEADERS_9;
	  } elsif (length($mrna_acc) == 1 and $mrna_acc =~ /^[ACGT]$/) {
	    # version 2
	    @headers = HEADERS_9_MK2;
	  }
	  %test = ();
	  @test{@headers} = @{$row};
	  die "test failed" unless $test{mRNA_acc} =~ /_/;
	} else {
	  confess "unknown postprocessed variant report column count in $file: " . scalar @{$row};
	}
	$self->headers(\@headers);
      }

      @r{@headers} = @{$row};
      $self->format_check(\%r);
      push @{$rows}, \%r;
    }
  }

  if ($self->skip_indels()) {
    my @filtered;
    foreach my $row (@{$rows}) {
      next if $row->{Flanking} =~ /\(/;
      push @filtered, $row;
    }
    $self->rows(\@filtered);
  }

  return $self->rows();

}

sub format_check {
  my ($self, $row) = @_;
  my $broken_field;
  $broken_field = "SJQuality" unless $row->{SJQuality} =~ /^SJ/;
#  $broken_field = F_GENE unless $row->{F_GENE()} =~ /^[\w\-]+$/;
  $broken_field = F_GENE unless $row->{F_GENE()} =~ /^[\w\-\.]+$/;
  # 3/2014: allow HGC6.3
  $broken_field = "Sample" unless $row->{Sample} =~ /^SJ/;
  unless ($row->{Chr} =~ /^chr/) {
    if ($self->standardize_chr()) {
      $row->{Chr} = "chr" . $row->{Chr};
    } else {
      $broken_field = "Chr";
    }
  }

  $broken_field = "WU_HG19_Pos" unless $row->{WU_HG19_Pos} =~ /^\d+$/;

  unless ($self->skip_format_check->{ProteinGI}) {
    $broken_field = "ProteinGI" unless $row->{ProteinGI} =~ /^\d+$/ or $row->{ProteinGI} eq "-1";
  }

  unless ($self->skip_format_check->{Flanking}) {
    $broken_field = "Flanking" unless $row->{Flanking} =~ /^[ACGT]/;
  }

  foreach my $f (qw(ReferenceAllele MutantAllele)) {
    my $v = $row->{$f} || die;
    if ($v eq "-") {
      # indel, OK
    } elsif ($v =~ /^[acgt]+$/i) {
      # ok
    } else {
      $broken_field = $f;
    }
  }


  confess sprintf "broken field %s, value %s, row=%s", $broken_field, $row->{$broken_field}, join ",", values %{$row} if $broken_field;
  return 1;
}

sub skip_format_check_for_field {
  my ($self, $field) = @_;
  $self->skip_format_check()->{$field} = 1;
}

sub get_variant_string {
  my ($self, $row) = @_;
  return sprintf '%s.%d.%s.%s', @{$row}{qw(Chr WU_HG19_Pos ReferenceAllele MutantAllele)};
}

sub is_snv {
  my ($self, $row) = @_;
  return $row->{Flanking} =~ /\[/ ? 1 : 0;
}

sub get_ordered_headers {
  my ($self) = @_;
  return $self->headers();
}

sub tolerate_check {
  # STATIC
  my ($row, @allowed) = @_;
  my $ok;
  foreach (@allowed) {
    $ok = 1 if exists $row->{$_};
  }
  die sprintf 'need either %s', join " / ", @allowed unless $ok;
}

sub repair_swapped_alleles {
  my ($self, $row) = @_;
  if ($row->{Flanking} and $row->{Flanking} =~ /\[(\w)\/(\w)\]/) {
    my ($actual_reference, $actual_variant) = ($1, $2);
    my $claimed_reference = $row->{ReferenceAllele} || die;
    my $claimed_variant = $row->{MutantAllele} || die;
    if ($claimed_reference ne $actual_reference or
	$claimed_variant ne $actual_variant) {
      # mismatch
      if ($claimed_reference eq $actual_variant and
	  $claimed_variant eq $actual_reference) {
	# pure swap: repair
#	printf STDERR "WARNING: fixing swapped ReferenceAllele and VariantAllele\n";
	$row->{ReferenceAllele} = $actual_reference;
	$row->{MutantAllele} = $actual_variant;
      } else {
	die "unhandled mismatch";
      }
    }
  }
  die "FAIL" if $row->{ReferenceAllele} eq $row->{MutantAllele};
}

sub get_headers {
  # if the raw report appeared to have compatible headers, return those.
  # otherwise (or if headerless) return the autodetected detected headers.
  my ($self) = @_;
  return $self->headers_raw || $self->headers();
}

1;
