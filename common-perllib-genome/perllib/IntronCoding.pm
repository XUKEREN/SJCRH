package IntronCoding;
# generate synthetic amino acid sequence into intronic sequence
# - useful for in-frame checking (see Cicero sv_inframe.pl)

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(dump_die);
use GenomeUtils qw(reverse_complement);

@IntronCoding::ISA = qw(Configurable Exporter);
@IntronCoding::EXPORT_OK = qw();

use MethodMaker qw(
	r2g
strand
generated
target
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  return $self;
}

sub generate_intron_coding {
  my ($self, %options) = @_;
  my $r2g = $self->r2g() || die;
  my $hit = $options{"-hit"} || die "-hit";
  my $chrom = $hit->{chrom} || die "no chrom in hit";
  my $base_exon_edge = $options{"-edge"} || die "-edge";
  my $base_target = $options{"-target"} || die "-target";
  $self->target($base_target);

  my $map = $hit->{codon_map} || die;
  my $genome = $r2g->fai()->get_sequence("-id" => $chrom);

  my %map;
  my $strand = $map->[1]->{ref_base_num} > $map->[0]->{ref_base_num} ? "+" : "-";
  $self->strand($strand);
  foreach my $r (@{$map}) {
    $map{$r->{ref_base_num} || die} = $r;
  }
  my $genomic_direction = $base_exon_edge < $base_target ? 1 : -1;

  my $edge = $map{$base_exon_edge} || die "$base_exon_edge not in exon";
  my $build_direction;

  if ($strand eq "+") {
    # transcript on + strand
    if ($genomic_direction > 0) {
      # target is downstream, e.g. possible intron retention of geneA
      # ichack -genome GRCh37-lite -refflat NM_003128.txt -nm NM_003128 -chrom 2 -target 54839486 -exon-edge 54839471
      $build_direction = 1;
    } else {
      # target is upstream, e.g. possibly inframe via coding->genomic hit?
      $build_direction = -1;
      # ichack -genome GRCh37-lite -refflat NM_003128.txt -nm NM_003128 -chrom 2 -target 54826200 -exon-edge 54826229
      #   this example is a split exon
    }
  } elsif ($strand eq "-") {
    if ($genomic_direction > 0) {
      # target is upstream
      $build_direction = -1;
      # ichack -genome GRCh37-lite -refflat NM_004304.txt -nm NM_004304 -chrom 2 -target 29447684 -exon-edge 29446394
    } else {
      $build_direction = 1;
      # ichack -genome GRCh37-lite -refflat NM_004304.txt -nm NM_004304 -chrom 2 -target 29446200 -exon-edge 29446208
    }
  } else {
    die "unhandled";
  }

  my $pos = $base_exon_edge + $genomic_direction;
  die "1st neighbor base is not intronic" if $map{$pos};
  my $codon_base = $edge->{codon_base} || die;
  my $codon_number = $edge->{codon_number} || die;
  my %encode;
  my %generated;
  while (!$map{$pos}) {
    #
    # generate artificial codons from the intronic sequence until
    # we reach the next exon
    #
    my %r;
    $r{ref_base_num} = $pos;
    $r{ref_base_genome} = substr($$genome, $pos - 1, 1);
    my $tbase = $r{ref_base_genome};
    if ($strand eq "+") {
    } else {
      $tbase = reverse_complement($tbase);
    }
    $r{ref_base_transcript} = $tbase;

    $codon_base += $build_direction;
    if ($codon_base > 3) {
      $codon_base = 1;
      $codon_number++;
    } elsif ($codon_base < 1) {
      $codon_base = 3;
      $codon_number--;
    }
    $r{codon_base} = $codon_base;
    $r{codon_number} = $codon_number;

    $encode{$codon_number}{$codon_base} = $tbase;

#    printf STDERR "%s\n", join " ", @r{qw(ref_base_num ref_base_genome ref_base_transcript codon_base codon_number)};

    $generated{$pos} = \%r;
    $pos += $genomic_direction;
  }

  #
  #  associate codons with amino acids:
  #
  my $ct = new Bio::Tools::CodonTable();
  my %codon2aa;
  foreach my $codon_number (sort {$a <=> $b} keys %encode) {
    if (scalar(keys %{$encode{$codon_number}}) == 3) {
      my $codon = join "", @{$encode{$codon_number}}{qw(1 2 3)};
      $codon2aa{$codon_number} = $ct->translate($codon);
    }
  }

  foreach my $pos (sort keys %generated) {
    my $ref = $generated{$pos};
    $ref->{AA} = $codon2aa{$ref->{codon_number}} || undef;
    printf STDERR "%s\n", join " ", @{$ref}{qw(ref_base_num codon_number AA ref_base_genome ref_base_transcript codon_base)};
  }

  $self->generated(\%generated);
}

sub get_target_chunk {
  # get a chunk of generated AA sequence in the vicinity of the
  # target site
  my ($self, %options) = @_;
  my $transcript_direction = $options{"-direction"} || die "-direction";
  # this the direction from the TRANSCRIPT's perspective,
  # NOT the genomic perspective
  my $length = $options{"-length"} || die "-length";
  my $skip_first = $options{"-skip-first"};
  die "specify -skip-first [0|1]" unless defined $skip_first;

  my $target = $self->target() || die;
  my $strand = $self->strand() || die;
  my $genomic_direction;
  if ($strand eq "+") {
    if ($transcript_direction == -1) {
      $genomic_direction = -1;
    } elsif ($transcript_direction == 1) {
      $genomic_direction = 1;
    } else {
      die;
    }
  } elsif ($strand eq "-") {
    if ($transcript_direction == -1) {
      $genomic_direction = 1;
    } elsif ($transcript_direction == 1) {
      $genomic_direction = -1;
    } else {
      die;
    }
  } else {
    die;
  }

  my @aa;
  my %saw_codons;
  my $pos = $target;
  my $generated = $self->generated();
  my $first = 1;
  while (@aa < $length) {
    my $r = $generated->{$pos} || last;
    my $cnum = $r->{codon_number};
    unless ($saw_codons{$cnum}) {
      my $usable = 1;
      $usable = 0 if $first and $skip_first;
      my $aa = $r->{AA} || last;
      push @aa, $aa if $usable;
      $saw_codons{$cnum}++;
    }
    $pos += $genomic_direction;
    $first = 0;
  }

  return join "", @aa;
}



1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
