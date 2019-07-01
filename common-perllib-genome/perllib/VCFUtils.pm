package VCFUtils;
# utility routines to complement vcftools' Vcf.pm
# MNE 3/2015

use strict;
use Configurable;
use MiscUtils qw(dump_die);

use Exporter;

@VCFUtils::ISA = qw(Configurable Exporter);
@VCFUtils::EXPORT_OK = qw();

use MethodMaker qw(
vcf	
agtags
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  return $self;
}

sub get_agtags {
  my ($self) = @_;
  my $agtags = $self->agtags();
  $agtags = $self->agtags($self->vcf->has_AGtags()) unless $agtags;
  return $agtags;
}

sub get_alt_rows {
  # generate multiple output rows, one for each alternate allele,
  # with appropriately parsed tag info
  my ($self, %options) = @_;
  my $hash = $options{"-hash"} || die "-hash";
  # Vcf.pm next_data_hash();

  my @out_rows;
  my $agtags = $self->get_agtags();
  my $a_fields = $agtags->{infoA};
  # may not be present for all VCFs
  $a_fields = [] unless $a_fields;

  my %split;
  my $info = $hash->{INFO} || die;
#  dump_die($info, "debug", 1);

  my $alt_alleles = $hash->{ALT} || die;
  my $alt_allele_count = scalar @{$alt_alleles};

  foreach my $af (@{$a_fields}) {
    my $raw = $info->{$af};
    if (defined $raw) {
      my @v = split /,/, $raw;
#      printf STDERR "ERROR: %s has %d columns, expected %d (%s)!\n",
#      $af, scalar(@v), $alt_allele_count, $raw unless @v == $alt_allele_count;
      # appears to be out of spec (v 4.1, section 1.2.2)
      # something I'm missing here??
      #
      # ExAC:
      # X.153008476, reference is T, variant is either C or A.
      #
      # Dale: "the two alternates allow for the creation of 3 possible
      #       heterozygote genotypes: TC, TA, CA"
      #
      # so, we don't have to worry about additional records if
      # we only need data for the reference allele vs. one of the
      # alternates.
      
      $split{$af} = \@v;
    }
  }
#  dump_die($info);

  for (my $idx = 0; $idx < $alt_allele_count; $idx++) {
    my $alt = $alt_alleles->[$idx];
#    dump_die($hash, "Debug", 1);
    my %new;

    # global info:
    $new{vcf_hash} = $hash;
    # main hash containing all info
    $new{CHROM} = $hash->{CHROM} || die;
    $new{POS} = $hash->{POS} || die;
    $new{REF} = $hash->{REF} || die;

    # information split out for this ALT allele:
    $new{ALT} = $alt;
    $new{INFO} = {};

    $new{vcf_alt_entry} = $idx + 1;
    # useful for genotype lookup

    foreach my $af (@{$a_fields}) {
      if (defined $split{$af}) {
	$new{INFO}{$af} = $split{$af}->[$idx];
      }
    }
#    dump_die(\%new, "debug2", 1);
#    dump_die($new{INFO}, "debug3", 1);
    push @out_rows, \%new;
  }
  return \@out_rows;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
