package PromoterAnnotation;
# promoter-related annotation

use strict;
use Exporter;

use Configurable;
use ConfigUtils qw(config_or_manual);
use ChrBucketMap;
use MiscUtils qw(dump_die);
use GenomeUtils qw(cook_chromosome_name);
use VariantOverlap;

@PromoterAnnotation::ISA = qw(Configurable Exporter);
@PromoterAnnotation::EXPORT_OK = qw();

use MethodMaker qw(
	genome
f_regions
f_positions
gene_list

cbm_regions
vo_known_sites
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
  my ($self, %options) = @_;

  my $f_positions = config_or_manual(
				   "-config-type" => "genome",
				   "-config-name" => $self->genome,
				   "-parameter" => "CLINCLS_PROMOTER_SITES",
				   "-manual" => $self->f_positions,
      );

  
  my $vo = new VariantOverlap();
  my $df = new DelimitedFile(
			     "-file" => $f_positions,
			     "-headers" => 1,
			    );
  while (my $row = $df->get_hash()) {
    my $chr = $row->{Chr} || die;
    my $pos = $row->{Pos} || die;
    $vo->add_site("-reference" => $chr, "-position" => $pos);
  }
  $self->vo_known_sites($vo);

  my $f_regions = config_or_manual(
				   "-config-type" => "genome",
				   "-config-name" => $self->genome,
				   "-parameter" => "CLINCLS_PROMOTER_REGIONS",
				   "-manual" => $self->f_regions,
      );

  my $gene_list = config_or_manual(
				   "-config-type" => "genome",
				   "-config-name" => $self->genome,
				   "-parameter" => "CLINCLS_PROMOTER_GENES",
				   "-manual" => $self->gene_list,
      );

  my %restrict_genes;
  my %restrict_needed;
  if ($gene_list) {
    %restrict_genes = map {$_, 1} split /,/, $gene_list;
    %restrict_needed = %restrict_genes;
  }

  $df = new DelimitedFile(
			  "-file" => $f_regions,
			  "-headers" => 1,
			 );
  my $pr = new ChrBucketMap(
			 "-f_chr" => "chrom",
			 "-f_start" => "txStart",
			 "-f_end" => "txEnd",
			);
  while (my $row = $df->get_hash()) {
    if (%restrict_genes) {
      my $gene = $row->{name2} || dump_die($row, "no name2");
      if ($restrict_genes{$gene}) {
	delete $restrict_needed{$gene};
      } else {
	next;
      }
    }

    $row->{txStart}++;
    # convert from UCSC interbase to in-base
    $pr->add_row("-row" => $row);
  }

  if (%restrict_needed) {
    die sprintf "gene restriction in effect but didn't see entries for %s", join ",", sort keys %restrict_needed;
  }

  $self->cbm_regions($pr);



}

sub in_promoter_region {
  my ($self, %options) = @_;
  return $self->cbm_regions->find(%options);
}

sub overlaps_known_promoter {
  my ($self, %options) = @_;
  return $self->vo_known_sites->overlaps(%options);
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
