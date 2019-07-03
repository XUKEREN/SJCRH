package TARGETBarcode;
# parse TARGET sample barcode

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);

@TARGETBarcode::ISA = qw(Configurable Exporter);
@TARGETBarcode::EXPORT_OK = qw();

use MethodMaker qw(
                   full_barcode

                   tumor_code
                   sample
                   tissue_code
                   aliquot
                   nucleic_code
		  );

my %TARGET2SJ = (
  "01" => "D",
  # Primary Tumor: Primary Solid Tumor
  "02" => "R",
  # Recurrent Tumor: Recurrent Solid Tumor
  "03" => "D",
  # Primary Blood Cancer: Primary Blood Derived Cancer – Peripheral blood
  "04" => "R",
  # Recurrent Blood Cancer: Recurrent Blood Derived Cancer - Bone Marrow
  "05" => "D",
  # Additional - New Primary: Additional  - New Primary
  "06" => "M",
  # Metastatic: Metastatic
  "07" => "M",
  # Additional Metastatic: Additional Metastatic

#    Post neo-adjuvant therapyTissue disease-specific post-adjuvant therapy08
  # ?

  "09" => "D",
  # Primary Blood Cancer BM: Primary Blood Derived Cancer – Bone Marrow

  "10" => "G",
  # Blood Derived Normal: Blood Derived Normal

  "11" => "G",
  # Solid Tissue Normal: Solid Tissue Normal

  "12" => "G",
  # Buccal Cell Normal: Buccal Cell Normal

#    EBV NormalEBV Immortalized Normal13
# ?cell line?

  "14" => "G",
  # BM Normal: Bone Marrow Normal

  "15" => "G",
  # Fibroblast Normal: Fibroblasts from Bone Marrow Normal
  
  "16" => "G",
  # Mononuclear Cell Normal: Mononuclear Cells from Bone Marrow Normal

#  "20" => "C",
  # Cell Line Control: Cell Line Control (Control Analyte)
  # does SJ "C" assume disease or normal??

  "40" => "R",
  # Recurrent Blood Cancer: Recurrent Blood Derived Cancer – Peripheral blood

#    Post treatment Blood Cancer Bone Marrow Blood Derived Cancer- Bone Marrow, Post-treatment41
#    Post treatment Blood Cancer BloodBlood Derived Cancer- Peripheral Blood, Post-treatment42
  # Cancer cell lineCell line from patient tumor50

  "60" => "X",
  # Xenograft, primary: Xenograft from patient not grown as intermediate on plastic tissue culture dish
  "61" => "X",
  # Xenograft, cell-line derived: Xenograft grown in mice from established cell lines61
  # GranulocytesGranulocytes after a Ficoll separation99
);
# see "OCG Sample Codes 060815 final.docx"

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  return $self;
}

sub parse {
  my ($self, %options) = @_;
  my $allow_sample_only = $options{"-allow-sample-only"};
  my $barcode = get_hash_option(\%options, "-barcode");
  my $parse_ok;
  my ($tumor_code, $sample, $tissue_code, $aliquot, $nucleic_code);
  foreach ($tumor_code, $sample, $tissue_code, $aliquot, $nucleic_code) {
    $_ = "";
  }
  my $full_barcode = "";

  if ($barcode =~ /TARGET\-(\d+)-([A-Z]{6})\-(\d\d)([A-Z])\-(\d\d[A-Z]?)/) {
    # full barcode, e.g. TARGET-10-PARTRW-04A-01D
    ($tumor_code, $sample, $tissue_code, $aliquot, $nucleic_code) = 
	($1, $2, $3, $4, $5);
    $full_barcode = 1;
  } elsif ($barcode =~ /TARGET\-(\d+)-([A-Z]{6})\-(\d\d)([A-Z])/) {
    # limited barcode
    # e.g. SJBALL020579_G1-TARGET-10-PANZPJ-14A.bam
    ($tumor_code, $sample, $tissue_code, $aliquot) = 
	($1, $2, $3, $4);
  } elsif ($allow_sample_only and $barcode =~ /[\-_]([A-Z]{6})[\-_\.]/) {
    # just the TARGET sample only, e.g.
    # /nfs_exports/genomes/1/projects/VALCAP/TARGET_ALL/BucketRaw/SJBALL/SJBALL021159_R1-PAPLDL_R.bam
    # - /nfs_exports/genomes/1/projects/VALCAP/TARGET_ALL/BucketRaw/SJCOGALL/SJCOGALL010221_D1-PAPAGK.bam
    $sample = $1;
  }
  # TO DO: alternate format
  $parse_ok = 1 if $sample;

  $self->tumor_code($tumor_code);
  $self->sample($sample);
  $self->tissue_code($tissue_code);
  $self->aliquot($aliquot);
  $self->nucleic_code($nucleic_code);
  $self->full_barcode($full_barcode);

  return $parse_ok;
}

sub get_sj_disease_code {
  my ($self, %options) = @_;
  my $tissue_code = $self->tissue_code;
  return $TARGET2SJ{$tissue_code} || die "can't map $tissue_code";
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
