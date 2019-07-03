package VCFBedOverlap;
# wrapper to bed overlap VCF filtering
# MNE 5/2015
# WARNING: this utility may strip INFO tags from output!

use strict;

use File::Basename;

use Configurable;
use Exporter;
use TdtConfig;
use Cluster;
use FileUtils qw(find_binary);

@VCFBedOverlap::ISA = qw(Configurable Exporter);
@VCFBedOverlap::EXPORT_OK = qw();

use constant RAM_REQUIRED => 15000;

use MethodMaker qw(
bed

command
ram
outfile
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
  $self->ram(RAM_REQUIRED);
}

sub overlap {
  my ($self, %options) = @_;
  my $vcf = $options{"-vcf"} || die "-vcf";
  my $bed = $self->bed || die "-bed";

  my %ref_bed;
  my %ref_vcf;

  find_binary("vcftools", "-die" => 1);

  unless ($options{"-no-qc"}) {
    die "fix me: use GenomeUtils::check_reference_name_compatibility";

    printf STDERR "checking VCF/bed compatibility...\n";
    printf STDERR "  scanning %s...", $bed;
    open(BEDTMP, $bed) || die;
    while (<BEDTMP>) {
      next if /^#/;
      /^(\S+)/ || die;
      printf STDERR "%s...", $1 unless $ref_bed{$1};
      $ref_bed{$1} = 1;
    }
    print STDERR "\n";

    open(VCFTMP, $vcf) || die;
    printf STDERR "  scanning %s...", $vcf;
    while (<VCFTMP>) {
      next if /^#/;
      /^(\S+)/ || die;
      printf STDERR "%s...", $1 unless $ref_vcf{$1};
      $ref_vcf{$1} = 1;
    }
    print STDERR "\n";

    foreach my $r (sort keys %ref_bed) {
      my $errors;
      unless ($ref_vcf{$r}) {
	printf STDERR "ERROR: VCF does not contain an entry for bed reference %s\n", $r;
	$errors = 1;
      }
      die if $errors;
    }
  }

  my $outfile = basename($vcf) . ".bed.recode.vcf";

  my $cmd = sprintf 'vcftools --bed %s --vcf %s --out %s --recode',
    $bed,
      $vcf,
	basename($vcf) . ".bed";

  $self->outfile($outfile);
  $self->command($cmd);
  $self->run_cluster() if $options{"-cluster"};

  return $cmd;
}

sub run_cluster {
  my ($self, %options) = @_;

  my $c = new Cluster(
		      "-outfile" => $self->outfile,
		      "-project" => "PCGP",
		     );

  $c->node_class("");
  $c->memory_reserve_mb($self->ram());
  $c->memory_limit_mb($self->ram());
  $c->command($self->command() || die);
  $c->run();
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
