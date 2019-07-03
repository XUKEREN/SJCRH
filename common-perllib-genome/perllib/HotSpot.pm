package HotSpot;
# hotspot builder/tracker

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);

@HotSpot::ISA = qw(Configurable Exporter);
@HotSpot::EXPORT_OK = qw();

use MethodMaker qw(
	warn
cooked2raw
track_codon
min_hotspot_count
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->warn({});
  $self->cooked2raw({});
  $self->track_codon({});
  $self->configure(%options);
  return $self;
}

sub add_gene_aa {
  my ($self, %options) = @_;
  my $gene = get_hash_option(\%options, "-gene");
  my $aa = get_hash_option(\%options, "-aa");
  my $track_codon = $self->track_codon();
  my $warn = $self->warn();
  my $cooked2raw = $self->cooked2raw();
  $aa =~ s/^p\.//;
  my $key;
  if ($aa =~ /^([A-Z\*]\d+)/) {
    $key = $1;
  } elsif ($aa =~ /[A-Z\*]+ins([A-Z\*]+\d+)$/) {
    # SJ: RinsT377, 377 = T
    $key = $1;
  } elsif ($aa =~ /^([A-Z\*]+)(\d+)delins[A-Z\*]+$/) {
    # SJ: IPR310delinsTTYML, 310 = I
    $key = substr($1, 0, 1) . $2;
  }

  if ($key) {
    $cooked2raw->{$gene}{$key}{$aa}++;
    $track_codon->{$gene}{$key}++;
    # any need to track raw row data??
  } else {
    printf STDERR "reject %s\n", $aa unless $warn->{$aa};
#    printf STDERR "reject %s\n", $aa;
    $warn->{$aa} = 1;
  }
}

sub find_compound {
  # debug
  my ($self) = @_;
  my $min = $self->min_hotspot_count() || die;
  my $cooked2raw = $self->cooked2raw();
  foreach my $gene (keys %{$cooked2raw}) {
    foreach my $key (keys %{$cooked2raw->{$gene}}) {
      my @aas = keys %{$cooked2raw->{$gene}{$key}};
      my $total = 0;
      foreach my $aa (@aas) {
	$total += $cooked2raw->{$gene}{$key}{$aa};
      }
      printf STDERR "%s\n", join ",", $total, $gene, $key, @aas if @aas > 1 and $total >= $min;
    }
  }
}

sub get_gene_hotspots {
  my ($self) = @_;
  my $min = $self->min_hotspot_count() || die;
  my $track_codon = $self->track_codon || die;
  my @results;
  foreach my $gene (keys %{$track_codon}) {
    foreach my $key (keys %{$track_codon->{$gene}}) {
      my $count = $track_codon->{$gene}{$key};
      if ($count >= $min) {
 	push @results, [$gene, $key];
      }
    }
  }
  return \@results;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
