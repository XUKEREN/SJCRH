package Chr2GeneParser;
# parse chr2gene files, e.g.
# /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/chr2gene//chr1.gene.txt
# Michael Edmonson 11/2013

use strict;
use FileHandle;

use Configurable;

@Chr2GeneParser::ISA = qw(Configurable Exporter);
@Chr2GeneParser::EXPORT_OK = qw();

use MethodMaker qw(
	all_rows
wanted_genes
wanted_isoforms
wanted_features

bucket_by_isoform
all_isoforms
saw_genes
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->all_rows([]);
  $self->all_isoforms({});
  $self->saw_genes({});
  $self->configure(%options);
  return $self;
}

sub restrict_genes {
  my ($self, $list) = @_;
  $self->wanted_genes({ map {$_, 1} @{$list} });
}

sub restrict_isoforms {
  my ($self, $list) = @_;
  $self->wanted_isoforms({ map {$_, 1} @{$list} });
}

sub restrict_feature_types {
  my ($self, $list) = @_;
  $self->wanted_features({ map {$_, 1} @{$list} });
}

sub parse {
  my ($self, %options) = @_;
  my $file = $options{"-file"} || die;
  my $wanted_genes = $self->wanted_genes();
  my $wanted_isoforms = $self->wanted_isoforms();
  my $wanted_features = $self->wanted_features();
  my $saw_genes = $self->saw_genes();

  my $fh = new FileHandle();
  $fh->open($file) || die;
  my @rows;
  my $VERBOSE = 0;
  $file =~ /(chr\w+)\./ || die "can't find chrom in $file";
  my $chrom = $1;
  while (<$fh>) {
    chomp;
    my @f = split /\t/, $_;
    die unless @f == 4;
    # AADACL3|NM_001103169|+|UTR3  UTR  12785964  12788726
    my ($tags, $feature_type, $start, $end) = @f;
    my ($gene, $accession, $strand, $feature_name) = split /\|/, $tags;

    next if $wanted_genes and not($wanted_genes->{$gene});
    next if $wanted_isoforms and not($wanted_isoforms->{$accession});
    next if $wanted_features and not($wanted_features->{$feature_type});
    # user restrictions

    $saw_genes->{$gene} = 1;
#    printf STDERR "debug gene %s\n", $gene;

    print STDERR "$_\n" if $VERBOSE;

    my %r;
    $r{feature_type} = $feature_type;
    $r{start_interbase} = $start;
    $r{start_inbase} = $start + 1;
    $r{end} = $end;
    $r{gene} = $gene;
    $r{accession} = $accession;
    $r{strand} = $strand;
    $r{feature_name} = $feature_name;
    die if ref $chrom;
    $r{chrom} = $chrom;

    if ($feature_name =~ /exon(\d+)/) {
      my $exno = $1;
#      printf STDERR "Exon %s %s\n", $feature_name, $exno;
      $r{exon_number} = $exno;
    }
    # TO DO: MOAR

    push @{$self->all_rows()}, \%r;
    push @rows, \%r;
    push @{$self->all_isoforms()->{$r{accession}}}, \%r
	if $self->bucket_by_isoform();
  }

  return \@rows;
}

sub qc_wanted_genes {
  my ($self, %options) = @_;
  my $ignore = $options{"-ignore"};

  my $wanted_genes = $self->wanted_genes();
  if ($wanted_genes) {
    my $saw_genes = $self->saw_genes();
    my @fail;
    foreach my $gene (sort keys %{$wanted_genes}) {
      unless ($saw_genes->{$gene}) {
	if ($ignore and $ignore->{$gene}) {
#	  die "hey now";
	} else {
	  push @fail, $gene;
	}
      }
    }
    printf STDERR "ERROR: %d genes not found!: %s\n", scalar(@fail), join ", ", @fail if @fail;
  }
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
