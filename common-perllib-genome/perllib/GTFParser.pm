package GTFParser;
# GTF parser; no doubt a bioperl version but simple enough & wanted
# a callback feature for huge files
# MNE 2/2015

use strict;
use Configurable;
use Exporter;
use FileHandle;

use MiscUtils qw(unquote);

@GTFParser::ISA = qw(Configurable Exporter);
@GTFParser::EXPORT_OK = qw();

my @GTF_FIELDS = qw(
		   seqname
		   source
		   feature
		   start
		   end
		   score
		   strand
		   frame
		   group
		);


use MethodMaker qw(
	
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  return $self;
}

sub parse {
  my ($self, %options) = @_;
  my $file = $options{"-file"} || die "-file";
  my $callback = $options{"-callback"} || die "-callback";
  my $limit = $options{"-limit"};
  my $batch = $options{"-batch"};
  my $fh = new FileHandle();
  $fh->open($file) || die;

  my $lines = 0;
  printf STDERR "parsing %s...", $file;
  my @batch;
  while (<$fh>) {
    next if /^#/;
    chomp;
    my @f = split /\t/, $_;
    die unless @f == @GTF_FIELDS;
    my %r;
    @r{@GTF_FIELDS} = @f;

    if (++$lines and $lines % 200000 == 0) {
      printf STDERR "%d%%...", tell($fh) * 100 / -s $file;
    }
    last if $limit and $lines > $limit;

    my %attr = map {split /\s+/, $_} split /;\s*/, $r{group};
    foreach (values %attr) {
      $_ = unquote($_);
    }
    $r{attributes} = \%attr;

    if ($batch) {
      push @batch, \%r;
      if ($lines % $batch == 0) {
	&$callback(\@batch);
	@batch = ();
      }
    } else {
      &$callback(\%r);
    }
  }
  &$callback(\@batch);
  print STDERR "done\n";
  
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
