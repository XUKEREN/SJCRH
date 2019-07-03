package BAMFlagStats;
# parse samtools flagstat files

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);

@BAMFlagStats::ISA = qw(Configurable Exporter);
@BAMFlagStats::EXPORT_OK = qw();

use MethodMaker qw(
	fields
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
  my $file = get_hash_option(\%options, "-file");
  open(FSTMP, $file) || die;
  my %info;

#die $file;
  my @ordered;
  while (<FSTMP>) {
    if (/^(\d+) in (total)$/) {
      $info{$2} = $1;
      push @ordered, $2;
    } elsif (/(\d+) (QC failure)$/) {
      my ($count, $k) = ($1, $2);
      $k =~ s/\s+/_/g;
      $info{$k} = $count;
      push @ordered, $k;
#      dump_die(\%info);
    } elsif (/^(\d+) (duplicates|read[12])$/) {
      my ($count, $k) = ($1, $2);
      $info{$k} = $count;
      push @ordered, $k;
      if ($k eq "duplicates") {
	my $k2 = $k . "_percent";
	my $pct = sprintf '%.2f', $count * 100 / $info{total};
	$info{$k2} = $pct;
	push @ordered, $k2;
      }
    } elsif (/^(\d+) (mapped|properly paired|singletons) \((\S+)%\)$/) {
      my ($count, $k, $pct) = ($1, $2, $3);
      $k =~ s/\s+/_/g;
      $info{$k} = $count;
      $info{$k . "_percent"} = $pct;
      push @ordered, $k;
      push @ordered, $k . "_percent";
    } elsif (/^(\d+) (paired in sequencing)$/) {
      my ($count, $k) = ($1, $2);
      $k =~ s/\s+/_/g;
      $info{$k} = $count;
      push @ordered, $k;
    } elsif (/^(\d+) with (itself and mate mapped|mate mapped to a different chr|mate mapped to a different chr \(mapQ>=5\))$/) {
      my ($count, $k) = ($1, $2);
      $k =~ s/\s+/_/g;
      $info{$k} = $count;
      push @ordered, $k;
    } elsif (/^(\d+) \+ (\d+) in total \(QC-passed reads \+ QC-failed reads\)$/) {
      # 14609280 + 0 in total (QC-passed reads + QC-failed reads)
      $info{total} = $1;
      $info{"QC failure"} = $2;
      push @ordered, "total";
      push @ordered, "QC failure";
    } elsif (/^(\d+) \+ (\d+) (duplicates|read[12]|paired in sequencing|with itself and mate mapped|with mate mapped to a different chr|with mate mapped to a different chr \(mapQ>=5\))$/) {
      my ($count, $count2, $k) = ($1, $2, $3);
      $k =~ s/^with //;
      # backwards compatibility
      $k =~ s/\s+/_/g;
      $info{$k} = $count;
      push @ordered, $k;
    } elsif (/^(\d+) \+ (\d+) (mapped|properly paired|singletons) \((\d+\.\d+)%/) {
      my ($count, $count2, $k, $pct) = ($1, $2, $3, $4);
      $k =~ s/\s+/_/g;
      $info{$k} = $count;
      $info{$k . "_percent"} = $pct;
      push @ordered, $k;
      push @ordered, $k . "_percent";
    } else {
      dump_die(\%info, "unhandled line: $_");
    }
  }
  $self->fields(\@ordered);

  return \%info;

}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
