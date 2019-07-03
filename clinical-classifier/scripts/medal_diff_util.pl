#!/bin/env perl
# parse medal diffs to identify differences

use strict;
use warnings;

use 5.10.1;
use Getopt::Long;

use MiscUtils qw(dump_die);
use DelimitedFile;

# use Reporter;

my %FLAGS;
GetOptions(\%FLAGS,
	   "-medals=s",
	   # example medal file to pull headers from
	   "-diff=s",
	   # diff output
	  );

my $medals = $FLAGS{medals} || die;
my $diff =  $FLAGS{diff} || die;

my $df = new DelimitedFile("-file" => $medals,
			   "-headers" => 1,
			  );
my $headers = $df->headers_raw();
open(IN, $diff) || die;
my %track;
while (<IN>) {
  chomp;
  if (/^([<>])\s+(.*)/) {
    my $where = $1;
    my $line = $2;
    die "argh" unless $line =~ /\t/;
    my @f = split /\t/, $line, scalar @{$headers};
    die unless @f == @{$headers};
    my %r;
    @r{@{$headers}} = @f;

    my $key = get_key(\%r);
    die if $track{$key}{$where};
    $track{$key}{$where} = \%r;
  }
}

foreach my $key (sort keys %track) {
  my $before = $track{$key}{"<"};
  my $after = $track{$key}{">"};

  my $b_reason = parse_reason($before);
  my $a_reason = parse_reason($after);

  my %all = map {$_ => 1} (keys %{$b_reason}, keys %{$a_reason});

  foreach my $r (keys %all) {
    if ($b_reason->{$r} and $a_reason->{$r}) {
      # same tag appears in both, remove
      delete $b_reason->{$r};
      delete $a_reason->{$r};
    }
  }

  printf "%s before:%s after:%s\n",
    $key,
      join(",", sort keys %{$b_reason}),
	join(",", sort keys %{$a_reason});

}

sub parse_reason {
  my ($row) = @_;
  my @f = split /;/, $row->{Reason};
  my %r = map {$_ => 1} @f;
  return \%r;
}


sub get_key {
  my ($row) = @_;
  my @f = @{$row}{qw(Chr WU_HG19_Pos ReferenceAllele MutantAllele)};
  foreach (@f) {
    die unless $_;
  }
  return join ".", @f;
}
