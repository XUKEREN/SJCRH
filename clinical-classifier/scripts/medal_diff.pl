#!/bin/env perl
use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list get_hash_option);
use DelimitedFile;
use Reporter;

my %FLAGS;

my @IGNORE_FIELDS;
my @REPORT_FIELD_BEFORE;

my @clopts = (
	      "-before=s",
	      "-after=s",

	      "-ignore-field=s" => \@IGNORE_FIELDS,
	      "-report-field-before=s" => \@REPORT_FIELD_BEFORE,

	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $f_before = $FLAGS{before} || die "-before";
my $f_after = $FLAGS{after} || die "-after";
$f_after = basename($f_before) if $f_after eq ".";

my $m_before = parse_medal_info($f_before);
my $m_after = parse_medal_info($f_after);

my %f_ignore = map {$_, 1} @IGNORE_FIELDS;
foreach my $v (sort keys %{$m_before}) {
  my $bef = $m_before->{$v} || die;
  my $aft = $m_after->{$v} || die;

  my %all = (%{$bef}, %{$aft});

  foreach my $k (sort keys %all) {

    next if $f_ignore{$k};

    my $v_before = $bef->{$k};
    my $v_after = $aft->{$k};
    foreach ($v_before, $v_after) {
      $_ = "" unless defined $_;
    }

    unless ($v_before eq $v_after) {
      if ($k eq "Reason" or $k eq "Somatic_Reason") {
	my %before = map {$_, 1} split /;/, $v_before;
	my %after = map {$_, 1} split /;/, $v_after;

	my %rall = (%before, %after);

	foreach my $rk (sort keys %rall) {
	  my $rv_before = exists $before{$rk} ? $rk : "";
	  my $rv_after = exists $after{$rk} ? $rk : "";
	  foreach ($rv_before, $rv_after) {
	    $_ = "" unless defined $_;
	  }
	  unless ($rv_before eq $rv_after) {
	    printf "%s\t", join "\t", map {$bef->{$_}} @REPORT_FIELD_BEFORE if @REPORT_FIELD_BEFORE;
	    printf "%s\t%s\t<%s\t>%s\n", $v, $k, $rv_before, $rv_after;
	  }
	}
      } else {
	printf "%s\t%s\t<%s\t>%s\n", $v, $k, $v_before, $v_after;
      }
    }
  }

}


sub parse_medal_info {
  my ($file) = @_;
  my $df = new DelimitedFile("-file" => $file,
			     "-headers" => 1,
			     );

  my @f_pos = qw(Chr WU_HG19_Pos ReferenceAllele MutantAllele);
  my %f_pos = map {$_, 1} @f_pos;

  my %track;
  while (my $row = $df->get_hash()) {
    my $key = join ".", @{$row}{@f_pos};
    my %info;
    foreach my $k (keys %{$row}) {
      if ($f_pos{$k}) {
	# ignore
      } else {
	$info{$k} = $row->{$k};
      }
    }
    die "dup" if $track{$key};
    $track{$key} = \%info;
  }
  return \%track;
}



