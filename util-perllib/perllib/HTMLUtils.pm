package HTMLUtils;
# HTML utilites

use strict;
use warnings;

use Configurable;
use Exporter;

use HTML::TableExtract;
use HTML::Entities;
use Getopt::Long;

@HTMLUtils::ISA = qw(Configurable Exporter);
@HTMLUtils::EXPORT_OK = qw(
parse_html_tables
);

sub parse_html_tables {
  my (%options) = @_;
  my $blob = $options{"-blob"};
  my $f = $options{"-file"};
  my $dump = $options{"-dump"};
  # report
  my $require_header = $options{"-require-header"};
  my $header_mode = defined $options{"-headers"} ? $options{"-headers"} : 1;

  my @results;

  if ($f) {
    local $/ = undef;
    open(HTMLTMP, $f) || die;
    $blob = <HTMLTMP>;
    close HTMLTMP;
  } elsif ($blob) {
  } else {
    die "specify -file or -blob";
  }

  die "no HTML blob" unless $blob;
  $blob = decode_entities($blob);

  my @options;
  push @options, ("keep_html" => 1) if $options{"-keep-html"};
  my $te = new HTML::TableExtract(@options);
  $te->parse($blob);

  my @tables = $te->table_states();

  my $ti = 0;
  foreach my $ts (@tables) {
    my (@rows) = $ts->rows;
    printf STDERR "NEW TABLE (index %d):\n", $ti++ if $dump;
    my @rows_out;
    my @headers;
    my $table_usable = 1;
    foreach (@rows) {
      my @cells = @{$_};
      foreach (@cells) {
	$_ = "" unless defined $_;
	s/^\s+//;
	s/\s+$//;
#	$_ = "" unless /\w/;
        # LOVD databases: blanks out "Path." field, e.g. "+/+"!
	$_ = "" unless /\S/;
      }
      next unless grep {/\w/} @cells;
      # empty row
      printf STDERR "%s\n", join "\t", @cells if $dump;
      if ($header_mode) {
	if (@headers) {
	  my %r;
	  @r{@headers} = @cells;
	  push @rows_out, \%r;
	} else {
	  @headers = @cells;
	  if ($require_header) {
	    $table_usable = 0 unless grep {$_ eq $require_header} @headers;
	  }
	}
      } else {
	push @rows_out, \@cells;
      }
    }
    print STDERR "\n" if $dump;
    push @results, \@rows_out if $table_usable;
  }

  return \@results;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
