package EUtilities;
# "lite" e-utilities implementation
# http://www.ncbi.nlm.nih.gov/books/NBK25499/
# MNE 10/2014

use strict;

use Configurable;
use MiscUtils qw(split_list dump_die shell_cmd);
use UniqueCacheFile;

use LWP::Simple;
use HTTP::Request;
use LWP::UserAgent;
use URI::URL;

use Exporter;

use constant MAINTAINER_EMAIL => 'Michael.Edmonson@stjude.org';
use constant REQUEST_SLEEP_TIME => 1;
# In order not to overload the E-utility servers, NCBI recommends that
# users post no more than three URL requests per second and limit
# large jobs to either weekends or between 9:00 PM and 5:00 AM Eastern
# time during weekdays.

use constant URL_BASE => 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';

@EUtilities::ISA = qw(Configurable Exporter);
@EUtilities::EXPORT_OK = qw();

use MethodMaker qw(
        database
	retmode
	rettype
	retmax

	cache
        max_query_records
        ucf
	result_files
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->cache(1);
#  $self->max_query_records(100);
  $self->max_query_records(1000);
  $self->retmax(500);
#  $self->max_query_records(20);
  $self->ucf(new UniqueCacheFile("-unique" => 1,
				 "-md5_only" => 1,
				 "-prefix" => "eutilities",
	     ));
  $self->configure(%options);
  return $self;
}

sub fetch_genbank {
  # fetch a set of genbank records by accession
  my ($self, %options) = @_;
  $self->database("nucleotide");
#  $self->database("nuccore");
  $self->retmode("gb");
#  $self->rettype("text");
  $self->rettype("gb");
  return $self->fetch(%options);
}

sub fetch {
  # see http://www.ncbi.nlm.nih.gov/books/NBK25498/#chapter3.Application_3_Retrieving_large

  my ($self, %options) = @_;
  my $ids = $options{"-ids"} || die "-ids";
  # list of query items

  if (0) {
    print STDERR "****DEBUG: truncated list!\n";
    $ids = [ @{$ids}[0..10] ];
  }

  my $database = $self->database() || die "database";
  my $retmax = $self->retmax || die;
  my $retmode = $self->retmode || die "retmode";
  my $rettype = $self->rettype || die "rettype";
  my $retmax_query = $self->max_query_records();

  my %pairs_base;
  $pairs_base{db} = $database;
  $pairs_base{tool} = "EUtilities.pm";
  $pairs_base{email} = MAINTAINER_EMAIL;
  $pairs_base{usehistory} = "y";
  # use history server so we don't have to fetch and resubmit IDs

  my $query_url = URL_BASE . "esearch.fcgi";
  my $fetch_url = URL_BASE . "efetch.fcgi";

  my @result_files;
  $self->result_files(\@result_files);

  foreach my $set (split_list($ids, $retmax_query)) {
    my %pairs_q = %pairs_base;
    $pairs_q{term} = join ",", @{$set};
    $pairs_q{retmax} = $retmax_query;
    # 2/2019: now needed, otherwise default is now apparently 20
    my $output_ref = $self->polite_request("-url" => $query_url,
					   "-pairs" => \%pairs_q,
					   "-content" => 1,
					   "-post" => 1);

    my $web = $1 if ($$output_ref =~ /<WebEnv>(\S+)<\/WebEnv>/);
    my $key = $1 if ($$output_ref =~ /<QueryKey>(\d+)<\/QueryKey>/);
    my $count = $1 if ($$output_ref =~ /<Count>(\d+)<\/Count>/);

    for (my $retstart = 0; $retstart < $count; $retstart += $retmax) {
      my %pairs_r = %pairs_base;
      $pairs_r{WebEnv} = $web || die;
      $pairs_r{query_key} = $key || die;
      $pairs_r{retstart} = $retstart;
      $pairs_r{retmax} = $retmax;
#      $pairs_r{RetMax} = $retmax;
      $pairs_r{retmode} = $retmode;
      $pairs_r{rettype} = $rettype;

      my $of = $self->polite_request("-url" => $fetch_url,
				     "-pairs" => \%pairs_r
				    );
      push @result_files, $of;
    }
  }

  return \@result_files;
}

sub polite_request {
  # submit query/fetch requests to NCBI servers, sleeping between
  # submissions and and caching results
  my ($self, %options) = @_;
  my $url = $options{"-url"} || die;
  my $pairs = $options{"-pairs"};
  my $ucf = $self->ucf();

  $ucf->reset();
  $ucf->add($url);
  my @pairs;
  if ($pairs) {
    @pairs = map {$_, $pairs->{$_}} sort keys %{$pairs};
    $ucf->add(\@pairs)
  }

  my $cache_file = $ucf->get_file();
  unless (-s $cache_file) {
    if ($pairs) {
      # pairs of values given: use POST
      # (I think GET will eventually break for very long URLs)

      if ($options{"-post"}) {
	# might need for original query as might be too long for GET
	my $url_o = new URI::URL($url);

	my $request = new HTTP::Request(POST => $url);
	$request->content_type("application/x-www-form-urlencoded");
	$url_o->query_form(@pairs);
	$request->content($url_o->equery());

	my $ua = new LWP::UserAgent();
	my $response = $ua->request($request, $cache_file);
	die "error" if $response->is_error;
      } else {
	# larger requests seem to get truncated while fetching results, ???
	# use curl to copy URL instead
	my $url_o = new URI::URL($url);
	$url_o->query_form(@pairs);

	my $cmd = sprintf 'curl -o %s "%s"', $cache_file, $url_o->as_string();
	shell_cmd("-cmd" => $cmd);
      }
    } else {
      getstore($url, $cache_file);
    }
    die "ERROR fetching $url" unless -s $cache_file;
    sleep REQUEST_SLEEP_TIME;
  }

  my $result;
  if ($options{"-content"}) {
    local $/ = undef;
    open(EUTMP, $cache_file) || die;
    my $blob = <EUTMP>;
    $result = \$blob;
    close EUTMP;
  } else {
    $result = $cache_file;
  }
  return $result;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
