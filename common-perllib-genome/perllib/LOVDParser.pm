package LOVDParser;
# generic parser for HTML databases using LOVD
# - RB1, APC, ...?
# MNE 2/2014
#
# this is meant for the "all contents" view rather than the "unique variants"
# view.  The complete list contains duplicate entries, but also 
# also additional PubMed IDs which we are very interested in.
#

use strict;

use Configurable;
use Exporter;

use HTMLUtils qw(parse_html_tables);
use MiscUtils qw(dump_die);
use Reporter;

use HTML::TableExtract;
#              my $stripper = HTML::TableExtract::StripHTML->new;
#              $target = $stripper->strip($item);

@LOVDParser::ISA = qw(Configurable Exporter);
@LOVDParser::EXPORT_OK = qw();

use PostprocessedVariantReport qw(
F_GENE
F_REFSEQ
F_AACHANGE
);

my $FIELD_PUBMED = "PubMed";

my $FIELD_PATHOGENIC_CLEANED_SUBMITTED = "Pathogenicity_submitted_cleaned";
my $FIELD_PATHOGENIC_CLEANED_CONCLUDED = "Pathogenicity_concluded_cleaned";

my %PATH2SJ = (
	       "+" => "P",
	       "+?" => "LP",
	       "?" => "U",
	       "-?" => "LB",
	       "-" => "B",
	      );
# map LOVD codes to SJ 5-tier

use MethodMaker qw(
		    general_information
		    headers
		    rows

		    summary_field_genomic_refseq

		    field_reference
		    standardize_fields
		    field_aachange
		 );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->standardize_fields(1);
  $self->summary_field_genomic_refseq("Genomic refseq ID");
  $self->configure(%options);
  return $self;
}

sub parse {
  my ($self, %options) = @_;
  my $directory = $options{"-directory"} || die;
  my $limit = $options{"-limit"};

  my @files = glob(sprintf "%s/*.htm", $directory);

  my $f_homepage;
  my @f_data;
  foreach my $file (@files) {
    if ($file =~ /homepage/) {
      $f_homepage = $file;
    } else {
      push @f_data, $file;
    }
  }
  die "can't find homepage" unless $f_homepage;
  die "can't find data files" unless @f_data;

  #
  #  parse general info from homepage:
  #
  my $tables = parse_html_tables(
				 "-file" => $f_homepage,
				 "-headers" => 0,
				);

  my @wanted = grep {$_ and @{$_} and $_->[0]->[0] =~ /general information/i} @{$tables};
  die unless @wanted == 1;
  # general info table

  my %info;
  foreach my $row (@{$wanted[0]}) {
    $info{$row->[0]} = $row->[1];
  }
  $self->general_information(\%info);

  #
  #  parse data rows:
  #
  my @final_rows;
  my $total_entries;
  my $final_headers;

  my $reference_header = $self->field_reference || die "-field_reference";

  INFILE:
  foreach my $f (@f_data) {
    printf STDERR "parsing %s...\n", $f;
    my @po = (
	      "-file" => $f,
	      "-headers" => 0,
	      # won't work here
	      #	"-require-header" => "g-position"
	     );

    my $tables = parse_html_tables(@po);

    #
    #  extract headers
    #
    my @headers;
    my $enabled;

    foreach my $rows (@{$tables}) {
      next unless @{$rows};
      if (@{$rows} == 1) {
	my $thing = $rows->[0]->[0];
	$enabled = 1 if $thing =~ /Path\./;
	if ($thing =~ /^(\d+) public entries/) {
	  # "all contents" view: appears first
	  $total_entries = $1;
	  last;
	}

	if ($enabled) {
	  my @words = split /\s+/, $thing;
	  if (@words == 1 and $thing =~ /\w/) {
	    $thing =~ s/[^\w\-\#\.]+/_/g;
#	    printf STDERR "parsing header row: %d, %s\n", scalar(@words), $thing;
	    push @headers, $thing;
	  } else {
#	    printf STDERR "skipping row with %d words\n", scalar @words;
	  }
	}
      } else {
#	printf STDERR "skipping row of %d, first=%s second=%s\n", scalar @{$rows}, join(",", @{$rows->[0]}), join(",", @{$rows->[1]});
      }
    }

    unless ($total_entries) {
      # might not be in a table, e.g. LOVD MSH2
      open(SCANTMP, $f) || die;
      while (<SCANTMP>) {
	if (/(\d+) public entries/) {
	  $total_entries = $1;
	  last;
	}
      }
      close SCANTMP;
    }

    die "can't identify total entries" unless $total_entries;

    printf STDERR "headers: %s\n", join ",", @headers;

    unless ($final_headers) {
      $final_headers = [ @headers ];
    }

    my $tables_html = parse_html_tables(@po, "-keep-html" => 1);
    printf STDERR "parsing %s (HTML)...\n", $f;

    my $data_table;
    foreach my $rows (@{$tables_html}) {
      next unless @{$rows};
      if (@{$rows} > 3) {
	die "multiple data tables" if $data_table;
	$data_table = $rows;
      }
    }

    my $parsed = 0;
    foreach my $r (@{$data_table}) {
      if (@headers == @{$r}) {
	# column counts match
	my %r;
	@r{@headers} = @{$r};

#	printf STDERR "RAW: %s\n", join ",", @{$r};

	my $usable = 1;

#	printf STDERR "RAW:%s\n", join "\t", @{$r};

	die "no field $reference_header " . join "\n", @headers unless exists $r{$reference_header};
	$r{$FIELD_PUBMED} = $r{$reference_header} =~ /nih\.gov\/pubmed\/(\d+)/ ? $1 : "";
	my $blank = 1;
#	printf STDERR "PATH=%s\n", $r{"Path."};

	foreach my $h (@headers) {
	  if ($r{$h} =~ />/) {
	    my $raw = $r{$h};
#	    printf STDERR "possible HTML: before=%s, ", $raw;
	    if (1) {
	      my $stripper = HTML::TableExtract::StripHTML->new();
	      $r{$h} = $stripper->strip($raw);
	      # huzzah
	    } else {
	      $r{$h} =~ s/<[^>]+>//g;
	      # crude HTML strip
	      # FAIL:
	      # <a href="https://googledrive.com/host/0B8HVsr5izQxJUi1XTzEtWFlRc00/index.html?gene=MSH2&variant=c.-433T>G" target="_blank">Class 1</a>
	      #  (dies because > appears in a string)
	    }
#	    printf STDERR "after=%s\n", $r{$h};
	  }

	  $r{$h} =~ s/^\s+//;
	  $r{$h} =~ s/\s+$//;
	  $blank = 0 if $r{$h};
	}
	$usable = 0 if $blank;

	if ($usable) {
	  push @final_rows, \%r;
	  $parsed++;
	  last INFILE if $limit and $parsed >= $limit;
	}
      } else {
	die sprintf "header/data sync problem: expected=%d found=%d",
	scalar(@headers), scalar @{$r};
      }
    }
  }

  #  die unless $genbank_id;

  die sprintf "expected %d entries, parsed %d", $total_entries, scalar @final_rows unless $total_entries == @final_rows or $limit or $total_entries == 0;

  $self->headers($final_headers);
  $self->rows(\@final_rows);
}

sub export_to_file {
  #
  #  export database to flatfile.  Generally a passthrough,
  #  but with some additions/renaming of known SJ fieldnames.
  #
  my ($self, %options) = @_;
  my $info = $self->general_information();
  my $outfile = $options{"-file"};
  unless ($outfile) {
    my @things;
    foreach my $f (
      "Gene symbol",
      "Database location",
      "Version"
	) {
      my $v = $info->{$f} || die "no value for $f";
      $v =~ s/\s+/_/g;
      # version may contain gene name and a space, e.g. MSH2
      push @things, $v;
    }
    $outfile = sprintf 'export_%s.tab', join "_", @things;
  }

  my $field_genomic_refseq = $self->summary_field_genomic_refseq || die;

  my $gene_sym = $info->{"Gene symbol"} || die "no gene symbol";
  my $nm = $info->{"Transcript refseq ID"} || die "no transcript refseq id";
  my $genomic_refseq = $info->{$field_genomic_refseq} || die;
  # may need field abstraction

  my @labels;
  push @labels, F_GENE, F_REFSEQ, $field_genomic_refseq;
  my $local_aa = $self->field_aachange();
  if ($self->standardize_fields) {
    die "standardizing headers: neede -field_aachange" unless $local_aa;
  }
  my %label_map;

  my $standardize = $self->standardize_fields();

  foreach my $h (@{$self->{headers}}) {
    if ($standardize) {
      if ($h eq $local_aa) {
	$label_map{$h} = F_AACHANGE;
	push @labels, F_AACHANGE;
      } else {
	push @labels, $h;
      }
    } else {
      push @labels, $h;
    }
  }
  push @labels, $FIELD_PUBMED,
    $FIELD_PATHOGENIC_CLEANED_SUBMITTED,
      $FIELD_PATHOGENIC_CLEANED_CONCLUDED;

  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => \@labels,
			 "-auto_qc" => 1
			);

  foreach my $row (@{$self->rows}) {
    $row->{F_GENE()} = $gene_sym;
    $row->{F_REFSEQ()} = $nm;
    $row->{$field_genomic_refseq} = $genomic_refseq;
    $self->add_cooked_pathogenicity($row);

    if ($standardize) {
      die "ERROR: field_aachange $local_aa doesn't exist" unless exists $row->{$local_aa};
    }

    foreach my $from (keys %label_map) {
      my $to = $label_map{$from};
      $row->{$to} = $row->{$from};
    }

    foreach (values %{$row}) {
      s/\r?\n/ /g;
    }

    $rpt->end_row($row);
  }
  $rpt->finish();

  return $outfile;
}

sub add_cooked_pathogenicity {
  my ($self, $row) = @_;
  my $p_raw = $row->{"Path."} || die;
  my @p = split /\//, $p_raw;
  die unless @p == 2;
  $row->{$FIELD_PATHOGENIC_CLEANED_SUBMITTED} = cook_code($p[0]);
  $row->{$FIELD_PATHOGENIC_CLEANED_CONCLUDED} = cook_code($p[1]);
}

sub cook_code {
  my ($code_raw) = @_;
  die "no raw code" unless $code_raw;
  my $cooked = $PATH2SJ{$code_raw} || die "can't translate path code $code_raw";
  return $cooked;
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
