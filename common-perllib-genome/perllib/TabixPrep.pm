package TabixPrep;
# helper to prepare lar

use strict;
use Configurable;
use Exporter;
use DelimitedFile;
use MiscUtils qw(log_message dump_die);
use WorkingFile;
use FileUtils qw(find_binary);

@TabixPrep::ISA = qw(Configurable Exporter);
@TabixPrep::EXPORT_OK = qw(tabix_concatenate);

use MethodMaker qw(
	outfile
headers
header_chr
header_start
header_end

cook_chrom_names
tfw
spill
compress_intermediate_files
use_scratch
delete_tempfiles
input_sorted

vcf_headers
		  );
use GenomeUtils qw(cook_chromosome_name);
use FileUtils qw(find_binary);
use TemporaryFileWrangler;
use TabixFile;
use Reporter;
use FileUtils qw(universal_open);

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->compress_intermediate_files(1);
  # on the fence
  $self->use_scratch(1);
  $self->delete_tempfiles(1);
  $self->configure(%options);
  find_binary("bgzip", "-die" => 1);
  # hacky but for whatever reason seems difficult to detect
  # pipe errors if this command is not on path (??)
  find_binary("tabix", "-die" => 1);
  my $tfw = new TemporaryFileWrangler();
  unless ($self->delete_tempfiles) {
    printf STDERR "DEBUG: verbose / no auto-unlink tempfiles\n";
    $tfw->verbose(1);
    $tfw->auto_unlink(0);
  }
  $self->tfw($tfw);
  $self->spill({});

  my $outfile = $self->outfile || die "outfile";
  die "ERROR: outfile $outfile must end in .gz" unless $outfile =~ /gz$/;

  return $self;
}

sub add_row {
  my ($self, %options) = @_;
  my $spill = $self->spill();

  my $chr_cooked;
  # for tracker
  my $line = $options{"-line"};
  my $row = $options{"-row"};

  if ($row) {
    my $f_chr = $self->header_chr() || die "header_chr";
    if ($self->cook_chrom_names()) {
      $row->{$f_chr} = cook_chromosome_name($row->{$f_chr});
    }
    $chr_cooked = cook_chromosome_name(($row->{$f_chr} || dump_die($row, "no $f_chr")), "-return-unknown" => 1);
  } elsif ($line) {
    $chr_cooked = $options{"-chrom"} || die "-chrom";
  } else {
    die "need -row or -line";
  }

  unless ($spill->{$chr_cooked}) {
    my $suffix = sprintf ".tabix_prep.%s", $chr_cooked;
    my $cif = $self->compress_intermediate_files();
    $suffix .= ".gz" if $cif;
    my $outfile;
    if ($self->use_scratch) {
      $outfile = $self->tfw->get_tempfile("-append" => $suffix);
    } else {
      $outfile = sprintf 'tabix_prep.%s%s', $$, $suffix;
      $self->tfw->add_tempfile($outfile) if $self->delete_tempfiles;
    }
    my $rpt;
    if ($row) {
      my @ro = (
	"-file" => $outfile,
	"-delimiter" => "\t",
	"-labels" => ($self->headers || die "headers"),
	  );
      push @ro, "-compress" => "gz" if $cif;
      $rpt = new Reporter(@ro);
      $rpt->headers_done(1);
      # skip headers so we can use external sort.
      # might be more efficient for very large datasets.
    } elsif ($line) {
      my @wo;
      push @wo, "-compress" => "gz" if $cif;
      $rpt = new WorkingFile($outfile, @wo);
    } else {
      die;
    }

    $spill->{$chr_cooked}{file} = $outfile;
    $spill->{$chr_cooked}{rpt} = $rpt;
    printf STDERR "outfile for %s: %s\n", $chr_cooked, $outfile;
  }
  my $rpt = $spill->{$chr_cooked}{rpt} || die;
  if ($row) {
    $rpt->end_row($row);
  } elsif ($line) {
    my $fh = $rpt->output_filehandle;
    print $fh $$line;
  } else {
    die;
  }
}

sub get_wf_outfile {
  my ($self) = @_;
  my $outfile = $self->outfile || die;
  my $wf = new WorkingFile(
			   $outfile,
			   "-compress" => "bgzip"
			  );
  my $fh = $wf->output_filehandle();

  if (my $vh = $self->vcf_headers()) {
    # output in VCF format: write main header first
    foreach (@{$vh}) {
      printf $fh $_;
    }
  }

  my $headers = $self->headers || die;
  my $need_comment = $headers->[0] =~ /^#/ ? 0 : 1;

  printf $fh "%s%s\n", ($need_comment ? "#" : ""),
  join "\t", @{$headers};
  # comment header for tabix compatibility if req'd
  return $wf;
}

sub finish {
  my ($self) = @_;
  my $spill = $self->spill();
  my $wf = $self->get_wf_outfile();
  my $fh = $wf->output_filehandle();
  my $f_start = $self->header_start || die;
  my $f_end = $self->header_end;

  foreach my $chr (sort keys %{$spill}) {
    $spill->{$chr}{rpt}->finish();

    my $sort_cmd = sprintf "%s %s|", $self->compress_intermediate_files ? "gzip -dc" : "cat", $spill->{$chr}{file};
    unless ($self->input_sorted) {
      my $fno_start = $self->get_field_number($f_start) || die;
      if (0 and $f_end) {
	# don't think -k1,2 will work unless fields are adjacent and
	# in same order, and this is not required of input
	die "multi key sort, implement me, start key = $fno_start";
      } else {
#	$sort_cmd .= sprintf 'sort -nk%d|', $fno_start;
	# sort each chrom one at a time rather than simultaneously
	$sort_cmd .= sprintf "sort -nk%d -t \"\t\"|", $fno_start;
	# explicitly use tab as the field separator.
	# if we don't some files (e.g. ExAC vcf2tab in mode 2)
	# will fail due to field separator confusion (presence of 
	# some characters breaks sort's field separator assumptions??)
#	die $sort_cmd;
      }
    }

    printf STDERR "%s: processing %s (%s)...\n", scalar(localtime), $chr, $sort_cmd;

    open(TPTMP, $sort_cmd) || die "can't open $sort_cmd: $!";
    while (<TPTMP>) {
#      printf STDERR "sorted: %s\n", $_;
      print $fh $_;
    }
    close TPTMP;
  }
  $wf->finish();

  # index:
  $self->build_tabix_index();

}

sub build_tabix_index {
  my ($self) = @_;
  my $outfile = $self->outfile || die "outfile";
  unlink $outfile . ".tbi";
  # refresh
  my $f_chr = $self->header_chr || die;
  my $f_start = $self->header_start || die;
  my $f_end = $self->header_end();

  if ($self->vcf_headers) {
    my $cmd = sprintf 'tabix -p vcf %s', $outfile;
    system $cmd;
    die "error running $cmd" if $?;
    die unless -s $outfile . ".tbi";
  } else {
    my $tf = new TabixFile(
			   "-file" => $outfile,
			   "-index" => 1,
			   "-f_chr" => $f_chr,
			   "-f_start" => $f_start,
			   "-f_end" => $f_end,
			  );
  }
}



sub get_field_number {
  my ($self, $name) = @_;
  my $fno;
  my $headers = $self->headers;
  for (my $i = 0; $i < @{$headers}; $i++) {
    if ($headers->[$i] eq $name) {
      $fno = $i + 1;
      last;
    }
  }
  die "can't find header number for $name" unless $fno;
  return $fno;
}

sub init_headers_from_file {
  my ($self, $file) = @_;
  my $df = new DelimitedFile(
			     "-file" => $file,
			     "-headers" => 1,
			     );
  $self->headers($df->headers_raw);
}

sub tabix_concatenate {
  my (%options) = @_;
  my $infiles = $options{"-in"} || die "-in";
  my $outfile = $options{"-out"} || die "-out";
  my $sorted = $options{"-sorted"};
  die "-sorted" unless defined $sorted;
  my $max = $options{"-max"};
  my $ping = $options{"-ping"};

  my $h_chr = $options{"-header_chr"} || die "-header_chr";
  my $h_pos = $options{"-header_start"} || die "-header_start";

  if ($sorted) {
    # files are sorted, possibly already by TabixPrep (e.g. split processing).
    # very slow to parse rows and copy to intermediate files again.
    my $wf = new WorkingFile($outfile, "-compress" => "bgzip");
    my $fh_out = $wf->output_filehandle();
    my $first_file = 1;
    my $total_count = 0;
    foreach my $infile (@{$infiles}) {
      my $fh = universal_open($infile);
      if ($first_file) {
	# keep header line in first file...
	$first_file = 0;
      } else {
	# ...skip in second and later files
	my $hdr = <$fh>;
#	printf STDERR "skipping header %s\n", $hdr;
	
      }

      my $count = 0;
      while (<$fh>) {
	$count++;
	$total_count++;
	if ($ping and $total_count % $ping == 0) {
	  log_message(sprintf "processed %d, at %s", $count, $_);
	}
	last if $max and $count >= $max;
	print $fh_out $_;
      }
    }
    $wf->finish();

    my $tp = new TabixPrep(
      "-outfile" => $outfile,
      "-header_chr" => $h_chr,
      "-header_start" => $h_pos
	);
    $tp->init_headers_from_file($infiles->[0]);
    $tp->build_tabix_index();

  } else {
    my $tp = new TabixPrep(
      "-outfile" => $outfile,
      "-header_chr" => $h_chr,
      "-header_start" => $h_pos
	);
    $tp->init_headers_from_file($infiles->[0]);
    $tp->input_sorted($sorted);

    my $total_count = 0;
    foreach my $f (@{$infiles}) {
      my $count = 0;
      printf STDERR "%s...\n", $f;
      
      my $df = new DelimitedFile(
	"-file" => $f,
	"-headers" => 1,
	  );
#    $df->headers_raw->[0] =~ s/^#//;
      # undo tabix comment
      while (my $row = $df->get_hash()) {
	$tp->add_row("-row" => $row);
	$total_count++;
	$count++;
	if ($ping and $total_count % $ping == 0) {
	  log_message(sprintf "processed %d, at %s.%s", $count, @{$row}{$h_chr, $h_pos});
	}
	last if $max and $count >= $max;
      }
    }
    $tp->finish();
  }

  
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
