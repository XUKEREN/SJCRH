#!/bin/env perl
# helper for reclassifying data exported from csith
# MNE 12/2014

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;

use MiscUtils qw(dump_die);
use DelimitedFile;
use Reporter;
use AtomicOutfile;
use DelimitedMux;

my %FLAGS;

my $DIFF_MEDAL_FIELD_ONLY = 0;
# if set, deltas only produced if medal field is different
# (reason/evidence fields ignored)

my @options = (
	       "-genome=s",

	       "-prep",
	       # step 1:
	       # parse csith-exported files and extract into files
	       # suitable for reclassification runs

	       "-csith-snv=s",
	       "-csith-cnv=s",
	       "-csith-sv=s",

	       # step 2: reclassify prep'd files
	       "-reclassify-all",
	       "-reclassify-cnv",
	       "-reclassify-snv",
	       "-reclassify-sv",

	       "-compare-all",
	       # step 3: compare results, reporting differences
	       "-compare-snv",
	       "-compare-sv",
	       "-compare-cnv",

	       "-extend",
	       # add extended info to delta report:
	       # breaks format but adds info useful for debugging

	       "-snv-medal-clash",
	       # debug:
	       # if we bucket SNVs by genomic information only (i.e.
	       # excluding Class annotation), how many inconsistent
	       # medals do we find?
	       "-compare-hack",
	       # debug

	      );
GetOptions(\%FLAGS, @options);

my $MUX_SPLIT_LINES = 3000;
my $SNV_SOMATIC_RAM = 4000;
my $SNV_GERMLINE_RAM = 2500;

my @MEDAL_COLUMN_REMAP = (
			  "medal" => "csith_medal",
			  "medal_reason" => "csith_medal_reason",
			  "medal_evidence" => "csith_medal_evidence",
			 );

my %SITH2MEDAL_SNV = (
		      "gene" => "GeneName",
		      "chr" => "Chr",
		      "position" => "WU_HG19_Pos",
		      "reference_allele" => "ReferenceAllele",
		      "mutant_allele" => "MutantAllele",
		      "class" => "Class",
		      "variant" => "AAChange",
		      "accession" => "mRNA_acc",
		      @MEDAL_COLUMN_REMAP,
		     );
my @UNIQUE_FIELDS_SNV = (
			"Chr",
			"WU_HG19_Pos",
			"ReferenceAllele",
			"MutantAllele",

			 "Class",
			 # there are some uniqueness problems!
			 # e.g. chr1.120612005.GG.-
			 # is sometimes annotated as a frameshift (gold)
			 # and sometimes an exon (bronze)!
		       );
# this "should" be enough, but who knows, maybe there are additional
# variants in Class annotation, etc.  auto-QC should detect medal
# discrepancies anyway.


my @UNIQUE_FIELDS_SV = (
			"Fusion Gene",
			"Usage",
		       );
# classifier fields required to uniquely identify a SV


my %SITH2MEDAL_SV = (
		     "fusion_gene" => "Fusion Gene",
		     "usage" => "Usage",
		      @MEDAL_COLUMN_REMAP,
		    );
# translate csith column names to medal input column names
# (or rename so not to conflict)

my @UNIQUE_FIELDS_CNV = qw(
			    chrom
			    loc.start
			    loc.end
			    seg.mean
			    LogRatio
		       );
my %SITH2MEDAL_CNV = (
		      "logratio" => "LogRatio",
		      @MEDAL_COLUMN_REMAP,
		    );

my $RC_SNV = "reclassify_snv.tab";
my $RC_CNV_SOMATIC = "reclassify_cnv_somatic.tab";
my $RC_CNV_GERMLINE = "reclassify_cnv_germline.tab";
my $RC_SV = "reclassify_sv.tab";

if ($FLAGS{prep}) {
  prep_files();
} elsif ($FLAGS{"reclassify-cnv"}) {
  reclassify_cnv();
} elsif ($FLAGS{"reclassify-snv"}) {
  reclassify_snv();
} elsif ($FLAGS{"reclassify-sv"}) {
  reclassify_sv();
} elsif ($FLAGS{"reclassify-all"}) {
  reclassify_all();
} elsif ($FLAGS{"compare-snv"}) {
  compare_snv();
} elsif ($FLAGS{"compare-sv"}) {
  compare_sv();
} elsif ($FLAGS{"compare-cnv"}) {
  compare_cnv();
} elsif ($FLAGS{"compare-hack"}) {
  compare_hack();
} elsif ($FLAGS{"compare-all"}) {
  compare_all();
} elsif ($FLAGS{"snv-medal-clash"}) {
  snv_medal_clash();
} else {
  printf STDERR "ERROR: option?:\n";
  foreach (@options) {
    printf STDERR "  %s\n", $_;
  }
  exit(1);
}

sub prep_files {
  # step 1:
  # parse csith-exported files and convert into the format
  # required by classifier

  prep_file(
	    "-param" => "csith-snv",
	    "-map" => \%SITH2MEDAL_SNV,
	    "-out" => $RC_SNV,
	   );

  prep_file(
	    "-param" => "csith-sv",
	    "-map" => \%SITH2MEDAL_SV,
	    "-out" => $RC_SV,
	   );

  prep_file(
	    "-param" => "csith-cnv",
	    "-map" => \%SITH2MEDAL_CNV,
	    "-out" => $RC_CNV_SOMATIC,
	    "-somatic" => 1,
	   );

  prep_file(
	    "-param" => "csith-cnv",
	    "-map" => \%SITH2MEDAL_CNV,
	    "-out" => $RC_CNV_GERMLINE,
	    "-germline" => 1,
	   );

}

sub prep_file {
  my (%options) = @_;
  my $param = $options{"-param"} || die;
  my $map = $options{"-map"} || die;
  my $infile = get_file_param($param);
  my $outfile = $options{"-out"} || die;
  my $somatic = $options{"-somatic"};
  my $germline = $options{"-germline"};

  printf STDERR "processing %s...\n", $infile;

  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			     );

  my @h2 = map {$map->{$_} || $_} @{$df->headers_raw()};
  # rename medal column names to not conflict with columns added
  # by classifier

  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-labels" => \@h2,
			 "-delimiter" => "\t",
			 "-auto_qc" => 1,
			);

  while (my $row = $df->get_hash()) {
    next if somatic_germline_filter(%options, "-row" => $row);

    foreach my $k (sort keys %{$row}) {
      if (my $k2 = $map->{$k}) {
	$row->{$k2} = $row->{$k};
      }
    }
    $rpt->end_row($row);
  }
  $rpt->finish();
}

sub compare_sv {
  compare_file(
	       "-param" => "csith-sv",
	       # param for original csith file
	       "-reclassified" => "reclassify_sv.tab.medals.tab",
	       # reclassified results
	       "-unique-fields" => \@UNIQUE_FIELDS_SV,
	       # fields required to bucket variants uniquely
	       "-out" => "reclassify_sv_delta.tab",
	       # output file (delta)
	       "-sith2medal" => \%SITH2MEDAL_SV,
	       # convert csith report column names to medal column names
	      );
}

sub compare_cnv {
  compare_file(
	       "-param" => "csith-cnv",
	       "-reclassified" => "reclassify_cnv_somatic.tab.medals.tab",
	       "-unique-fields" => \@UNIQUE_FIELDS_CNV,
	       "-out" => "reclassify_cnv_somatic_delta.tab",
	       "-sith2medal" => \%SITH2MEDAL_CNV,
	       "-somatic" => 1,
	      );

  compare_file(
	       "-param" => "csith-cnv",
	       "-reclassified" => "reclassify_cnv_germline.tab.medals.tab",
	       "-unique-fields" => \@UNIQUE_FIELDS_CNV,
	       "-out" => "reclassify_cnv_germline_delta.tab",
	       "-sith2medal" => \%SITH2MEDAL_CNV,
	       "-germline" => 1,
	      );
}

sub compare_hack {
  compare_file(
	       "-param" => "csith-cnv",
#	       "-reclassified" => "reclassify_cnv_germline.tab.medals.tab",
	       "-reclassified" => "example_germline_cnv_was_silver_now_unknown.tab",
	       "-unique-fields" => \@UNIQUE_FIELDS_CNV,
	       "-out" => "reclassify_cnv_germline_delta.tab",
	       "-sith2medal" => \%SITH2MEDAL_CNV,
	       "-germline" => 1,
	      );
}

sub compare_snv {
  compare_file(
	       "-param" => "csith-snv",
	       "-reclassified" => "reclassify_snv.tab.medals.tab",
	       # somatic medal calls
	       "-unique-fields" => \@UNIQUE_FIELDS_SNV,
	       "-out" => "reclassify_snv_somatic_delta.tab",
	       "-sith2medal" => \%SITH2MEDAL_SNV,
	       "-somatic" => 1,
	      );

  compare_file(
	       "-param" => "csith-snv",
	       "-reclassified" => "reclassify_snv.tab.medals.tab.medals.tab",
	       # germline medal calls (GSBClass now contains germline)
	       "-unique-fields" => \@UNIQUE_FIELDS_SNV,
	       "-out" => "reclassify_snv_germline_delta.tab",
	       "-sith2medal" => \%SITH2MEDAL_SNV,
	       "-germline" => 1,
	      );


}

sub compare_all {
  compare_sv();
  compare_cnv();
  compare_snv();
}

sub compare_file {
  my (%options) = @_;
  my $param = $options{"-param"} || die;
  my $fn_csith = get_file_param($param);
  my $fn_reclassified = $options{"-reclassified"} || die;
  my $unique_fields = $options{"-unique-fields"} || die;
  my $sith2medal = $options{"-sith2medal"} || die;
  my $outfile = $options{"-out"} || die;
  printf STDERR "generating %s...\n", $outfile;

  my $extend_mode = $FLAGS{extend};

  $outfile .= ".extended" if $extend_mode;

  die "where is $fn_reclassified?" unless -s $fn_reclassified;

  #
  # bucket new calls:
  #
  my $new_calls = bucket_calls(
			       "-file" => $fn_reclassified,
			       "-fields" => $unique_fields
			      );

  #
  # parse through csith calls and compare:
  #
  my $df = new DelimitedFile(
			     "-file" => $fn_csith,
			     "-headers" => 1,
			     );

  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-auto_qc" => 1,
			      "-extra" => $extend_mode ?
			      [qw(csith_medal csith_medal_reason csith_medal_evidence)] : []
			     );

  my $medal2sith = invert_hash($sith2medal);

  while (my $row = $df->get_hash()) {
    next if somatic_germline_filter(%options, "-row" => $row);

    my $key = build_key(
			"-row" => $row,
			"-fields" => $unique_fields,
			"-map" => $medal2sith
		       );
#    dump_die($row, $key);
    my $new_set = $new_calls->{$key} || die "no new calls for $key";
    my $new_row = $new_set->[0];
    # all rows guaranteed by bucketing QC to have same medal call

    die "WTF" unless exists $row->{medal} and exists $new_row->{GSBClass};
    my $medal_csith = $row->{medal};
    my $medal_reclassify = $new_row->{GSBClass};

    $medal_reclassify = "" if $medal_reclassify eq "Unknown";
    $medal_reclassify = uc($medal_reclassify);
    # update medal call to look like it does in csith

    my $changed;
    $changed = 1 if $medal_csith ne $medal_reclassify;
    # reclassification produces different medal

    unless ($DIFF_MEDAL_FIELD_ONLY) {
      # also compare the Reason and Evidence fields.

#      dump_die($new_row);

      my %re;
      $re{medal_reason} = "Reason";
      $re{medal_evidence} = "Evidence";

      foreach my $f_csith (keys %re) {
	my $f_new = $re{$f_csith};
	die unless exists $row->{$f_csith} and exists $new_row->{$f_new};

	my $v_csith = lc($row->{$f_csith});
	my $v_new = lc($new_row->{$f_new});

	if ($v_csith ne $v_new) {
	  # change
	  dump_die($new_row, "non-medal change in $f_new: csith:$v_csith new:$v_new", 1) unless $changed;
	  $changed = 1;
	}

      }


    }

    if ($changed) {
      if ($extend_mode) {
	foreach my $k (keys %{$row}) {
	  if ($k =~ /medal/) {
	    my $k2 = "csith_" . $k;
	    $row->{$k2} = $row->{$k};
	  }
	}
      }

      $row->{medal} = $medal_reclassify;
      $row->{medal_reason} = $new_row->{Reason};
      # DO NOT UPPERCASE:
      # reason fields may contain case-sensitive gene symbols
      $row->{medal_evidence} = $new_row->{Evidence};
      # not sure whether case ever matters here
      $rpt->end_row($row);
    }
  }
  $rpt->finish();

}

sub get_file_param {
  my ($param) = @_;
  my $file = $FLAGS{$param} || die "-" . $param;
  die "where is $file" unless -s $file;
  return $file;
}

sub bucket_calls {
  my (%options) = @_;
  my $file = $options{"-file"} || die;
  my $clash_mode = $options{"-clash"};

  my %bucket;

  my $df = new DelimitedFile("-file" => $file,
			     "-headers" => 1,
			     );
  while (my $row = $df->get_hash()) {
    next if somatic_germline_filter(%options, "-row" => $row);
    my $key = build_key(
			%options,
			"-row" => $row,
		       );
    push @{$bucket{$key}}, $row;
  }

  # QC check: ensure medals are consistent for the set of fields
  # that are supposedly to uniquely define a variant:
  my %clash;
  foreach my $key (keys %bucket) {
    my $set = $bucket{$key};
    if (@{$set} > 1) {
      my %medals = map {($_->{GSBClass} || die), 1} @{$set};
      unless (scalar keys %medals == 1) {
	if ($clash_mode) {
	  $clash{$key} = 1;
	} else {
	  die "medal ambiguity for $key!";
	}
      }
    }
  }

  return $clash_mode ? \%clash : \%bucket;
}

sub build_key {
  my (%options) = @_;
  my $row = $options{"-row"} || die;
  my $fields = $options{"-fields"} || die;
  my $map = $options{"-map"};

  my @tokens;
  foreach my $field (@{$fields}) {
    my @try = $field;
    push @try, $map->{$field} if $map and $map->{$field};

    my $value;
    foreach my $f (@try) {
      $value = $row->{$f};
      last if defined $value;
    }
    dump_die($row, sprintf "no value for %s", join ",", @try) unless defined $value;
    push @tokens, $value;
  }
  my $key = join "_", @tokens;
  return $key;
}

sub invert_hash {
  my ($hash) = @_;
  my %inv;
  foreach my $key (keys %{$hash}) {
    my $value = $hash->{$key};
    die unless defined $value;
    die "duplicate" if exists $inv{$value};
    $inv{$value} = $key;
  }
  return \%inv;
}

sub somatic_germline_filter {
  my (%options) = @_;
  my $row = $options{"-row"} || die;
  my $somatic = $options{"-somatic"};
  my $germline = $options{"-germline"};

  my $filter;
  if ($somatic or $germline) {
    my $origin = $row->{origin} || die;
    if ($origin eq "GERMLINE") {
      $filter = 1 unless $germline;
    } elsif ($origin eq "SOMATIC") {
      $filter = 1 unless $somatic;
    } else {
      die;
    }
  }
  return $filter;
}

sub reclassify_all {
  reclassify_sv();
  reclassify_cnv();
  reclassify_snv();
}

sub reclassify_sv {
  my $genome = $FLAGS{genome} || die "-genome";
  my $cmd = sprintf 'medal_ceremony_using_configs.sh %s -single-sv %s', $genome, $RC_SV;
  run_cmd($cmd);
}

sub reclassify_cnv {
  my $genome = $FLAGS{genome} || die "-genome";

  my @cmds;
  push @cmds, sprintf 'medal_ceremony_using_configs.sh %s -single-cnv %s',
    $genome, $RC_CNV_SOMATIC;

  push @cmds, sprintf 'medal_ceremony_using_configs.sh %s -single-cnv %s -cnv-germline',
    $genome, $RC_CNV_GERMLINE;

  foreach my $cmd (@cmds) {
    run_cmd($cmd);
  }
}

sub reclassify_snv {
  my $genome = $FLAGS{genome} || die "-genome";

  my $base = sprintf 'medal_ceremony_using_configs.sh %s -no-pvr', $genome;

  my $somatic_template = sprintf '%s -single-snv-indel %%s', $base;
  my $gl_template = sprintf '%s -single-gl %%s', $base;

  printf STDERR "running somatic SNV/indel classification...\n";
  my $somatic_out = mux_medals(
			       "-in" => $RC_SNV,
			       "-template" => $somatic_template,
			       "-ram" => $SNV_SOMATIC_RAM,
			      );

  printf STDERR "running germline SNV/indel classification...\n";
  my $gl_out = mux_medals(
			  "-in" => $somatic_out,
			  "-template" => $gl_template,
			  "-ram" => $SNV_GERMLINE_RAM,
			 );


}

sub mux_medals {
  my (%options) = @_;
  my $infile = $options{"-in"} || die;
  my $template = $options{"-template"} || die;
  my $ram = $options{"-ram"} || die;

  my $dm = new DelimitedMux();
  $dm->subdir_mode(1);
  # classifier wrapper can 'splode if multiple run attempts
  # in the same directory

  my $split_files = $dm->split_file(
				    "-file" => $infile,
				    "-lines" => $MUX_SPLIT_LINES,
				   );
  # demux

  $dm->run_jobs(
		"-template" => $template,
		"-out-suffix" => ".medals.tab",
		"-ram" => $ram,
		"-wait" => 1,
		# wait for all jobs to complete
	       );

  my $joined = $infile . ".medals.tab";

  $dm->join_files(
		  "-out" => $joined
		 );
  # mux

  return $joined;
}


sub run_cmd {
  my ($cmd) = @_;
  printf STDERR "running %s...\n", $cmd;
  system $cmd;
  die "ERROR" if $?;
}

sub snv_medal_clash {

  my @genomic_fields = (
			"Chr",
			"WU_HG19_Pos",
			"ReferenceAllele",
			"MutantAllele",
		       );
  # just genomic information, exclude Class annotation

  my $somatic = bucket_calls(
			     "-file" => "reclassify_snv.tab.medals.tab",
			     "-somatic" => 1,
			     "-fields" => \@genomic_fields,
			     "-clash" => 1
			    );
  printf STDERR "somatic: %d\n", scalar keys %{$somatic};

  my $gl = bucket_calls(
			"-file" => "reclassify_snv.tab.medals.tab.medals.tab",
			"-germline" => 1,
			"-fields" => \@genomic_fields,
			"-clash" => 1
		       );
  printf STDERR "germline: %d\n", scalar keys %{$gl};
  foreach (sort keys %{$gl}) {
    printf STDERR "  %s\n", $_;
  }


}
