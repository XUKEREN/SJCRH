#!/bin/env perl
# split HGVS annotation into component chr/pos/ref/alt fields

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
use NucleotideSubstitutionParser;
use FAI;
use TdtConfig;
use RefFlatFile;
use TdtConfig;
use ConfigUtils qw(config_or_manual);
use GenomeUtils qw(reverse_complement);

my $F_OUT_CHR = "sj_chr";
my $F_OUT_POS = "sj_pos";
my $F_OUT_RA = "sj_ref_allele";
my $F_OUT_VA = "sj_alt_allele";
my $F_OUT_TYPE = "sj_var_type";
my $F_OUT_NOTE = "sj_conversion_note";
my $F_OUT_TRANSCRIPT = "sj_mrna_acc";
my $F_OUT_AA = "sj_aachange";

my $VERIFY_GENOMIC = 1;

my %FLAGS;
my @clopts = (
	      "-file=s",
	      "-split",
	      "-split-c",

	      "-field-hgvs=s",
	      "-field-hgvs-c=s",
	      "-field-hgvs-p=s",

	      "-field-chr=s",
	      "-field-chr-pos=s",
	      "-field-ra=s",
	      "-field-va=s",

	      "-chr=s",

	      "-extract-gpos",
	      "-strip-trailing-comment",
	      "-genome=s",
	      "-fasta=s",
	      "-patch-deletion-sequence",
	      "-strip-whitespace",
	      "-strip-non-ascii",
	      "-reformat1",
	      "-reformat-slash",
	      "-refflat=s",

	      "-rc",
	      "-extract-hgvs-chr",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

if ($FLAGS{split}) {
  hgvs_split_g();
} elsif ($FLAGS{"split-c"}) {
  hgvs_split_c();
}

sub hgvs_split_c {
  # HGVS c. (transcript-oriented alleles)
  my $infile = get_hash_option(\%FLAGS, "file");
  my $f_hgvs = get_hash_option(\%FLAGS, "field-hgvs-c");
  my $f_hgvs_p = $FLAGS{"field-hgvs-p"};
  my $f_chr_pos = $FLAGS{"field-chr-pos"};
  my $strip_whitespace = $FLAGS{"strip-whitespace"};
  my $strip_non_ascii = $FLAGS{"strip-non-ascii"};
  my $reformat1 = $FLAGS{reformat1};
  my $reformat_slash = $FLAGS{"reformat-slash"};
  my $f_ra = $FLAGS{"field-ra"};
  my $f_va = $FLAGS{"field-va"};

  my $genome = $FLAGS{genome} || die "-genome";
  my $f_refflat = config_or_manual(
				   "-config-type" => "genome",
				   "-config-name" => $genome,
				   "-parameter" => "REFSEQ_REFFLAT",
				   "-manual" => $FLAGS{refflat}
				  );

  my $f_fa = config_or_manual(
			      "-config-type" => "genome",
			      "-config-name" => $genome,
			      "-parameter" => "FASTA",
			      "-manual" => $FLAGS{fasta}
			     );
  my $fai = new FAI("-fasta" => $f_fa);

  my $rf = new RefFlatFile("-strip_sharp_annotations" => 1);
  $rf->parse_file(
		  "-refflat" => $f_refflat,
		  "-type" => "refgene",
		 );

  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  my $outfile = basename($infile) . ".hgvs_split.tab";
  my @extra;
  push @extra, $F_OUT_CHR, $F_OUT_POS if $f_chr_pos;
  # optional
  push @extra, $F_OUT_TRANSCRIPT;
  # from HGVSc
  push @extra, $F_OUT_AA if $f_hgvs_p;
  # if HGVSp available
  push @extra, $F_OUT_RA, $F_OUT_VA;
  # from HGVSc
  push @extra, $F_OUT_NOTE;

  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => \@extra,
  			      "-auto_qc" => 1,
			     );

  my $count_broken = 0;
  my $count_unknown = 0;
  my $count_ok = 0;

  my $nsp = new NucleotideSubstitutionParser();
  while (my $row = $df->get_hash()) {
    my $hgvs = $row->{$f_hgvs};
    $hgvs =~ s/\s+//g if $strip_whitespace;

    my $cpos;
    if ($hgvs =~ /(c\..*)/) {
      $cpos = $1;
    }

    my ($nm) = split /:/, $hgvs;
    $nm = "" unless defined $nm;

    my $parse_ok = 0;
    my $ra = "";
    my $va = "";
    my $chr = "";
    my $pos = "";
    my $aa = "";
    if ($f_hgvs_p) {
      my $hp = $row->{$f_hgvs_p};
      if ($hp) {
	my @f = split /:/, $hp;
	if (@f == 2) {
	  $aa = $f[1];
	} else {
	  $aa = $hp;
	}
	printf STDERR "WARNING: extracted %s from %s\n", $aa, $hp unless $aa =~ /^p\./;
      }
    }

    my @notes;
    if ($f_chr_pos) {
      if (my $v = $row->{$f_chr_pos}) {
	($chr, $pos) = split /:/, $v;
      }
      $row->{$F_OUT_CHR} = $chr;
      $row->{$F_OUT_POS} = $pos;
    }

    if ($cpos) {
      $cpos =~ tr/\x80-\xFF//d if $strip_non_ascii;
      # 128 and above

      if ($reformat1 and $cpos !~ />/ and
	  $cpos =~ /c\.([ACGT]+)(\d+)([ACGT]+)$/) {
	# NM_020919.3:c.G2002T => NM_020919.3:c.2002G>T
	my $new = sprintf 'c.%d%s>%s', $2, $1, $3;
	my $orig = $hgvs;
	$cpos =~ s/c\..*$/$new/ || die;
      }

      if ($reformat_slash and $cpos !~ />/ and $cpos =~ /\//) {
	$cpos =~ s/\//>/ || die;
      }

#      printf STDERR "parse: %s\n", $cpos;
      if ($nsp->parse($cpos)) {
	$count_ok++;
	$parse_ok = 1;

	$ra = $nsp->reference_sequence();
	$va = $nsp->variant_sequence();

	#
	# check for strand flip:
	#
	die "where is NM" unless $nm;
	my $strand = $rf->get_strand_for_accession($nm) || die "can't find strand for $nm";
	if ($strand eq "+") {
	  # OK as is
	} elsif ($strand eq "-") {
	  # bases are in transcript orientation, so if transcript is on -,
	  # reverse-complement to get genomic bases
	  push @notes, "transcript_on_antisense";
	  foreach ($ra, $va) {
	    $_ = reverse_complement($_);
	  }
	} else {
	  die;
	}

	if ($chr and $pos and $ra and $ra ne "-") {
	  #
	  # verify genomic sequence for non-insertions:
	  #
	  my $chunk = $fai->get_chunk(
				      "-id" => $chr,
				      "-start" => $pos,
				      "-length" => length($ra),
				     );
	  $ra = $chunk if $ra =~ /^N+$/;
	  # details of sequence not specified

	  if ($ra eq $chunk) {
#	    print STDERR "reference sequence verification OK\n";
	  } else {
	    push @notes, "reference verification failed ($ra vs $chunk)";
	  }
	}

      } else {
	$count_broken++;
	print STDERR "ERROR: can't parse \"$cpos\"\n";
      }
    } else {

      if ($f_ra) {
	my $v = $row->{$f_ra};
	if ($v eq "" or $v =~ /^[ACGT]+$/) {
	  $ra = $v;
	  push @notes, "rescue ref allele from separate column";
	} else {
	  push @notes, "rescue ref allele failed";
	}
      }

      if ($f_va) {
	my $v = $row->{$f_va};
	if ($v eq "" or $v =~ /^[ACGT]+$/) {
	  $va = $v;
	  push @notes, "rescue alt allele from separate column";
	} else {
	  push @notes, "rescue alt allele failed";
	}
      }
      $count_unknown++ unless $f_ra or $f_va;
    }

    if ($ra !~ /\w/ and $va !~ /\w/) {
      $ra = $va = "";
      push @notes, "broken allele pairing set to blank";
    }

    $row->{$F_OUT_TRANSCRIPT} = $nm;
    $row->{$F_OUT_AA} = $aa;
    $row->{$F_OUT_RA} = $ra;
    $row->{$F_OUT_VA} = $va;
    $row->{$F_OUT_NOTE} = join ",", @notes;
    $rpt->end_row($row);
  }

  $rpt->finish();

}

sub hgvs_split_g {
  # HGVS g. (genomic position and alleles)
  my $infile = get_hash_option(\%FLAGS, "file");
  my $f_hgvs = get_hash_option(\%FLAGS, "field-hgvs");
  my $f_chr = $FLAGS{"field-chr"};
  my $fixed_chr = $FLAGS{chr};
  my $extract_hgvs_chr = $FLAGS{"extract-hgvs-chr"};
  my $extract_gpos = $FLAGS{"extract-gpos"};
  my $strip_trailing_comment = $FLAGS{"strip-trailing-comment"};
  my $patch_deletion_sequence = $FLAGS{"patch-deletion-sequence"};

  my $fai;
  if ($VERIFY_GENOMIC or $patch_deletion_sequence) {
    $fai = new FAI("-genome" => $FLAGS{genome} || die "-genome");
  }


  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  my $outfile = basename($infile) . ".hgvs_split.tab";

  my @extra;
  push @extra, $F_OUT_CHR if $extract_hgvs_chr;
  push @extra, (
	       $F_OUT_POS,
	       $F_OUT_RA,
	       $F_OUT_VA,
	       $F_OUT_TYPE
	      );

  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => \@extra,
  			      "-auto_qc" => 1,
			     );

  my $nsp = new NucleotideSubstitutionParser();

  if ($patch_deletion_sequence) {
    die "need -field-chr or -chr" unless $f_chr or $fixed_chr;
  }

  my $count_broken = 0;
  my $count_unknown = 0;
  my $count_ok = 0;

  my %unknown = map {$_, 1} (
			     "g.?",
			     # IARC
			     "-",
			     # NHGRI
			    );

  my $hgvs_chr;

  while (my $row = $df->get_hash()) {
    my $hgvs = $row->{$f_hgvs};
    my $gpos = $hgvs;
    if ($gpos) {
      if ($extract_gpos) {
	# e.g. IARC TP53 HGVS
	if ($gpos =~ /(g\..*)/) {
	  $gpos = $1;
	} elsif ($gpos =~ /g;(.*)/) {
	  # hack for IARC TP53 somatic: g;7577595_7577598del4
	  $gpos = "g." . $1;
	} else {
	  die "can't extract gpos from $gpos";
	}
      } elsif ($gpos =~ /:/) {
	die "possible 2-field HGVS-style annotation found in $gpos, use -extract-gpos";
      }

      if ($extract_hgvs_chr) {
	my @f = split /:/, $hgvs;
	die unless @f == 2;
	$hgvs_chr = $f[0];
	$row->{$F_OUT_CHR} = $hgvs_chr;
      }

    }

    $gpos =~ s/\([^\(]+\)$// if $strip_trailing_comment;
    # IARC TP53: g.7577157_7577498del342(del intron7)


    my $sj_pos = "";
    my $sj_ref_allele = "";
    my $sj_alt_allele = "";
    my $sj_type = "";
    my $vm;
    my $is_exception;
    my $is_unknown;

    if (not($gpos) or $unknown{$gpos}) {
      $is_unknown = 1;
    } else {
      $nsp->parse($gpos) || dump_die($row, "can't parse \"$gpos\"");
      $sj_pos = $nsp->start();
      if ($nsp->is_substitution()) {
	$sj_type = "sub";
	$sj_ref_allele = $nsp->reference_sequence();
	$sj_alt_allele = $nsp->variant_sequence();
      } elsif ($nsp->is_deletion()) {
	$sj_type = "del";
	$sj_ref_allele = $nsp->reference_sequence();
	unless ($sj_ref_allele) {
	  # hack: might not have sequence present, e.g.
	  # g.41258548_41267744delinsA
	  my $length = $nsp->end - $nsp->start + 1;
	  $sj_ref_allele = "N" x $length;
	}
	die "broken ref allele in $gpos: $sj_ref_allele" unless $sj_ref_allele =~ /^[acgtn]+$/i;
	$sj_alt_allele = "-";
      } elsif ($nsp->is_insertion()) {
	# simple insertion
	$sj_type = "ins";
	$sj_ref_allele = "-";
	$sj_alt_allele = $nsp->variant_sequence();
	if ($sj_alt_allele) {
	  die $sj_alt_allele unless $sj_alt_allele =~ /^[acgtn]+$/i;
	} else {
	  printf STDERR "ERROR: insertion with no variant sequence: %s\n", $gpos;
	  $is_exception = 1;
	}
      } elsif ($nsp->is_complex_indel()) {
	$sj_type = "complex";
	$sj_ref_allele = $nsp->reference_sequence();
	$sj_alt_allele = $nsp->variant_sequence();
	unless ($sj_ref_allele) {
	  # hack: might not have sequence present, e.g.
	  # g.41258548_41267744delinsA
	  my $length = $nsp->end - $nsp->start + 1;
	  $sj_ref_allele = "N" x $length;
	}
	if ($sj_alt_allele) {
	  die $sj_alt_allele unless $sj_alt_allele =~ /^[acgtn]+$/i;
	} else {
	  printf STDERR "ERROR: complex with no alt allele: %s\n", $gpos;
	  $is_exception = 1;
	}
      } else {
	die sprintf "FIX ME: gpos $gpos unhandled";
      }
    }

    if ($patch_deletion_sequence and $sj_ref_allele =~ /^n+$/i) {
      my $len = length($sj_ref_allele);
      my $chr = $fixed_chr || $hgvs_chr || $row->{$f_chr} || die;
      my $chunk = $fai->get_chunk(
				  "-id" => $chr,
				  "-start" => $sj_pos,
				  "-length" => $len
				 );
      printf STDERR "patched %s %d %s => %s\n", $gpos, $sj_pos, $sj_ref_allele, $chunk;
      $sj_ref_allele = $chunk;
    }

    if ($FLAGS{rc}) {
      # reverse complement broken data, e.g. Tierry TP53 functional
      # data includes transcript-oriented alleles even in HGVS g.*
      # format, which is supposed to use genome-oriented alleles (+)
      $sj_ref_allele = reverse_complement($sj_ref_allele);
      $sj_alt_allele = reverse_complement($sj_alt_allele);
    }

    if ($VERIFY_GENOMIC and $sj_ref_allele ne "-") {
      my $chunk = $fai->get_chunk(
				  "-id" => ($hgvs_chr || die "no chr"),
				  "-start" => $sj_pos,
				  "-length" => length($sj_ref_allele)
				 );
      if (uc($chunk) ne uc($sj_ref_allele)) {
	$is_exception = 1;
	printf STDERR "reference sanity check failed for %s %s %s, actual=%s\n", $hgvs_chr, $sj_pos, $sj_ref_allele, $chunk;
      }
    }

    if ($is_exception) {
      $count_broken++;
      $sj_pos = $sj_ref_allele = $sj_alt_allele = $sj_type = "";
    } elsif ($is_unknown) {
      $count_unknown++;
      $sj_pos = $sj_ref_allele = $sj_alt_allele = $sj_type = "";
    } else {
      $count_ok++;
    }

    $row->{$F_OUT_POS} = $sj_pos;
    $row->{$F_OUT_RA} = $sj_ref_allele;
    $row->{$F_OUT_VA} = $sj_alt_allele;
    $row->{$F_OUT_TYPE} = $sj_type;

    $rpt->end_row($row);
  }
  $rpt->finish();

  printf STDERR "OK:%d unknown:%d exception:%d\n", $count_ok, $count_unknown, $count_broken;


}

