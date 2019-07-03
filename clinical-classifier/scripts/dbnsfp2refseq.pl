#!/bin/env perl
# convert isoform lists in dbNSFP to RefSeq where possible, via UniProt, etc.
# MNE 1/2018
# EXPERIMENTAL/UNFINISHED
# TO DO:
# - convert to use UniProtIDMapping.pm

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list get_hash_option);
use FileUtils qw(universal_open);
use DelimitedFile;
use Reporter;
use TdtConfig;
use DelimitedFileHP;
use MiscUtils qw(float_trim);

use constant DBNSFP_NULL => ".";
#use constant DBNSFP_DELIM => ";";
my $DBNSFP_DELIM = ";";
# constant doesn't work within split// even with () for some reason

my $FILTER_NULL_TO_BLANK = 1;
my $FILTER_FLOAT = 1;

my $F_OUT_UNIPROT_REFSEQ = "Uniprot_acc_Polyphen2_refseq";
my $F_OUT_ENSEMBL_REFSEQ = "Ensembl_transcriptid_refseq";
my $F_OUT_MA_REFSEQ = "MutationAssessor_UniprotID_refseq";

my %FLAGS;
my @clopts = (
	      "-dbnsfp=s",
	      "-uniprot=s",
	      # ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
	      # *** TO DO: try to get the version used by specific dbNSFP build?
	      # dbNSFP 3.5a: Uniprot, release 2017_07
	      # => it seems archived releases don't save this
	      "-genome=s",

	      "-biomart=s",
	      # TSV from biomart query output

	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $genome = $FLAGS{genome} || die "-genome";

my $f_uniprot = $FLAGS{uniprot};
unless ($f_uniprot) {
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  $f_uniprot = $config_genome->{CLINCLS_UNIPROT_ID_MAPPING_GZIP_FILE} || die;
}

my $f_dbnsfp = $FLAGS{dbnsfp};
unless ($f_dbnsfp) {
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  $f_dbnsfp = glob(($config_genome->{DBNSFP3_TABIX_DIR} || die) . "/*gz");
}

my $f_biomart = $FLAGS{biomart} || die "-biomart";

my %ENS2REF;
my $df = new DelimitedFile("-file" => $f_biomart,
			   "-headers" => 1,
			  );
while (my $row = $df->get_hash()) {
  my $ens = $row->{"Transcript stable ID"} || die;
  my $nm = $row->{"RefSeq mRNA ID"};
  push @{$ENS2REF{$ens}}, $nm;
  # can be more than one
}


my %UNIPROT2NM;
my %UNIPROT_COUNTS;
my %UNIPROT2KBID;

printf STDERR "parsing %s...", $f_uniprot;
my $fh = universal_open($f_uniprot) || die;
my $saw_subisoforms;
while (<$fh>) {
  chomp;
  my ($uniprot, $type, $id) = split /\t/, $_;
  $saw_subisoforms = 1 if $uniprot =~ /\-\d+$/;

  $UNIPROT_COUNTS{$uniprot}++;
  # hack to help resolve hits to multiple uniprot entries, e.g.
  # "-chr" => 5,
  # "-pos" => 112173767,
  # "-reference-base" => "T",
  # "-variant-base" => "G",
  # "-nm" => "NM_000038",

  if ($type eq "RefSeq_NT") {
    # may map to multiple entries, e.g.
    # P31946  RefSeq_NT       NM_003404.4
    # P31946  RefSeq_NT       NM_139323.3
    #	printf STDERR "map uniprot %s => %s\n", $uniprot, $id;
    push @{$UNIPROT2NM{$uniprot}}, $id;
  } elsif ($type eq "UniProtKB-ID") {
    die "duplicate $uniprot $id" if $UNIPROT2KBID{$id};
    $UNIPROT2KBID{$id} = $uniprot;
#    print STDERR "$id => $uniprot\n";
  }
}

close $fh || die "i/o error $! $?";
die "id parsing error" unless scalar keys %UNIPROT2NM;
die "ERROR: didn't find any sub-isoform identifiers (required for refSeq mapping).  Newer UniProt ID mapping file needed?" unless $saw_subisoforms;
print STDERR "done\n";

my %extra;
$df = new DelimitedFileHP(
			     "-file" => $f_dbnsfp,
			     "-headers_extra" => [
						  $F_OUT_UNIPROT_REFSEQ,
						  $F_OUT_ENSEMBL_REFSEQ,
						  $F_OUT_MA_REFSEQ
						 ]
			    );

foreach my $r (qw(
		   ExAC
		   hg18
		   gnomAD
		   1000Gp3
		   ESP6500
		   clinvar_
		   Ensembl_
		   TWINSUK_
		   ALSPAC_

		   Uniprot_
		)) {
  $df->add_header_exclude_regexp($r);
}

foreach my $f (qw(
		   aaref
		   aaalt
rs_dbSNP150
genename
cds_strand
refcodon
codonpos
codon_degeneracy
Ancestral_allele
AltaiNeandertal
Denisova

		)) {
  $df->add_header_exclude_field($f);
}

#
#  remove main chr/pos fields for all genomes.  The genome-specific ones
#  we actually want to use will be renamed to use common
#
my $already_tabix = $df->has_header("#chr");

my $F_HG19_CHR = "hg19_chr";
my $F_HG19_POS = "hg19_pos(1-based)";
my $F_HG38_CHR = $already_tabix ? "#chr" : "chr";
my $F_HG38_POS = "pos(1-based)";
my $F_COOKED_CHR = "genome_chr";
my $F_COOKED_POS = "genome_pos";

my %rename;
if ($genome eq "GRCh37-lite") {
  # use hg19, remove hg38
  $rename{$F_HG19_CHR} = $F_COOKED_CHR;
  $rename{$F_HG19_POS} = $F_COOKED_POS;
  $df->add_header_exclude_field($F_HG38_CHR);
  $df->add_header_exclude_field($F_HG38_POS);
} elsif ($genome eq "GRCh38") {
  # use hg38, remove hg19
  $rename{$F_HG38_CHR} = $F_COOKED_CHR;
  $rename{$F_HG38_POS} = $F_COOKED_POS;
  $df->add_header_exclude_field($F_HG19_CHR);
  $df->add_header_exclude_field($F_HG19_POS);
}
$df->header_out_rename(\%rename);

$df->write_init();

while ($df->next_row()) {
  #
  # attempt to map Polyphen2 Uniprot accessions to refSeq IDs:
  #
  my @up_refseqs;
  my $list = get_nsfp_list($df, "Uniprot_acc_Polyphen2");
  foreach my $id (@{$list}) {
    my $match_nms;
    if ($id ne DBNSFP_NULL) {
      if ($match_nms = $UNIPROT2NM{$id}) {
	# exact match
      } else {
	my @try;
	if ($id =~ /\-\d+$/) {
	  # dbNSFP ID contains a sub-version number, try without.
	  my $stripped = $id;
	  $stripped =~ s/\-\d+$//;
	  push @try, $stripped;
	}

	foreach my $try (@try) {
	  if ($match_nms = $UNIPROT2NM{$try}) {
	    # e.g. dbNSFP Q96NU1-1 doesn't match, but Q96NU1 does
	    # (NM_152486.2).  Is this too aggressive?  Does it
	    # ever introduce ambiguity?
	    last;
	  }
	}
      }
    }

    my $refseq = $match_nms ? join(",", @{$match_nms}) : "";
    # there can be ambiguity so each entry may itself be a list

    push @up_refseqs, $refseq;
  }
  $extra{$F_OUT_UNIPROT_REFSEQ} = join ";", @up_refseqs;

  #
  #  convert MutationAssessor UniProt ID to refseq:
  #
  my $uid = $df->get_value("MutationAssessor_UniprotID");
  my $ma_refseq = "";
  if ($uid ne DBNSFP_NULL) {
    if (my $acc = $UNIPROT2KBID{$uid}) {
      if (my $nms = $UNIPROT2NM{$acc}) {
	$ma_refseq = join ",", @{$nms};
      } else {
	printf STDERR "ERROR: can't find uniprot refseq for %s\n", $acc;
      }
    } else {
      printf STDERR "ERROR: can't find uniprot accession for %s\n", $uid;
    }
  }
  $extra{$F_OUT_MA_REFSEQ} = $ma_refseq;

  #
  #  convert Ensembl_transcriptid to refseq:
  #
  my $list_enst = get_nsfp_list($df, "Ensembl_transcriptid");
  my @ens_refseqs;
  foreach my $tid (@{$list_enst}) {
    my $match_nms = $ENS2REF{$tid};
    unless ($match_nms) {
      printf STDERR "ERROR: no ENSEMBL entry for %s\n", $tid;
    }
    push @ens_refseqs, $match_nms ? join(",", @{$match_nms}) : "";
  }
  $extra{$F_OUT_ENSEMBL_REFSEQ} = join ";", @ens_refseqs;

  my $row = $df->current_line();

  if ($FILTER_NULL_TO_BLANK) {
    foreach (@{$row}) {
      $_ = "" if $_ eq DBNSFP_NULL;
    }
  }

  if ($FILTER_FLOAT) {
    foreach (@{$row}) {
      $_ = float_trim($_,
#		      "-precision" => 5
		     );
    }
  }


  $df->write_row("-extra" => \%extra);
}

$df->write_finish();


sub get_nsfp_list {
  my ($df, $key) = @_;
  my $v = $df->get_value($key);
#  print STDERR "raw=$v\n";
  return [ split /$DBNSFP_DELIM/, $v ];
}
