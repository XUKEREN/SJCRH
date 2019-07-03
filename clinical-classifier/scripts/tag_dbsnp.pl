#!/bin/env perl
# add annotations from dbSNP
#
# Hi Mike,
#
# Would you be able to help me figure out why these three snps were not picked up by Bambino?
#
# rs66650371
# rs7775698
# rs94941442
#
# I’m working in
# /rgs01/project_space/weissgrp/SCD/partners/CMPB/SGP1000/BAMBINO_MEDAL
#
# And I grepped for these in the *germline.out.liftover.tab files.

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use MiscUtils qw(dump_die build_argv_list get_hash_option log_message);
use FileUtils qw(read_simple_file);
use DelimitedFile;
use Reporter;
use TdtConfig;
use TabixFile;
use Variant;
use TabixBatchAnnotation;
use VariantParsePolicy;
use DelimitedFileHP;

my $TABIX_INDEL_EQUIVALENCE_DISTANCE = 25;
#my $TABIX_BATCH_SIZE_DBNSFP = 200;
use constant TABIX_BATCH_SIZE_DBSNP => 1000;

my $QUEUE_SIZE = 10000;
#my $QUEUE_SIZE = 500;

my $FIELD_DBSNP = "dbSNP";
my $FIELD_DBSNP_HITS = "_dbsnp_hits";

my %FLAGS;
my @clopts = (
	      "-file=s",
	      "-files=s",

	      "-genome=s",
	      "-dbsnp=s",

	      "-dbsnp-only",
	      # filter output to only variants found in dbSNP

	      "-rs-to-variant=s",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts,
	   VariantParsePolicy::get_command_line_parameters()
	  );

if ($FLAGS{"rs-to-variant"}) {
  rs2variant();
  exit(0);
}

my $VPP = new VariantParsePolicy("-flags" => \%FLAGS);

my $infiles = build_argv_list(
			      "-flags" => \%FLAGS,
			      "-single" => "file",
			      "-set" => "files"
			     );

my $genome = $FLAGS{genome} || die "-genome";
my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";

#    add_batch_dbsnp(
#		    "-rows" => $rows
#		   );

unless ($FLAGS{dbsnp}) {
  $FLAGS{dbsnp} = config2tabix($config_genome, "DBSNP_TABIX_DIR");
}

foreach my $infile (@{$infiles}) {
  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			    );

  my $outfile = basename($infile) . ".dbsnp.tab";
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [ $FIELD_DBSNP ],
			      "-clobber" => 1,
#			      "-auto_qc" => 1,
			     );
  my @queue;
  my $flush = sub {
    if (@queue) {
      add_batch_dbsnp(
		      "-rows" => \@queue,
		     );
      foreach my $row (@queue) {
	my $usable = 1;
	if ($FLAGS{"dbsnp-only"}) {
	  $usable = 0 unless $row->{$FIELD_DBSNP};
#	  dump_die($row, "no dbSNP hit") unless $usable;
	}
	$rpt->end_row($row) if $usable;
      }
    }
    @queue = ();
  };

  while (my $row = $df->get_hash()) {
    push @queue, $row;
    &$flush() if @queue >= $QUEUE_SIZE;
  }
  &$flush();

  $rpt->finish();
}

sub add_batch_dbsnp {
  my (%options) = @_;
  my $rows = $options{"-rows"} || die "-rows";
  my $start_time = time;
  log_message("batch dbSNP start");
  my $tabix = get_tabix_dbsnp();
  my $f_row_key = "_user_row";
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $twobit = $config_genome->{TWOBIT} || die "no TWOBIT";

  my @query_variants;
  foreach my $r (@{$rows}) {
    my $v = $VPP->get_variant_from_row($r);
    $v->{$f_row_key} = $r;
    push @query_variants, $v;
  }

  my $tba = new TabixBatchAnnotation(
				     "-tabix" => $tabix,
				     "-twobit" => $twobit,
				     "-split_count" => TABIX_BATCH_SIZE_DBSNP,
				     "-f_tabix_chr" => "Chr",
				     "-f_tabix_pos" => "WU_HG19_Pos",
				     "-f_tabix_ref_allele" => "ReferenceAllele",
				     "-f_tabix_var_allele" => "MutantAllele",

				     "-user_row_key" => $f_row_key,
				     "-store_hits" => $FIELD_DBSNP_HITS,
#				     "-annotation_map" => \%map,
				    );

  $tba->query(
	      "-query" => \@query_variants,
	     );

  foreach my $r (@{$rows}) {
    my @dbsnp;
    if (my $hl = $r->{$FIELD_DBSNP_HITS}) {
      @dbsnp = map {$_->{rs} || dump_die($_)} @{$hl};
    }
    $r->{$FIELD_DBSNP} = join ",", @dbsnp;
  }


  log_message(sprintf "batch dbSNP annotation for %d rows took %d",
	      scalar(@{$rows}), time - $start_time);


}


sub config2tabix {
  # find tabix file from a genome config directory (i.e. tartan output)
  my ($config_genome, $cv) = @_;
  my $tf;
  if (my $dir = $config_genome->{$cv}) {
    if (-d $dir) {
      my @gz = glob($dir . "/*.gz");
      if (@gz == 1) {
	($tf) = @gz;
	my $tbi = $tf . ".tbi";
	die "where is $tbi" unless -s $tbi;
      } else {
	confess "ERROR: not exactly one .gz file in $dir for $cv";
      }
    }
  }
  return $tf;
}

sub get_tabix_dbsnp {
  my $f_tabix = $FLAGS{dbsnp} || die;
  printf STDERR "dbSNP: %s\n", $f_tabix;
  return get_batch_tabix($f_tabix);
}

sub get_batch_tabix {
  my ($f_tabix) = @_;
  return new TabixFile(
		       "-file" => $f_tabix,
		       "-indel_wiggle_bases" => $TABIX_INDEL_EQUIVALENCE_DISTANCE
		      );
}

sub rs2variant {
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";

  #    add_batch_dbsnp(
  #		    "-rows" => $rows
  #		   );

  unless ($FLAGS{dbsnp}) {
    $FLAGS{dbsnp} = config2tabix($config_genome, "DBSNP_TABIX_DIR");
  }

  my $f_dbsnp = $FLAGS{dbsnp} || die;
  printf STDERR "dbsnp: %s\n", $f_dbsnp;

  my $infile = $FLAGS{"rs-to-variant"} || die "-rs-to-variant";
  my $rsids = read_simple_file($infile);
  die "no IDs" unless @{$rsids};
  my %wanted;
  foreach my $rs (@{$rsids}) {
    die $rs unless $rs =~ /^rs\d+$/;
    $wanted{$rs} = 1;
  }

  my $rpt = new Reporter(
			 "-file" => basename($infile) . ".rs.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   rs
					   found
					   chr
					   pos
					   ref
					   alt
					)
				      ],
			 "-auto_qc" => 1,
			);


  my $df = new DelimitedFileHP(
			       "-file" => $f_dbsnp,
			      );
  $df->prepare_query("-fields" => ["rs"]);
  my $count = 0;
  my $found = 0;
  my $wanted_count = scalar keys %wanted;
  while ($df->next_row()) {
    my $rs = $df->get_query()->[0];
    if ($wanted{$rs}) {
      my %r;
      $r{rs} = $rs;
      $r{found} = 1;
      $r{chr} = $df->get_value("#Chr");
      $r{pos} = $df->get_value("WU_HG19_Pos");
      $r{ref} = $df->get_value("ReferenceAllele");
      $r{alt} = $df->get_value("MutantAllele");
      $rpt->end_row(\%r);

      delete $wanted{$rs};
      $found++;
      last unless %wanted;
      # quit if found everything
    }

    if (++$count % 1000000 == 0) {
      printf "%d, found %d/%d, at %s.%s\n", $count, $found, $wanted_count, $df->get_value("#Chr"), $df->get_value("WU_HG19_Pos");
    }

  }

  if (%wanted) {
    foreach my $rs (sort keys %wanted) {
      my %r;
      $r{rs} = $rs;
      $r{found} = 0;
      $r{chr} = $r{pos} = $r{ref} = $r{alt} = "";
      $rpt->end_row(\%r);
    }
  }

  $rpt->finish();

}
