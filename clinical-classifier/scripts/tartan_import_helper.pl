#!/bin/env perl
# helper for updating clinical flatfiles via tartan
# MNE 2/2015

use strict;
use warnings;

use 5.10.1;

use Getopt::Long;
use File::Spec;
use File::Basename;
use File::Copy;
use File::Path;
use Cwd;
#use Archive::Zip qw(:ERROR_CODES :CONSTANTS);
# not available on clinical cluster  :/

use MiscUtils qw(dump_die);
use FileUtils qw(read_simple_file write_simple_file);
use TdtConfig;
use TARTANUtils qw(tartan_genome_put_helper);

my %FLAGS;

GetOptions(
	   \%FLAGS,
	   "-genome=s",

	   "-deploy-to=s",
	   # generate commands to import files from dev to
	   # some other config for -param parameter (e.g. "prod")
	   # ...unfortunately this no longer seems to work, tartan
	   # seem to croak if input is a tartan run in a different environment

	   # import a single file to an existing config parameter:
	   "-param=s",
	   "-params=s",

	   "-find-param-name=s",
	   "-find-param-value=s",

	   "-single-import",
	   "-new-file=s",
	   # new file to import

	   "-tartan",

	   "-export",
	   # export config data and generate commands to import
	   "-genome-target=s",

	   "-create-ok",

	   "-put=s",
	   # put an existing tartan output directory into index

	   "-clone=s",
	   # copy contents of specified dir into a new ad-hoc run,
	   # then put it
	   "-diff",

	   "-dir-zip",
	   "-dir-link",
	  );

my $genome = $FLAGS{"genome"} || die "-genome";

if ($FLAGS{"deploy-to"}) {
  generate_release_commands();
  # deploy an entry from dev to another config (e.g. prod)
} elsif ($FLAGS{"export"}) {
  export_config();
} elsif ($FLAGS{"single-import"}) {
  import_single_file();
} elsif ($FLAGS{put}) {
  put_helper();
} elsif ($FLAGS{clone}) {
  clone_helper();
} else {
  die "specify -single-import or -deploy-to [release]\n";
}

sub import_single_file {
  my $new_file = $FLAGS{"new-file"};

  unless ($new_file) {
    my @files = glob("*");
    die "need exactly one file" unless @files == 1;
    $new_file = $files[0];
    printf STDERR "file: %s\n", $new_file;
  }

  die "new file $new_file does not exist" unless -e $new_file;

  $new_file = File::Spec->rel2abs($new_file);

  my $ROOT_KEY = "GENOME_ROOT";

  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $genome_root = $config_genome->{$ROOT_KEY} || die "no $ROOT_KEY in $genome config";

  if ($FLAGS{tartan}) {
    # use tartan index
    $genome_root = get_genome_root() || die;
  }

  my $config_file = get_wc_config();

  my $param = find_param_name();

  unless ($param) {
    printf "parameter?: ";
    $param = <STDIN>;
    chomp $param;
  }

  open(IN, $config_file) || die;
  my $where;
  while (<IN>) {
    next if /^#/ or not /\w/;
    my @f = split /\t/, $_;
    if ($f[0] eq $param) {
      $where = $f[1];
      last;
    }
  }
  die "can't find entry for $param in $config_file" unless $where;

  $where =~ s/^\$$ROOT_KEY/$genome_root/ || die;
  $where =~ s/\/$//;
  # make sure target is a symlink rather than a dir

  printf "import:\n  %s\nto:\n  %s\n", $new_file, $where;

  if ($FLAGS{"diff"} and -f $new_file) {
    chomp $where;
    chomp $new_file;
    my $cmd = sprintf 'diff %s %s', $where, $new_file;
    printf "%s\n", $cmd;
    system $cmd;
    print "\n";
  }

  printf STDERR "confirm? [y/N]: ";
  my $resp = <STDIN>;

  chomp $resp;
  if ($resp and lc($resp) eq "y") {
    my $cmd = sprintf 'tartan import %s %s', $new_file, $where;
    printf STDERR "running: %s\n", $cmd;
    system $cmd;
    die "$cmd failed" if $?;
  } else {
    print STDERR "quitting.\n";
  }
}

sub find_param_name {
  my %options = @_;
  my $param = $FLAGS{param};
  my $find_param_name = $FLAGS{"find-param-name"};
  my $find_param_value = $FLAGS{"find-param-value"};

  my $config_file = get_wc_config();
  open(IN, $config_file) || die;
  my $where;
  my %maybe;
  while (<IN>) {
    next if /^#/ or not /\w/;
    my @f = split /\t/, $_;
    if ($param and $f[0] eq $param) {
      $maybe{$f[0]} = $f[1];
    } elsif ($find_param_name and $f[0] =~ /$find_param_name/i) {
      $maybe{$f[0]} = $f[1];
    } elsif ($find_param_value and $f[1] =~ /$find_param_value/i) {
      $maybe{$f[0]} = $f[1];
    }
  }
  close IN;

  unless ($param) {
    if (scalar keys %maybe == 1) {
      ($param) = (keys %maybe);
      printf STDERR "found: param %s (%s)\n", $param, $maybe{$param};
    } elsif (%maybe) {
      dump_die(\%maybe, "ambiguous");
    } else {
      die "no param found, config file is $config_file\n";
    }
  }

  if ($options{"-value"}) {
    return $maybe{$param};
  } else {
    return $param;
  }
}

sub generate_release_commands {
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  die "current env must be dev" unless $ENV{TARTAN_ROOT} =~ /\/dev\//;

  my $param = find_param_name() || die "no param";
  my $current = $config_genome->{$param} || die "no data for $param";
  # this is the RESOLVED config, loses symlink from standardized
  # config dir entry
  my $target_env = $FLAGS{"deploy-to"} || die;

  my $config_file = get_wc_config();
  open(IN, $config_file) || die;
  my $where;
  my %maybe;
  my $find_param_name = $param;
  while (<IN>) {
    chomp;
    next if /^#/ or not /\w/;
    my @f = split /\t/, $_;

    if ($param and $f[0] eq $param) {
      $maybe{$f[0]} = $f[1];
    } elsif ($find_param_name and $f[0] eq $find_param_name) {
      $maybe{$f[0]} = $f[1];
    }
  }
  close IN;

  my $param_value;
  if (scalar keys %maybe == 1) {
    ($param) = (keys %maybe);
    $param_value = $maybe{$param};
  } elsif (%maybe) {
    dump_die(\%maybe, "ambiguous");
  } else {
    die "no param found, config file is $config_file\n";
  }

  my $target_root = get_genome_root(
			   "-env" => $target_env,
			  );

#  my $target_root = sprintf '/rgs01/resgen/%s/tartan/index/reference/Homo_sapiens/%s', $target_env, $genome;
  # HACK: only works on research!
  die "where is $target_root" unless -d $target_root;
  # any way to derive programmatically?

  my $target_file = $param_value;
  $target_file =~ s/\$GENOME_ROOT/$target_root/ || die "fail";

  unless (-l $target_file or -d $target_file) {
    unless ($FLAGS{"create-ok"}) {
      die "$target_file is not a link already, specify -create-ok to proceed\n";
    }
  }

  if (1 or -f $current) {
    my $cmd = sprintf 'tartan import %s %s',
      $current, $target_file;
    printf "%s\n", $cmd;
    # command only executable in target environment
  } elsif (-d $current) {
    printf "tartan_import_helper.pl -genome %s -param %s -clone %s\n", $genome, $param, $current;
  } else {
    die;
  }

}

sub get_wc_config {
  my $base = sprintf '%s/wc/cluster_code/trunk', $ENV{HOME};

  my $config_file = sprintf '%s/configs/data/genome/%s.config.txt',
    $base, $genome;
  die "where is $config_file" unless -s $config_file;
  return $config_file;
}

sub export_config {
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  die "current env must be dev" unless $ENV{TARTAN_ROOT} =~ /\/dev\//;

  my @params;
  if (my $p = $FLAGS{param}) {
    push @params, $p;
  } elsif (my $f = $FLAGS{params}) {
    my $set = read_simple_file($f);
    @params = @{$set};
  } else {
    die "specify -param PARAM or -params FILE\n";
  }

  my $start_dir = getcwd();
  my $genome_target = $FLAGS{"genome-target"} || $genome;

  foreach my $param (@params) {
    chdir($start_dir) || die;
    my $current = $config_genome->{$param} || die "no data for $param";

    my $subdir = $param;
    mkpath($subdir) unless -d $subdir;
    chdir($start_dir . "/" . $subdir) || die;

    my $tartan_cmd;
    if (-f $current) {
      # simple file
      my $local = basename($current);
      copy($current, $local) || die;
#      $tartan_cmd = sprintf "tartan_import_helper.pl -genome %s -param %s -single-import -new-file %s -tartan \$*\n", $genome_target, $param, $local;
      $tartan_cmd = sprintf "tartan_import_helper.pl -genome %s -param %s -single-import -new-file %s \$*\n", $genome_target, $param, $local;
    } elsif ($FLAGS{"dir-zip"}) {
      # directory: zip contents
      my @files = glob($current . "/*");
      die unless @files;
      my $zipfile = $param . ".zip";
      unlink $zipfile;

      my $cmd = sprintf 'zip -j0 %s %s', $zipfile, join " ", @files;
      system $cmd;
      die $cmd if $?;
      die unless -s $zipfile;

      # not available on clinical cluster  :/
      #
      # my $zip = Archive::Zip->new();
      # foreach my $f (@files) {
      # 	my $zm = $zip->addFile($f);
      # 	$zm->desiredCompressionMethod(COMPRESSION_STORED);
      # }
      # unless ($zip->writeToFileNamed($zipfile) == AZ_OK) {
      # 	die 'write error';
      # }

      $tartan_cmd = sprintf 'tartan_import_helper.pl -genome %s -param %s -put $*' . "\n", $genome_target, $param;
    } elsif ($FLAGS{"dir-link"}) {
      my $local = basename($current);
      $tartan_cmd = sprintf 'tartan_import_helper.pl -genome %s -param %s -single -new %s', $genome_target, $param, $local;
      unlink $local;
      symlink($current, $local);
    } else {
      die "specify -dir-zip or -dir-link\n";
    }


    write_simple_file([$tartan_cmd], "import.sh");

  }
}

sub put_helper {
  my $tartan_out = $FLAGS{put} || die;

  my $path = find_param_name("-value" => 1);
  my $subdir = $path;
  $subdir =~ s/^\$\w+\/// || die;

  tartan_genome_put_helper(
			   "-out" => $tartan_out,
			   "-genome" => $genome,
			   "-subdir" => $subdir
			  );
}

sub clone_helper {
  my $in_dir = $FLAGS{clone} || die;
  # tartan output run being cloned (e.g. from rdev to rprod)
  die unless -d $in_dir;

  my @infiles = glob($in_dir . "/*");

  $in_dir =~ /ad_hoc\/(\w+)\-/ || die "$in_dir not ad-hoc?";
  my $prefix = $1 || die;
  my $cmd_new = sprintf 'tartan_run.pl new type=ad_hoc name=%s prefix=1', $prefix;
  # create an ad-hoc run with the same prefix

  my $run_dir = `$cmd_new`;
  # bleah; learn TartanRun module?
  printf STDERR "run dir: %s\n", $run_dir;
  chomp $run_dir;

  system "tartan_run.pl $run_dir startWork";

  my $out_dir = join "/", $run_dir, "output";
  die "where is $out_dir" unless -d $out_dir;

  my @tbi;
  foreach my $infile (@infiles) {
    my $outfile = sprintf '%s/%s', $out_dir, basename($infile);
    printf STDERR "copy %s => %s\n", $infile, $outfile;
    if (-d $infile) {
      my $cmd = sprintf 'cp -rp %s %s', $infile, $outfile;
      system($cmd);
      die "$cmd exited with $?" if $?;
    } else {
      copy($infile, $outfile) || die "can't copy $infile to $outfile";
      push @tbi, $outfile if $outfile =~ /\.tbi/;
    }
  }

  foreach my $tbi (@tbi) {
    # ensure timestamp is newer for tabix index files
    printf STDERR "touch %s\n", $tbi;
    system "touch $tbi";
  }

  system "tartan_run.pl $run_dir endWork";

  $FLAGS{put} = $run_dir . "/output";
  put_helper();

}

sub get_genome_root {
  my %options = @_;

  my $tr = $ENV{TARTAN_ROOT} || die;
  #    $genome_root = sprintf '%s/index/reference/Homo_sapiens/%s', $tr, $genome;
  # FAIL: GRCh38 on rdev is deployed to GRCh38_no_alt
  my $key;
  my $wc_root;
  if ($tr eq "/research/rgs01/resgen/dev/tartan") {
    $key = "rdev";
    if (my $alt = $options{"-env"}) {
      $key = "r" . $alt;
    }
    $wc_root = sprintf '/research/rgs01/resgen/dev/wc/%s', getpwuid($<) || die;
  } elsif ($tr eq "/cgs01/clingen/dev/tartan") {
    $key = "cdev";
    if (my $alt = $options{"-env"}) {
      $key = "c" . $alt;
    }

    $wc_root = sprintf '/cgs01/clingen/dev/wc/%s', getpwuid($<) || die;
  } else {
    die "unhandled root $tr";
  }


  die $wc_root unless $wc_root;
  my $genome_root;
  my $f_config = sprintf '%s/cluster_code/trunk/configs/data/genome/%s/%s.config.txt', $wc_root, $key, $genome;
  die "where is $f_config" unless -s $f_config;
  open(CFG, $f_config) || die;
  while (<CFG>) {
    chomp;
    my @f = split /\t/, $_;
    if ($f[0] eq "GENOME_ROOT") {
      $genome_root = $f[1];
      last;
    }
  }
  return $genome_root;
}

