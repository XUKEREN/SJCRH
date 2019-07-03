package SimpleVariantDB;
# simple automatic SQL database, originally intended for variants
# (genomic and protein)
# MNE 11/2017
# TO DO:
# - required properties:
#   - whether db should be auto-lifted or not (e.g. if native in both genomes)
# - auto-generate and save MD5

use strict;
use Exporter;
use Cwd;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use DelimitedFile;
use DBTools qw(selectall_hashref);
use Carp qw(confess);

use constant F_VDB_CHR => "svdb_chr";
use constant F_VDB_POS => "svdb_pos";
use constant F_VDB_RA => "svdb_ra";
use constant F_VDB_VA => "svdb_va";
use constant F_VDB_GENE => "svdb_gene";
use constant F_VDB_AA => "svdb_aachange";
use constant F_VDB_CDS => "svdb_cds";
use constant F_VDB_TRANSCRIPT => "svdb_transcript";
use constant F_VDB_PMID => "svdb_pmid";

use constant TABLE_CONSTANTS => "svdb_constant_fields";
use constant TABLE_PROPERTIES => "svdb_properties";

my @ARGV_RAW;
BEGIN {
  @ARGV_RAW = @main::ARGV;
}

# fields should be unique-ish so as to co-exist with raw row data

@SimpleVariantDB::ISA = qw(Configurable Exporter);
@SimpleVariantDB::EXPORT_OK = qw(
F_VDB_CHR
F_VDB_POS
F_VDB_RA
F_VDB_VA
F_VDB_GENE
F_VDB_AA
F_VDB_CDS
F_VDB_TRANSCRIPT
F_VDB_PMID
);
# TO DO: maybe break these into a separate package?
# most of this code is generic, will work with any format

use MethodMaker qw(
	dbi
	blankify
trim_flanking_whitespace
strip_non_ascii
unquote
convert_european_decimal
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->blankify({});
  $self->trim_flanking_whitespace(1);
  $self->strip_non_ascii(1);
  $self->unquote(1);
  $self->convert_european_decimal(1);
  $self->configure(%options);
  return $self;
}

sub init_table {
  # create table for database
  my ($self, %options) = @_;
  my $name = $options{"-name"} || die "-name";
  my $dbi = $self->dbi || die;
  my $rows = $options{"-rows"};
  my $blankify = $self->blankify();
  my $strip_non_ascii = $self->strip_non_ascii();
  my $unquote = $self->unquote();
  my $map = $options{"-map"};
  my $map_only = $options{"-map-only"};

  if (my $infile = $options{"-file"}) {
    die "-file requires -map" unless $map;
    $rows = [];
    my $df = new DelimitedFile("-file" => $infile,
			       "-headers" => 1,
			      );
    while (my $row = $df->get_hash()) {
      if ($map_only) {
	# only keep the mapped columns
	my %row;
	$self->map_row(
		       "-row-in" => $row,
		       "-map" => $map,
		       "-row-out" => \%row,
		      );
	push @{$rows}, \%row;
      } else {
	die "supplemental map not implemented yet";
      }
    }
  }

  if ($self->table_exists($name)) {
    my $cmd = sprintf "drop table %s", $name;
    printf STDERR "%s\n", $cmd;
    $dbi->do($cmd);
  }

  if ($options{"-auto-configure"}) {
    # automatically detect data types
    die "-auto-configure requires -rows" unless $rows;
    my @columns = sort keys %{$rows->[0]};
    my %types;
    foreach my $column (@columns) {
      $types{$column} = $self->detect_sql_type(
					       "-rows" => $rows,
					       "-column" => $column,
					      );
    }

    # create table:
    my $sql = sprintf "CREATE TABLE %s (\n", $name;
    $sql .= sprintf "  %s\n", join ",\n  ", map {"$_ " . $types{$_}} @columns;
    $sql .= ")\n";
    $dbi->do($sql) || die "exec failed for $sql";

    printf STDERR "%s\n", $sql;
  } else {
    die "manual config";
  }

  if ($options{"-load"}) {

    $self->set_constants(
		     "-name" => $name,
		     "-values" => ($options{"-constants"} || die "-constants")
		    );

    my $properties = $options{"-properties"} || die "-properties";
    $properties->{"cwd"} = getcwd();
    $properties->{"load_date"} = scalar localtime;
    $properties->{"command_line"} = $0 . " " . join " ", @ARGV_RAW;
    # may not always be sufficient but likely close enough

    $self->set_properties(
		     "-name" => $name,
		     "-values" => $properties,
		    );

    my @columns = sort keys %{$rows->[0]};

    my $sql = sprintf 'insert into %s (%s) values (%s)',
      $name,
	join(",", @columns),
	  join(",", split //, "?" x @columns);

    my $sth = $dbi->prepare($sql) || die "can't prepare $sql";

    printf STDERR "inserting %d rows...", scalar @{$rows};
    my $ping = scalar @{$rows} / 50;
    my $inserted = 0;
    foreach my $row (@{$rows}) {
      my @v = map {$row->{$_}} @columns;

      if (%{$blankify}) {
	foreach (@v) {
	  $_ = "" if $blankify->{$_};
	}
      }

      foreach (@v) {
	$_ = $self->clean_value($_);
      }
      
      $sth->execute(@v) || die "can't insert";
      print STDERR "." if ++$inserted % $ping == 0;
    }
    print STDERR "\n";
  }
}

sub get_table_names {
  my ($self) = @_;
  my $dbi = $self->dbi || die;
  my $sth = $dbi->table_info;
  my %tables;
  while (my $row = $sth->fetchrow_hashref) {
    my $table = $row->{TABLE_NAME};
    $tables{$table} = 1 if $table;
  }
  $sth->finish();
  return \%tables;
}

sub table_exists {
  my ($self, $table) = @_;
  my $tables = $self->get_table_names();
  return $tables->{$table};
}

sub map_row {
  my ($self, %options) = @_;
  my $row_in = $options{"-row-in"} || die "-row-in";
  my $row_out = $options{"-row-out"} || $row_in;
  my $map = $options{"-map"} || die "-map";
  my $whitespace_trim = $self->trim_flanking_whitespace();
  foreach my $f_to (keys %{$map}) {
    my $f_from = $map->{$f_to};
    dump_die($row_in, "no field $f_from") unless exists $row_in->{$f_from};
    my $v = $row_in->{$f_from};
    $v = "" unless defined $v;
    if ($whitespace_trim and $v =~ /\S/) {
      $v =~ s/^\s+//;
      $v =~ s/\s+$//;
    }

#    printf STDERR "%s => %s\n", $f_to, $v;
    $row_out->{$f_to} = $v;
  }
}

sub detect_sql_type {
  # TO DO:
  #  - float
  #  - signed numeric values
  #  - char
  my ($self, %options) = @_;
  my $rows = $options{"-rows"} || die;
  my $col = $options{"-column"} || die;

  my $max_length = 0;
  my $is_int = 1;
  foreach my $row (@{$rows}) {
    dump_die($row, "where is $col") unless defined $row->{$col};
    my $v = $row->{$col};
    if (defined $v) {
      $v = $self->clean_value($v);
      my $len = length $v;
      $max_length = $len if $len > $max_length;

      if (length($v) == 0 or $v =~ /^\d+$/) {
	# FIX ME: sign
	# FIX ME: float
      } else {
	printf STDERR "%s is not integer: %s\n", $col, $v unless $is_int == 0;
	$is_int = 0;
      }
    }
  }

  my $type;
  if ($is_int) {
    $type = "INT";
  } else {
    $type = sprintf 'VARCHAR2(%d)', $max_length;
  }
  return $type;
}

sub add_blankify {
  my ($self, $v) = @_;
  $self->blankify->{$v} = 1;
}

sub set_constants {
  # set values that are the same for all rows
  my ($self, %options) = @_;
  my $name = $options{"-name"} || die;
  my $values = $options{"-values"} || die;
  my $dbi = $self->dbi || die;
  my $sql;

  if ($self->table_exists(TABLE_CONSTANTS)) {
#    my $cmd = sprintf "drop table %s", TABLE_CONSTANTS;
#    printf STDERR "%s\n", $cmd;
#    $dbi->do($cmd);
# WRONG: don't want to trash entries for other databases
    $sql = sprintf "delete from %s where table_name='%s'",
    TABLE_CONSTANTS, $name;
    $dbi->do($sql) || die $sql;
  } else {
    $sql = sprintf "CREATE TABLE %s (\n", TABLE_CONSTANTS;
    $sql .= "  table_name VARCHAR2(100),\n";
    $sql .= "  field VARCHAR2(100),\n";
    $sql .= "  value VARCHAR2(1000)\n";
    $sql .= ")\n";
    $dbi->do($sql) || die "can't exec $sql";
  }

  $sql = sprintf 'insert into %s (%s) values (?,?,?)',
    TABLE_CONSTANTS, join ",", qw(table_name field value);
  my $sth = $dbi->prepare($sql) || die "can't prepare $sql";

  foreach my $k (keys %{$values}) {
    $sth->execute($name, $k, $values->{$k}) || die "can't insert";
  }
}

sub set_properties {
  # database-specific properties set at load time
  my ($self, %options) = @_;
  my $name = $options{"-name"} || die;
  my $values = $options{"-values"} || die;
  my $dbi = $self->dbi || die;
  my $sql;

  if ($self->table_exists(TABLE_PROPERTIES)) {
    $sql = sprintf "delete from %s where table_name='%s'",
    TABLE_PROPERTIES, $name;
    $dbi->do($sql) || die $sql;
  } else {
    $sql = sprintf "CREATE TABLE %s (\n", TABLE_PROPERTIES;
    $sql .= "  table_name VARCHAR2(100),\n";
    $sql .= "  field VARCHAR2(100),\n";
    $sql .= "  value VARCHAR2(1000)\n";
    $sql .= ")\n";
    $dbi->do($sql) || die "can't exec $sql";
  }

  $sql = sprintf 'insert into %s (%s) values (?,?,?)',
    TABLE_PROPERTIES, join ",", qw(table_name field value);
  my $sth = $dbi->prepare($sql) || die "can't prepare $sql";

  foreach my $k (keys %{$values}) {
    $sth->execute($name, $k, $values->{$k}) || die "can't insert";
  }
}

sub get_rows_with_constants {
  # get all rows from a table, populating any constant values
  my ($self, %options) = @_;
  my $name = $options{"-name"} || die;
  my $dbi = $self->dbi();

  my $rows = selectall_hashref($dbi, sprintf 'select * from %s', $name);

  my $r_const = selectall_hashref($dbi, sprintf 'select * from %s where table_name="%s"', TABLE_CONSTANTS, $name);
  my %constant;
  foreach my $r (@{$r_const}) {
    $constant{$r->{field}} = $r->{value};
  }
  if (%constant) {
    my @f = keys %constant;
    foreach my $r (@{$rows}) {
      @{$r}{@f} = map {$constant{$_}} @f;
    }
  }
  return $rows;
}

sub clean_value {
  my ($self, $v) = @_;
  $v =~ tr/\x80-\xFF//d if $self->strip_non_ascii();
  $v =~ s/^([\'\"])(.*)\1/$2/ if $self->unquote();
  # remove quotes around strings
  $v =~ s/^(\d+),(\d+)$/$1.$2/ if $self->convert_european_decimal();
  # convert e.g. 0,86701 => 0.86701 (TP53 functional data)

  return $v;
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
