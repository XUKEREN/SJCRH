#!/usr/bin/env perl
# Helper functions for working with databases
package DbUtil;

use strict;
use warnings;

use DBI;
use Term::ReadKey;
use TdtConfig;
use Carp;

# Connects to a database using information in a database connection file
# Config file should have the following values:
# * DSN (e.g. dbi:Pg:dbname=DB;host=HOST;port=PORT)
# * USER (user name, can be "SELF" to use whoami)
# * PASSFILE (a path to a file containing the db password)
# * PASS (can be "PROMPT" to prompt interactively)
#
# Note that these options are ALWAYS used no matter what you specify:
# {RaiseError => 1, AutoCommit => 0}
#
# Args:
# - a category to look in (should be 'dbconn')
# - the name of the config file to use
# - connect options hash ref
# - true if no prompting should be done (will fail if password is PROMPT)
sub connectUsingTdtConfig{
	my ($cat, $conf, $opts, $noPrompt) = @_;
	
	# Always RaiseError and don't AutoCommit
	# If other options are preferred, this class cannot be used. 
	$opts->{'RaiseError'} = 1; 
	$opts->{'AutoCommit'} = 0; 
	
	# Read config
	my $config = TdtConfig::readConfig($cat, $conf);
	my $dsn = $config->{'DSN'};
	
	# Handle special value for user
	my $user = $config->{'USER'};
	my $whoami = `whoami`;
	chomp $whoami;
	$user =~ s/SELF/$whoami/;
	#chomp($user = `whoami`) if($user eq 'SELF');
	
	# Handle special logic for password
	my $pass = undef;
	# Try file first
	if($config->{'PASSFILE'}) {
		# Note: abs_path doesn't work here (doesn't resolve ~ I think)
		chomp(my $passFile = `readlink -f $config->{'PASSFILE'}`);
		# FYI: Perl's stat returns the file permission in decimal. 0600 -> 33152
		if (-f "$passFile") {
			if ((stat($passFile))[2] eq "33152"){
				$pass = `cat $passFile`;
				chomp $pass;
			}
			else{
				die "password file: $passFile exists but permissions are not set to 0600"; 
			}
		}
	}
	# If that doesn't work, then use PASS, but handle special value "PROMPT"
	unless($pass) {
		if($config->{'PASS'} eq "PROMPT"){
			# Fail if noPrompt
			if($noPrompt) {
				warn "Cannot connect to DB without prompting when password is PROMPT";
				return undef;
			}
			
			# Otherwise, prompt without echo
			ReadMode 'noecho';
			print STDERR "Enter password for DSN $dsn, user $user: ";
			chomp($pass = <STDIN>);
			ReadMode 0;
			print "\n";
		}
		else {
			$pass = $config->{'PASS'};
		}
	}

	# Connect
	my $conn = DBI->connect($dsn, $user, $pass, $opts);
	
	# Handle init command
	$conn->do($config->{'INIT_COMMAND'}) if($config->{'INIT_COMMAND'});
	
	return $conn;
}

=pod
=item C<loadHashFromStatement>
Loads a hash using values from a select query

The prepared statement must be for a select with two fields, where the first
is the key, and the second is the value

See also executeTwoColumnsOrDieHashref()

Arguments
- hash reference to load
- prepared statement to run
=cut
sub loadHashFromStatement {
	my ($cache, $sth) = @_;
	$sth->execute();
	while(my($key, $id) = $sth->fetchrow()) {
		$cache->{$key} = $id;
	}
	$sth->finish;
}

=pod
=item C<executeOrConfess>
Wrapper around execute that will confess if there is a problem.

Regular execute will either return undef, or will die without stack trace (if
RaiseError).  This one will either succeed or confess

All parameters are passed to execute, and return is value from execute.
=cut
sub executeOrConfess {
	my ($sth, @params) = @_;
	my $out = eval {
		return $sth->execute(@params);
	};
	if($@ || ! $out) {
		confess "Internal error running query: $@";
	}
	return $out;
}

=pod
=item C<executeOneRowOrDieHashref>
Executes a query, and returns the single row, asserting tha there is only one
row.

If the die message is false, then undef will be returned instead of an exception
being thrown if there are no results. 

Parameters:
1. Statement handle
2. Die message if no results
3... Statement parameters

Returns: the one row as a hash ref
=cut
sub executeOneRowOrDieHashref {
	my ($sth, $msg, @params) = @_;
	
	# Execute, die with internal error on failure
	DbUtil::executeOrConfess($sth, @params);
	
	# Fetch first row, die with custom message if none
	my $row = $sth->fetchrow_hashref();
	unless($row) {
		confess $msg if($msg);
		return undef;
	}
	
	# If you can fetch another, then die with internal error
	confess "Internal error (got too many rows)" if($sth->fetchrow_array);
	
	# Finish and return
	$sth->finish();
	return $row;
}

=pod
=item C<executeOneRowOrDieArray>
Executes a query, and returns the single row, asserting tha there is only one
row.

If the die message is false, then undef will be returned instead of an exception
being thrown if there are no results. 

Parameters:
1. Statement handle
2. Die message if no results
3... Statement parameters

Returns: the one row as an array
=cut
sub executeOneRowOrDieArray {
	my ($sth, $msg, @params) = @_;
	# Execute, die with internal error on failure
	DbUtil::executeOrConfess($sth, @params);
	
	# Fetch first row, die with custom message if none
	my @row = $sth->fetchrow_array();
	unless(@row) {
		confess $msg if($msg);
		return undef;
	}
	
	# If you can fetch another, then die with internal error
	confess "Internal error (got too many rows)" if($sth->fetchrow_array);
	
	# Finish and return
	$sth->finish();
	return @row;
}

=pod
=item C<executeOneValueOrDie>
Executes a query, and returns the single value, asserting that there is only one
row.

This will return the value from the first field in the query.

If the die message is false, then undef will be returned instead of an exception
being thrown if there are no results. 

Parameters:
1. Statement handle
2. Die message if no results
3... Statement parameters

Returns: the one value as a scalar
=cut
sub executeOneValueOrDie {
	my @row = executeOneRowOrDieArray(@_);
	return @row ? $row[0] : undef;
}

=pod
=item C<executeOneColumnOrDieArrayRef>
Executes a query, and returns all results from the first output column in an
array ref.

Parameters:
1. Statement handle
2... Statement parameters

Returns: all values from the first column as an array ref
=cut
sub executeOneColumnOrDieArrayRef {
	my ($sth, @params) = @_;
	
	# Execute, die with internal error on failure
	DbUtil::executeOrConfess($sth, @params);
	# Fetch the first column of each row
	my @out = ();
	while(my $row = $sth->fetchrow_arrayref()) {
		push @out, $row->[0];
	}
	
	# Finish and return
	$sth->finish();
	return \@out;
}

=pod
=item C<executeTwoColumnsOrDieHashRef>
Executes a query, and returns all results from the first two output columns in a
hash ref where the keys come from the first column and values from the second.

Parameters:
1. Statement handle
2... Statement parameters

Returns: all values from the first column as an array ref
=cut
sub executeTwoColumnsOrDieHashRef {
	my ($sth, @params) = @_;
	
	# Execute, die with internal error on failure
	DbUtil::executeOrConfess($sth, @params);
	
	# Fetch the first two columns of each row
	my %out = ();
	while(my $row = $sth->fetchrow_arrayref()) {
		$out{$row->[0]} = $row->[1];
	}
	
	# Finish and return
	$sth->finish();
	return \%out;
}

=pod

=item C<queryForPkey>

Executes a query and returns the primary key value, using convention by default
to determine the primary key.

Named parameters:
- dbh: Database handle
- table: table name
- pkey: primary key field name (default is <table>_id)
- lenient: boolean, if true, returns undef if not found; if false or not
    specified, then it confesses
- fields: array ref of field names (mut ex with 'data')
- values: array ref of values (in same order as fields)
- data: hash ref containing field => value pairs (mut ex with fields/values)

=cut
sub queryForPkey {
	my %args = @_;
	
	my $table = $args{'table'};
	my $pkey = $args{'pkey'} ? $args{'pkey'} : $table.'_id';
	my (@fields, @values);
	if($args{'data'}) {
		@fields = keys(%{$args{'data'}});
		@values = values(%{$args{'data'}});
	}
	else {
		@fields = @{$args{'fields'}};
		@values = @{$args{'values'}};
	}

	my $sth = $args{'dbh'}->prepare_cached(
		"SELECT $pkey FROM $table WHERE ".join(' = ? AND ', @fields)." = ?");
	return executeOneValueOrDie(
		$sth,
		$args{'lenient'} ?
			undef :
			"Could not find $pkey for (".join(',', @fields).") = (".
			join(',', @values).")",
		@values);
}

=pod

=item C<genericInsert>

Generic implementation of single-record INSERT.  Uses last_insert_id to attempt
to return the pkey just inserted (may not always work).

Named parameters:
- dbh: database handle
- table: table into which to insert
- pkey: primary key field name (default is <table>_id if id value given,
    otherwise not used)
- id: primary key value
- return_id: pass true to return pkey value
- fields: array ref of field names (mut ex with 'data')
- values: array ref of values (in same order as fields)
- data: hash ref containing field => value pairs (mut ex with fields/values)

=cut
sub genericInsert {
	my %args = @_;
	
	my $table = $args{'table'};
	my (@fields, @values);
	if($args{'data'}) {
		@fields = keys(%{$args{'data'}});
		@values = values(%{$args{'data'}});
	}
	else {
		@fields = @{$args{'fields'}};
		@values = @{$args{'values'}};
	}
	if($args{'id'}) {
		my $pkey = $args{'pkey'} ? $args{'pkey'} : $table.'_id';
		unshift @fields, $pkey;
		unshift @values, $args{'id'};
	}
	
	# Prepare
	my $dbh = $args{'dbh'};
	my $sth = $dbh->prepare_cached(
		"INSERT INTO $table(".join(', ', @fields).")\n".
		"  VALUES(".join(', ', ('?') x scalar(@fields)).")");
	
	# Execute
	push @values, (undef) x (@fields - @values) if(@fields > @values);
	$sth->execute(@values);
	
	return unless($args{'return_id'});
	return $args{'id'} if($args{'id'});
	return $dbh->last_insert_id(undef, undef, $table, undef);
}

1;
