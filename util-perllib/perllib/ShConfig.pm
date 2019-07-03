#!/usr/bin/env perl
# Helper functions for using config files in shell script format
package ShConfig;

use strict;
use warnings;
use File::Path qw/ make_path /;

# Helper function for reading config
# - $file: config file to read (in sh format)
# Returns hash ref of var => val pairs
sub readConfig {
  my($file) = @_;
  return {} unless(-r $file);
  open(IN, $file) or die "Could not open $file";
  my $out = {};
  while(my $line = <IN>) {
    chomp $line;
    next unless($line);
    next if(substr($line, 0, 1) eq '#');
    my($var, $val) = split(/=/, $line, 2);
    $out->{$var} = $val if($var);
  }
  close(IN);
  return $out;
}

# Helper function for writing config in sh format
# - $file: config file to read (in sh format)
# - $config: configuration as a hash ref of var => val pairs
# - $desc: optional description to write to the file
sub writeConfig {
  my($file, $config, $desc) = @_;
  open(OUT, ">$file") or die "Could not open $file";
  print OUT "##SIP $desc\n" if($desc);
  while(my ($var, $val) = each (%$config) ) {
    print OUT "$var=$val\n";
  }
  close(OUT);
}

1;