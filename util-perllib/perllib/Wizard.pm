#!/usr/bin/env perl
# Helper functions for wizard-like scripts
package Wizard;

use strict;
use warnings;
use File::Path qw/ make_path /;
use TdtConfig; 

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(prompt promptYN promptQueue systemVerbose);

# Helper function for reading config
# - $file: config file to read (in tab-delimited text format)
# Returns hash ref of var => val pairs
sub readConfig {
  return TdtConfig::readConfig("", @_);
}

# Helper function for writing sip config
# - $file: config file to read (in tab-delimited text format)
# - $config: configuration as a hash ref of var => val pairs
# - $desc: optional description to write to the file
sub writeConfig {
  return TdtConfig::writeConfig(@_);
}

# Helper function for prompting
# - $msg: prompt message to print
# - $default: default value if user just hits enter
# Returns whatever the user typed in (chomped) or default if nothing
sub prompt {
  my($msg, $default) = @_;
  $default = "" unless($default);
  print "$msg [$default]> ";
  chomp(my $line = <STDIN>);
  return $line ? $line : $default;
}

# Helper function for Y/N prompting
# - $msg: prompt message to print
# - $default: default value if user just hits enter (should be y or n)
# Returns y or n (lowercased)
sub promptYN {
  my($msg, $default) = @_;
  $msg = "$msg (Y/N)?";
  $default = $default ? substr(uc($default), 0, 1) : "";
  my $resp;
  do {
    $resp = substr(lc(prompt($msg, $default)), 0, 1);
  } while($resp ne "n" && $resp ne "y");
  return $resp;
}

# Helper function for queue prompting
# Returns whatever the user typed in (chomped) or default if nothing
sub promptQueue {
  return prompt("Enter the queue to use for job submission", $ENV{'AFC_DEFAULT_QUEUE'});
}


# Run a command, telling the user what command you are running
# - $desc: description of the command
# - $cmd: command itself
# Returns commands exitcode
sub systemVerbose {
  my($desc, $cmd) = @_;
  print "Running $desc command: $cmd\n";
  system $cmd;
  return $?;
}

# Checks for available configs.
# - directory containing configs (will be created if it does not exist)
# Returns number of configs found
sub countConfigs {
  my ($configDir) = @_;
  make_path $configDir unless(-e $configDir);
  chomp(my $numConfigs = `ls $configDir | wc -l`);
  print "".($numConfigs ? "OK" : "WARNING" ).
  " Found $numConfigs configuration files in $configDir\n";
  return $numConfigs;
}

# Prompt for config file
# - directory containing configs
# Returns the path to the chosen config file, may be empty
sub promptConfig {
  my $sipConfigFilename;
  my ($configDir) = @_;
  print "Which config file would you like to use (hit enter for none):\n";
  system "ls -1 -rt $configDir | while read file; do echo \"\$file \`grep -m 1 '^##DESC' $configDir/\$file\`\"; done | sed 's/##DESC/:/'";
  print "\n";
  chomp(my $def = `ls -1 -rt $configDir | tail -n 1`);
  $sipConfigFilename = prompt("Enter config name (or leave blank for none)", $def);
}

# Determine editor
# - default editor if not set in environment (default default is vi)
# Returns the editor
sub getEditor {
  my ($default) = @_;
  my $editor = $ENV{"EDITOR"};
  if($editor) {
    print "Using editor $editor from the EDITOR environment variable.\n";
  }
  else {
    $editor = $default ? $default : "vi";
    print "Did not find editor in EDITOR environment variable.  Using $editor\n";
  }
  return $editor;
}

1;
