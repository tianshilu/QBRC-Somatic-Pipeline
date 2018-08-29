#!/usr/bin/perl

# Usage: <STDIN> | prog | <STOUT>

use strict;
use warnings;

my %prog;
my %cnt;

while (<>) {
  chomp;
  my @e = split /\t/;
  my $k = join("\t", @e[1 .. 6]);
  push @{$prog{$k}}, $e[0];
  push @{$cnt{$k}}, join(',', @e[7 .. $#e]);
}
foreach my $s (keys %prog) {
  print join("\t", $s, join(';', @{$prog{$s}}), join(';', @{$cnt{$s}})), "\n";
}
