#!/usr/bin/perl

# Usage: <STDIN> | prog CHR POS REF ALT | <STDOUT>

use strict;
use warnings;

my $chr = shift;
my $pos = shift;
my $ref = shift;
my $alt = shift;

$chr--;
$pos--;
$ref--;
$alt--;

while (<>) {
  chomp;
  my @e = split;
  my $type = substr($e[$alt], 0, 1);
  my $len = length($e[$alt]);
  if ($type eq '+') {
    substr($e[$alt], 0, 1, $e[$ref]);
  }
  elsif ($type eq '-') {
    my $tmp = $e[$ref];
    substr($e[$alt], 0, 1, $e[$ref]);
    $e[$ref] = $e[$alt];
    $e[$alt] = $tmp;
  }
  print join("\t", @e), "\n";
}
