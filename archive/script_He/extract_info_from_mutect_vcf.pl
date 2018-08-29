#!/usr/bin/perl

# Usage: prog FILE_normal_tumor.vcf.gz | <STDOUT>

use strict;
use warnings;

my $file = shift;

$file =~ /([^\.|\/]+)\.([^\.]+).vcf/;
my $norm = $1;
my $tumor = $2;

open FILE, "gzip -cd $file |"
  || die "Cannot open $file: $!\n";
while (<FILE>) {
  next if /^#/;
  chomp;
  my @e = split /\t/;
  my $idx = -1;
  my @form = split /:/, $e[8];
  foreach my $i (0 .. $#form) {
    if ($form[$i] eq 'AD') {
      $idx = $i;
      last;
    }
  }
  my $n1 = 0;
  my $n2 = 0;
  my $t1 = 0;
  my $t2 = 0;
  if ($idx != -1) {
    ($t1, $t2) = split /,/, (split /:/, $e[9])[$idx];
    ($n1, $n2) = split /,/, (split /:/, $e[10])[$idx];
  }
  print join("\t", 'mutect', $norm, $tumor, $e[0], $e[1], $e[3], $e[4], $n1, $n2, $t1, $t2), "\n";
}
close FILE;
