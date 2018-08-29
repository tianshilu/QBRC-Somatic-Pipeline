#!/usr/bin/perl

# Usage: prog FILE_normal_tumor.vcf(or vcf.gz) | <STDOUT>

use strict;
use warnings;

my $file = shift;

$file =~ /([^\.|\/]+)\.([^\.|\/]+)\.[^\/]*vcf/;
my $norm = $1;
my $tumor = $2;

if ($file =~ /.gz$/) {
  open FILE, "gzip -cd $file |"
    || die "Cannot open $file: $!\n";
}
else {
  open FILE, "<$file"
    || die "Cannot open $file: $!\n";
}

while (<FILE>) {
  next if /^#/;
  chomp;
  my @e = split /\t/;
  my $idxref = -1;
  my $idxalt = -1;
  my @form = split /:/, $e[8];
  foreach my $i (0 .. $#form) {
    if ($form[$i] eq $e[3] . 'U') {
      $idxref = $i;
    }
    elsif ($form[$i] eq $e[4] . 'U') {
      $idxalt = $i;
    }
  }
  my $n1 = 0;
  my $n2 = 0;
  my $t1 = 0;
  my $t2 = 0;
  if ($idxref != -1) {
    $n1 = (split /,/, (split /:/, $e[9])[$idxref])[0];
    $n2 = (split /,/, (split /:/, $e[9])[$idxalt])[0];
  }
  if ($idxalt != -1) {
    $t1 = (split /,/, (split /:/, $e[10])[$idxref])[0];
    $t2 = (split /,/, (split /:/, $e[10])[$idxalt])[0];
  }
  print join("\t", 'strelka', $norm, $tumor, $e[0], $e[1], $e[3], $e[4], $n1, $n2, $t1, $t2), "\n";
}
close FILE;
