# create reverse complement of a fastq.gz file
#!/usr/bin/perl
use strict;
use warnings;

my ($input_R1,$output_R2)=@ARGV;
my ($c,$cc,$line);

open(FILE_IN,"zcat ".$input_R1." |");
open(FILE_OUT,"| gzip > ".$output_R2);

while ($line=<FILE_IN>)
{
  # read name
  print FILE_OUT $line;
 
  # sequence
  $line=<FILE_IN>;
  $line=~s/\n//;
  $line=uc(reverse($line));
  for $c (split //, $line)
  {
    if ($c eq "A") {
      $cc="T";
    }elsif ($c eq "T") {
      $cc="A";
    }elsif ($c eq "C") {
      $cc="G";
    }elsif ($c eq "G") {
      $cc="C";
    }else {
      $cc="N";
    }
    print FILE_OUT $cc;
  }
  print FILE_OUT "\n";

  # read name 2
  $line=<FILE_IN>;
  print FILE_OUT $line;

  # quality score
  $line=<FILE_IN>;
  $line=~s/\n//;
  print FILE_OUT reverse($line)."\n";  
}

close(FILE_IN);
close(FILE_OUT);

