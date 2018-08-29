#!/usr/bin/perl 
#migrate_db.pl

my $vcf = shift @ARGV;
my $outfile = $vcf;
$outfile =~ s/temp.vcf/union.vcf/;
open VCF, "<$vcf" or die $!;
open OUT, ">$outfile" or die $!;

while (my $line = <VCF>) {
  chomp($line);
  if ($line =~ m/#/) {
    print OUT $line,"\n";
    next;
  }
  my ($chrom, $pos,$id,$ref,$alt,$score,
      $filter,$annot,$format,@gts) = split(/\t/, $line);
  my %hash = ();
  foreach $a (split(/;/,$annot)) {
    my ($key,$val) = split(/=/,$a);
    $hash{$key} = $val;
  }
  my @deschead = split(/:/,$format);
  my @newgts;
  $newformat = 'GT:DP:AD:AO:RO';
  foreach my $idx (0..1) {
    my @gtinfo = split(/:/,$gts[$idx]);
    my %gtdata;
    foreach my $i (0..$#deschead) {
	$gtdata{$deschead[$i]} = $gtinfo[$i];
    }
    if ($i) {
	$gtdata{DP} = $gtdata{DDP} if ($gtdata{DDP});
	if (exists $gtdata{DAC}) {
	    $gtdata{AO} = $gtdata{DAC};
	    $gtdata{RO} = $gtdata{DP} - $gtdata{AO};
	}
    }else {
	$gtdata{DP} = $gtdata{NDP} if ($gtdata{NDP});
	if (exists $gtdata{NAC}) {
	    $gtdata{AO} = $gtdata{NAC};
	    $gtdata{RO} = $gtdata{DP} - $gtdata{AO};
      }
    }
    unless ($gtdata{AO}) {
	if ($gtdata{AD}){
	    ($gtdata{RO},$gtdata{AO}) = split(/,/,$gtdata{AD});
	}
    }
    unless ($gtdata{AD}) {
	if (exists $gtdata{RO} && exists $gtdata{AO}) {
	    $gtdata{AD} = join(",",$gtdata{RO},$gtdata{AO});
	}
    }
    if ($gtdata{DP} < $gtdata{AO}+$gtdata{RO}) {
	$gtdata{DP} = $gtdata{AO}+$gtdata{RO};
    }
    push @newgts, join(":",$gtdata{GT},$gtdata{DP},$gtdata{AD},$gtdata{AO},$gtdata{RO});
  }
  print OUT join("\t",$chrom,$pos,$id,$ref,$alt,$score,$filter,$annot,$newformat,@newgts),"\n";
}
