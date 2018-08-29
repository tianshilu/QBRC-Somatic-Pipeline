#!/usr/bin/perl -w
#add_readct_shimmer.pl

open CT, "<shimmer/som_counts.bh.txt" or die $!;
while (my $line = <CT>) {
  chomp($line);
  my ($chrom,$pos,$refnorm,$reftum,$refnormct,
      $reftumct,$alt,$altnormct,$alttumct,
      $pval,$qval) = split(/\t/,$line);
  $stats{$chrom}{$pos}{'tumor'} = {ao=>$alttumct,ro=>$reftumct,
				   ad=>join(",",$reftumct,$alttumct)};
  $stats{$chrom}{$pos}{'normal'} = {ao=>$altnormct,ro=>$refnormct,
				   ad=>join(",",$refnormct,$altnormct)};
}

open VCF, "<shimmer/somatic_diffs.vcf" or die $!;
open OUT, ">shimmer/somatic_diffs.readct.vcf" or die $!;
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
  my @samples = ('normal','tumor');
  my @newgts;
  foreach my $i (0..1) {
    %gtinfo = %{$stats{$chrom}{$pos}{$samples[$i]}};
    my @gtinfo = split(/:/,$gts[$i]);
    foreach my $i (0..$#deschead) {
      $gtinfo{$deschead[$i]} = $gtinfo[$i];
    }
    push @newgts, join(":",$gtinfo{GT},$gtinfo{DP},$gtinfo{ad},$gtinfo{ao},$gtinfo{ro});
  }
  $newformat = 'GT:DP:AD:AO:RO';
  print OUT join("\t",$chrom,$pos,$id,$ref,$alt,$score,$filter,$annot,$newformat,@newgts),"\n";
}
