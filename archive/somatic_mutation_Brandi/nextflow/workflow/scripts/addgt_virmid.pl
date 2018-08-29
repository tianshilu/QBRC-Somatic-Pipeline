#!/usr/bin/perl 
#migrate_db.pl

my $vcf = shift @ARGV;
my $out = $vcf;
$out =~ s/\.vcf/.gt.vcf/g;

open VCF, "<$vcf" or die $!;
open OUT, ">$out" or die $!;
while (my $line = <VCF>) {
    chomp($line);
    $line =~ s/ID:/ID=/g;
    if ($line =~ m/#CHROM/) {
	print OUT join("\t",$line,'FORMAT','NORMAL','TUMOR'),"\n";
    }elsif ($line =~ m/#/) {
	print OUT $line,"\n";
    }else {
	my ($chrom, $pos,$id,$ref,$alt,$score,
	    $filter,$annot) = split(/\t/, $line);
	print OUT join("\t",$line,'GT','0/0','0/1'),"\n";
    }
}
