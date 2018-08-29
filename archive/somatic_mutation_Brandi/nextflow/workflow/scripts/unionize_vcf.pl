#!/usr/bin/perl -w
#integration_vcf_snpindel_highcov.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'refdir|r=s','gver|g=s','outdir|o=s');

$fasta = $opt{refdir}."/genome.fa";

my @vcffiles = @ARGV;
chomp(@vcffiles);
open SH, ">integrate.sh" or die $!;
print SH "module load gatk/3.5 python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 platypus/gcc/0.8.1  bcftools/intel/1.3 samtools/intel/1.3 jags/4.2.0 vcftools/0.1.14\n";
#system("module load gatk/3.5 python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 platypus/gcc/0.8.1  bcftools/intel/1.3 samtools/intel/1.3 jags/4.2.0");

foreach $vcf (@vcffiles) {
    ($prefix,$caller,$ext) = split(/\./,$vcf);
    system("zgrep CHROM $vcf > $vcf.chrline.txt");
    if ($caller eq 'sspanel') {
	$shuff = "ssvar.shuff.vcf.gz";
	print SH "ln -s $vcf $shuff\n";
	print SH "tabix $shuff\n"; 
	push @nonevcfs, $shuff;
	push @callers, $caller;
	push @combvars, "--variant:".$caller." ".$shuff;
    }elsif (-s "$vcf.chrline.txt") {
	my $shuff = "$caller\.shuff.vcf.gz";
	if ($caller eq 'varscan') {
	    print SH "tabix $vcf\n";
	    print SH "bcftools annotate -x 'INFO/SSC' $vcf |vcf-shuffle-cols -t ssvar.shuff.vcf.gz | vcf-sort | bgzip > $shuff\n";
	}else {
	    print SH "vcf-shuffle-cols -t ssvar.shuff.vcf.gz $vcf | vcf-sort | bgzip > $shuff\n";
	}
	print SH "tabix $shuff\n";
	push @nonevcfs, $shuff;
	push @callers, $caller;
	push @combvars, "--variant:".$caller." ".$shuff;
    }
}
print SH "java -Xmx32g -jar \$GATK_JAR -R $fasta -T CombineVariants --filteredrecordsmergetype KEEP_UNCONDITIONAL -genotypeMergeOptions PRIORITIZE ".join(" ",@combvars)." -priority ".join(",",@callers)." -o ".$prefix.".int.vcf\n";
print SH "bedtools multiinter -i ".join(" ",@nonevcfs)." -names ".join(" ",@callers)." |cut -f 1,2,3,5 | bedtools sort -i stdin | bedtools merge -c 4 -o distinct >  $prefix\.integrate.bed\n";
print SH "bgzip $prefix\.integrate.bed\n";
print SH "tabix $prefix\.integrate.bed.gz\n";
print SH "bcftools annotate -a $prefix\.integrate.bed.gz --columns CHROM,FROM,TO,CallSet -h $opt{refdir}/CallSet.header $prefix\.int.vcf > $prefix.temp.vcf\n";
close SH;
