#!/usr/bin/perl -w
#count_venn.vcf

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'fasta|f=s','help|h','prefix|p=s','complex|c=s','execdir|e=s');
my @zipvcf = @ARGV;

unless($opt{fasta}) {
    $opt{fasta} = '/project/shared/bicf_workflow_ref/GRCh38/hs38DH.fa';
}
unless($opt{prefix}) {
    $opt{prefix} = 'baysic';
}
my $total = 0;
open FASTA, "<$opt{fasta}\.fai" or die "unable to find fai index for $opt{fasta}, please index with samtools";
while (my $line = <FASTA>) {
    chomp($line);
    my ($chr,$length,$offset,$linelen,$linebytes) = split(/\t/,$line);
    $total += $length;
}

my %ct = ();
my %snpbins = ();
open CTS, ">baysic.cts" or die $!;
print CTS "Estimated sum: ".$total,"\n";
open VC, "<vcf_compare.out" or die $!;

while (my $line = <VC>) {
	my @key = ();
    next if($line =~ m/#/);
	if ($line =~ m/^VN\s+(.+)/) {
		my ($ct,@flist) = split(/\s+/,$1);
		foreach my $dset (@zipvcf) {
			if (grep(/$dset/,@flist)) {
				push @key, 1;
			}else {
		push @key, 0;
			}
		}
	$key = join("",@key);
	$ct{$key} = $ct;
	$total -= $ct;
    }
}
my @g = bits(scalar(@zipvcf));
$ct{$g[0]} = $total;
foreach (@g) {
    $ct{$_} = 0 unless ($ct{$_});
    print CTS join("\t",$_,$ct{$_}),"\n";
}

system("Rscript $opt{execdir}\/scripts/lca.R -c baysic.cts -s baysic.stats");

my @key1 = split(//,$g[-1]);
my @key2;
foreach $i (0..$#key1) {
    push @key2, $i if ($key1[$i] > 0);
}
open STATS, "baysic.stats" or die $!;
my @keeppos;
while (my $line = <STATS>) {
    chomp($line);
    next unless ($line =~ m/postprobs\[(\d+)\]\s+(\S+)/);
    my ($index,$prob) = ($1,$2);
    next unless ($prob >= 0.8);
    my $key = $g[$index-1];
    @key1 = split(//,$key);
    @key2 = ();
    foreach $i (0..$#key1) {
	push @key2, $i if ($key1[$i] > 0);
    }
    my $subset =  'integrate'.join('_',@key2).'.vcf.gz';
    push @keeppos, $subset if (-e $subset);
}
push @keeppos, $opt{complex} if ($opt{complex});
system("vcf-concat ".join(" ",@keeppos)." |vcf-sort |bgzip > final.integrated.vcf.gz");

sub bits { glob join "", map "{0,1}", 1..$_[0] }
