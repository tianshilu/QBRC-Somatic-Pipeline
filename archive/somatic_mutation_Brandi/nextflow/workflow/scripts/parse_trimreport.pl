#!/usr/bin/perl -w
#parse_trimreport.pl

my $outfile = shift @ARGV;
my @trimreport = @ARGV;
chomp(@trimreport);

open OUT, ">$outfile" or die $!;
foreach $treport (@trimreport) {
    open RP, "<$treport" or die $!;
    while (my $line = <RP>) {
	chomp($line);
	if ($line =~ m/^Total reads processed:\s+(\S+)/) {
	    $total_reads = $1;
	    $total_reads =~ s/,//g;
	} elsif ($line =~ m/^Reads written \(passing filters\):\s+(\S+)/) {
	    $pass = $1;
	    $pass =~ s/,//g;
	} elsif ($line =~ m/^Reads with adapters:\s+(\S+)/) {
	    $withadapter = $1;
	    $withadapter =~ s/,//g;
	} elsif ($line =~ m/^Quality-trimmed:\s+\S+\s+bp\s+\((\S+)\%/) {
	    $percentqualtrim = $1;
	    $percentqualtrim =~ s/,//g;
	}
       
    }
    print OUT join("\t",$treport,$total_reads,$pass,
		   $withadapter,$percentqualtrim),"\n";
}
close OUT;
