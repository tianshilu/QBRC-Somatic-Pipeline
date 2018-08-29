#!/usr/bin/perl -w
#uploadqc.pl

open OUT, ">sequence.stats.txt" or die $!;
print OUT join("\t",'Sample','total.raw','total.trimmed','pairs','maprate',
	       'propair','ontarget','frac.dups','library.size','medinsert',
	       'avginsert','stdinsert','perc.10x','perc.20x','perc.50x',
	       'perc.100x','perc.200x','perc.500x'),"\n";

my @statfiles = @ARGV;

foreach $sfile (@statfiles) {
    $sfile =~ m/(\S+)\.libcomplex.txt/;
    my $prefix = $1;
    my %hash;
    open FLAG, "<$prefix\.trimreport.txt" or die $!;
    while (my $line = <FLAG>) {
	chomp($line);
	my ($file,$raw,$trim) = split(/\t/,$line);
	$hash{rawct} += $raw;
	$hash{trimct} += $trim;
    }
    open FLAG, "<$prefix\.flagstat.txt" or die $!;
    while (my $line = <FLAG>) {
	chomp($line);
	if ($line =~ m/(\d+) \+ \d+ in total/) {
	    $hash{total} = $1;
	    $hash{total} = $hash{trimct} if ($hash{trimct} && $hash{trimct} > $hash{total});
	}elsif ($line =~ m/(\d+) \+ \d+ read1/) {
	    $hash{pairs} = $1;
	}elsif ($line =~ m/(\d+) \+ \d+ mapped/) {
	    $hash{maprate} = 100*sprintf("%.4f",$1/$hash{total});
	}elsif ($line =~ m/(\d+) \+ \d+ properly paired/) {
	    $hash{propair} = 100*sprintf("%.4f",$1/$hash{total});
	}elsif ($line =~ m/(\d+) \+ \d+ properly paired/) {
	    $hash{propair} = 100*sprintf("%.4f",$1/$hash{total});
	}
    }
    open FLAG, "<$prefix\.ontarget.flagstat.txt" or die $!;
    while (my $line = <FLAG>) {
	chomp($line);
	if ($line =~ m/(\d+) \+ \d+ in total/) {
	    $hash{ontarget} = $1;
	}elsif ($line =~ m/(\d+) \+ \d+ read1/) {
	    $hash{ontarget_pairs} = $1;
	}elsif ($line =~ m/(\d+) \+ \d+ properly paired/) {
	    $hash{ontarget_propair} = 100*sprintf("%.4f",$1/$hash{total});
	}elsif ($line =~ m/(\d+) \+ \d+ properly paired/) {
	    $hash{ontarget_propair} = 100*sprintf("%.4f",$1/$hash{total});
	}
    }
    open DUP, "<$prefix\.libcomplex.txt" or die $!;
    while (my $line = <DUP>) {
	chomp($line);
	if ($line =~ m/## METRICS/) {
	    $header = <DUP>;
	    $nums = <DUP>;
	    chomp($header);
	    chomp($nums);
	    my @stats = split(/\t/,$header);
	    my @row = split(/\t/,$nums);
	    my %info;
	    foreach my $i (0..$#stats) {
		$info{$stats[$i]} = $row[$i];
	    }
	    $hash{percdups} = sprintf("%.4f",$info{PERCENT_DUPLICATION});
	    $hash{libsize} = sprintf("%.4f",$info{ESTIMATED_LIBRARY_SIZE});
	}
    }
    $hash{medinsert} = 0;
    $hash{avginsert} = 0;
    $hash{stdinsert} = 0;
    open DUP, "<$prefix\.hist.txt";
    while (my $line = <DUP>) {
	chomp($line);
	if ($line =~ m/## METRICS/) {
	    $header = <DUP>;
	    $nums = <DUP>;
	    chomp($header);
	    chomp($nums);
	    my @stats = split(/\t/,$header);
	    my @row = split(/\t/,$nums);
	    my %info;
	    foreach my $i (0..$#stats) {
		$info{$stats[$i]} = $row[$i];
	    }
	    $hash{medinsert} = sprintf("%.0f",$info{MEDIAN_INSERT_SIZE});
	    $hash{avginsert} = sprintf("%.0f",$info{MEAN_INSERT_SIZE});
	    $hash{stdinsert} = sprintf("%.0f",$info{STANDARD_DEVIATION});
	}
    }
    my %cov;
    open COV, "<$prefix\.genomecov.txt" or die $!;
    while (my $line = <COV>) {
      chomp($line);
      my ($all,$depth,$bp,$total,$percent) = split(/\t/,$line);
      $cov{$depth} = $percent;
    }
    my @depths = sort {$a <=> $b} keys %cov;
    my @perc = @cov{@depths};
    my @cum_sum = cumsum(@perc);
    $hash{'perc10x'} = 100*sprintf("%.4f",1-$cum_sum[10]);
    $hash{'perc20x'} = 100*sprintf("%.4f",1-$cum_sum[20]);
    $hash{'perc50x'} = 100*sprintf("%.4f",1-$cum_sum[50]);
    $hash{'perc100x'} = 100*sprintf("%.4f",1-$cum_sum[100]);
    $hash{'perc200x'} = 100*sprintf("%.4f",1-$cum_sum[200]);
    $hash{'perc500x'} = 100*sprintf("%.4f",1-$cum_sum[500]);
    print OUT join("\t",$prefix,$hash{rawct}, $hash{total},$hash{pairs},$hash{maprate},$hash{propair},
		   $hash{ontarget},$hash{percdups},$hash{libsize},$hash{medinsert},$hash{avginsert},
		   $hash{stdinsert},$hash{'perc10x'},$hash{'perc20x'},$hash{'perc50x'},
		   $hash{'perc100x'},$hash{'perc200x'},$hash{'perc500x'}),"\n";
  }

sub cumsum {
  my @nums = @_;
  my @cumsum = ();
  my $mid = 0;
  for my $i (0..$#nums) {
    $mid += $nums[$i];
    push(@cumsum,$mid);
  }
  return @cumsum;
}
