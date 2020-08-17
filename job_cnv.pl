###########  sbatch wrapper for CNV  #############
# for other job submission system, you should be able to easily change "sbatch" to the appropriate command
# all paths must be absolute
#
# The instructions must be followed exactly!!! 
# input format:
# jobs: the batch job design file, it has 7 columns separated by \t. Commented lines ("#" at the front) are skipped
#       The first four are fastq files or bam files (normal+tumor), 
#       Next is the somatic mutation calling output file
#       Next is the output folder.
#       Last is "PDX" or "human" or "mouse"
# example: the demo job submission shell script. A default one is in this folder
# thread,index,disambiguate: follow those in cnv.pl
# n: bundle $n CNV calling jobs into one submission
#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';

my ($jobs,$example,$thread,$index,$disambiguate,$n)=@ARGV;
my ($line,$line1,@items,$i,$job);

my $path=abs_path($0);
$path=~s/job_cnv\.pl//;

open(JOB,$jobs) or die "Cannot find the design file!\n";
$i=0;

while ($line=<JOB>)
{
  $line=~s/(\r|\n)//g;
  if ($line eq "" || $line=~/^#/) {next;}

  if ($i++ % $n==0)
  {
    $job="cnv_".$i.".sh";
    open(SCRIPT,">".$job) or die "Cannot write to the shell submission script!\n";

    # write header
    open(HEADER,$example) or die "Cannot find the example shell script!\n";
    while ($line1=<HEADER>)
    {
      if ($line1=~/JOBSTART/) {last;}
      print SCRIPT $line1;
    }
    close(HEADER);
  }

  # write submission job
  @items=split("\t",$line);
  print SCRIPT "perl ".$path."/cnv.pl ".$items[0]." ".$items[1]." ".$items[2]." ".$items[3].
    " ".$thread." ".$index." ".$items[4]." ".$items[5]." ".$items[6]." ".$disambiguate."\n";

  if ($i % $n==0)
  {
    close(SCRIPT); 
    system("sbatch ".$job);
    unlink($job);
    sleep(1);
  }
}

close(JOB);

if ($i % $n!=0)
{
  close(SCRIPT);
  system("sbatch ".$job);
  unlink($job);
}

#perl /home2/twang6/software/cancer/somatic/job_cnv.pl \
#design.txt \
#~/example.sh 32 \
#/project/shared/xiao_wang/data/hg38/hs38d1.fa \
#/project/shared/xiao_wang/software/disambiguate_pipeline 2
