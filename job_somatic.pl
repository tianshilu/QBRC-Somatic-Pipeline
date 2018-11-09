###########  nucleus-specific sbatch wrapper for calling somatic mutations  #############
# for other job submission system, you should be able to easily change "sbatch" to the appropriate command
#
# The instructions must be followed exactly!!! 
# input format:
# jobs: the batch job design file, it has 6 columns separated by \t, the first four are fastq files or bam files (normal+tumor), 
#       the fifth is the output folder, and the last is "PDX" or "". Commented lines ("#" at the front) are skipped
# example: the demo job submission shell script. A default one is in this folder
# thread,build,index,java17: follow those in somatic.pl
# n: bundle $n somatic calling jobs into one submission
#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';

my ($jobs,$example,$thread,$build,$index,$java17,$n)=@ARGV;
my ($line,$line1,@items,$i,$job);

my $path=abs_path($0);
$path=~s/job_somatic\.pl//;

open(JOB,$jobs) or die "Cannot find the design file!\n";
$i=0;

while ($line=<JOB>)
{
  $line=~s/(\r|\n)//g;
  if ($line eq "" || $line=~/^#/) {next;}

  if ($i++ % $n==0)
  {
    $job="somatic_calling_".$i.".sh";
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
  print SCRIPT "perl ".$path."/somatic.pl ".$items[0]." ".$items[1]." ".$items[2]." ".$items[3].
    " ".$thread." ".$build." ".$index." ".$java17." ".$items[4]." ".$items[5]."\n";

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

#perl /home2/twang6/software/cancer/somatic/job_somatic.pl \
#design.txt \
#/project/bioinformatics/Xiao_lab/shared/neoantigen/code/somatic/example/example.sh \
#32 hg38 \
#/project/bioinformatics/Xiao_lab/shared/neoantigen/data/hg38/hs38d1.fa \
#/cm/shared/apps/java/oracle/jdk1.7.0_51/bin/java \
#2

