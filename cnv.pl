#####################  CNV (and also quality check) analysis for exome-seq data  ########################
#
# The instructions must be followed exactly!!! 
# prerequisite in path: 
#   R (R.oo and R.utils packages updated to the latest version), bwa, sambamba, perl (need Parallel::ForkManager), 
#   samtools (>=1.6), cnvkit, fastqc, disambiguate (use conda env by Yunguan, contact yunguan.wang@utsouthwestern.edu)
# input format: 
# fastq files: (1) fastq1 and fastq2 of normal sample, fastq1 and fastq2 of tumor sample
#              default input is full path to the 4 fastq files for tumor and normal samples 
#              (2) if need to directly input bam files, use "bam path_to_bam.bam" in replace of the two corresponding fastq input files
#              can be a mixture of fastq and bam input
# thread: number of threads to use. Recommended: 32
# index: path (including file name) to the human or mouse reference genome fasta file
# somatic: somatic mutation calling output file. This is for adjusting CNV by somatic mutation VAF. Set to 1 to turn off this adjustment
# output: the output folder, it will be deleted (if pre-existing) and re-created during analysis
# pdx: "PDX" or "human" or "mouse". if "PDX", can only handle paired-end sequencing reads
# disambiguate: disambiguate path, prepared by Yunguan, /project/shared/xiao_wang/software/disambiguate_pipeline
# need at least 128GB of memory
#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';
use File::Copy;
use Parallel::ForkManager;

# hidden paramters
my $read_ct_max=50000000; 

my ($fastq1_normal,$fastq2_normal,$fastq1_tumor,$fastq2_tumor,$thread,$index,$somatic,$output,
  $pdx,$disambiguate)=@ARGV;
my ($line0,$read_name,$valid,$line1,$line2,$n_line,$picard,%fastq,$path,$normal_bam,$tumor_bam,$type);
my ($zcat,$bam2fastq,$pm,$ppm,$pid,$command,$normal_output,$tumor_output);

########################  prepare  ##########################

%fastq=("tumor"=>[$fastq1_tumor,$fastq2_tumor],"normal"=>[$fastq1_normal,$fastq2_normal]);
$path=abs_path($0);
$path=~s/cnv\.pl//;
$bam2fastq=$path."/somatic_script/bam2fastq.pl";
$picard=$path."/somatic_script/picard.jar";

$normal_output=$output."/normal";
$tumor_output=$output."/tumor";

system("date");
system_call("rm -f -r ".$output);
system_call("mkdir ".$output);
if (!-w $output) {die "Error: directory ".$output." is not writable!\n";}
system_call("mkdir ".$normal_output);
system_call("mkdir ".$tumor_output);

################################### analysis ###################################

# alignment
$pm=Parallel::ForkManager->new(2);
foreach $type (keys %fastq)
{
  $pid = $pm -> start and next;      #fork and skip to open a new fork
  alignment($type);
  $pm->finish;    
}
$pm->wait_all_children;

# path of bam file
$normal_bam=abs_path($normal_output."/normal.bam");
$tumor_bam=abs_path($tumor_output."/tumor.bam");
if ((! -e $normal_bam) || (! -e $tumor_bam))
{
  print "At least one bam file doesn't exist!\n";
  exit;
}

# CNVkit
system_call("cnvkit.py batch ".$tumor_bam." --normal ".$normal_bam.
  " --targets ".$index.".exon.bed --fasta ".$index." --annotate ".$index."_cnvkit.txt ".
  " --access ".$index."_access-excludes.bed --output-dir ".$output);
system("rm -f ".$output."/*.bed");
system("rm -f ".$output."/reference*");
system("rm -f ".$output."/*antitargetcoverage*");
system_call("Rscript ".$path."/somatic_script/clean_cnv.R ".$output."/tumor.cns ".
  $index."_cnvkit.txt ".$somatic." ".$output."/CNV_gene_level.txt");

# keep bam
if (1==1)
{
  system("rm -f ".$normal_bam);
  system("rm -f ".$tumor_bam);
  system("rm -f ".$normal_bam.".bai");
  system("rm -f ".$tumor_bam.".bai");
}

################  functions  ######################

sub alignment{
    my $type = $_[0];
    my $type_output = $output.'/'.$type;

    # convert bam files
    if (${$fastq{$type}}[0] eq "bam") 
    {
      system_call("perl ".$bam2fastq." ".${$fastq{$type}}[1]." ".$type_output." ".$thread);
    }else
    {
      system_call("zcat ".${$fastq{$type}}[0]." > ".$type_output."/fastq1.fastq");
      system_call("zcat ".${$fastq{$type}}[1]." > ".$type_output."/fastq2.fastq");
    }

    ${$fastq{$type}}[0]=$type_output."/fastq1.fastq";
    ${$fastq{$type}}[1]=$type_output."/fastq2.fastq";
  
    if (! -e ${$fastq{$type}}[0])
    { 
      print "Fastq file doesn't exist!\n";
      exit;
    }

    # qc
    system_call("mkdir ".$type_output."/fastqc");
    system_call("fastqc -o ".$type_output."/fastqc --extract -t ".$thread." -q -d ".$type_output."/fastqc ".
      ${$fastq{$type}}[0]." ".${$fastq{$type}}[1]);

    # number of lines in a file
    open(NL,"< ".${$fastq{$type}}[0]);
    $n_line=0;
    $n_line++ while <NL>;
    close(NL);

    # other filtering 
    open(FQ_IN1,${$fastq{$type}}[0]);
    open(FQ_OUT1,">".$type_output."/fastq1.fastq.tmp");
    open(FQ_IN2,${$fastq{$type}}[1]);
    open(FQ_OUT2,">".$type_output."/fastq2.fastq.tmp");

    $line1=$line2="";
    $valid=1;
      
    while ($line1=<FQ_IN1>)
    {
      $line1=~s/\s.*//; # keep only read name
      $line1=~s/\n//;
      $line1.="\n";
      $read_name=$line1;
      $line0=<FQ_IN1>;if (length($line0)<=5) {$valid=0;};$line1.=$line0;
      $line0=<FQ_IN1>;$line1.="+\n";
      $line1.=<FQ_IN1>;

      $line2=<FQ_IN2>;
      unless (defined $line2) {print "Error: fq2 shorter than fq1!\n";last;}
      $line2=$read_name; # force read name to match between R1 and R2
      $line0=<FQ_IN2>;if (length($line0)<=5) {$valid=0;};$line2.=$line0;
      $line0=<FQ_IN2>;$line2.="+\n";
      $line2.=<FQ_IN2>;

      if (rand($n_line/4/$read_ct_max)>1) {$valid=0;}

      if ($valid==1)
      {
        print FQ_OUT1 $line1;
        print FQ_OUT2 $line2;
      }
      $valid=1;
    }

    close(FQ_IN1);
    close(FQ_OUT1);
    close(FQ_IN2);
    close(FQ_OUT2);

    system_call("mv ".$type_output."/fastq1.fastq.tmp ".$type_output."/fastq1.fastq");
    system_call("mv ".$type_output."/fastq2.fastq.tmp ".$type_output."/fastq2.fastq");

    if ($pdx=~/PDX/)
    {
      system_call("source activate ".$disambiguate."/conda_env;".
        "python ".$disambiguate."/ngs_disambiguate.py -o ".$type_output.
        " -i ".$type_output."/disambiguate -a bwa -r \"".$index."|".$index."_mouse/mm10.fasta\"".
        " -n ".$thread." -b ".$path."/somatic_script/bam2fastq.pl ".
        ${$fastq{$type}}[0]." ".${$fastq{$type}}[1].";source deactivate");
      system_call("rm -f -r ".$type_output."/alignment_human");
      system_call("rm -f -r ".$type_output."/alignment_mouse");
      system_call("rm -f -r ".$type_output."/fastq1.fastq.mouse");
      system_call("rm -f -r ".$type_output."/fastq2.fastq.mouse");
      system_call("mv ".$type_output."/fastq1.fastq.human ".${$fastq{$type}}[0]);
      system_call("mv ".$type_output."/fastq2.fastq.human ".${$fastq{$type}}[1]);
      system_call("rm -f -r ".$type_output."/disambiguate");
    }

    # align
    system_call("bwa mem -v 1 -t ".$thread." -a -M ".$index." ".${$fastq{$type}}[0]." ".${$fastq{$type}}[1]." > ".$type_output."/alignment.sam");
    system_call("rm -f ".$type_output."/fastq1.fastq*");
    system_call("rm -f ".$type_output."/fastq2.fastq*");

    # remove decoy
    open(FILE_SAM,$type_output."/alignment.sam");
    open(FILE_SAM1,">".$type_output."/alignment.sam1");

    while (<FILE_SAM>)
    {
      if ($_=~/^@/ || $_=~/chr[0-9XY]/) {print FILE_SAM1 $_;}
    }

    close(FILE_SAM);
    close(FILE_SAM1);

    # sam to bam
    system_call("mv ".$type_output."/alignment.sam1 ".$type_output."/alignment.sam");
    system_call("sambamba view -f bam -h -v -S -l 0 -t ".$thread." ".$type_output."/alignment.sam -o ".$type_output."/alignment.bam");
    unlink_file($type_output."/alignment.sam");

    # sort
    system_call("sambamba sort --memory-limit=64GB --tmpdir=".$type_output."/sambamba_tmp -o ".$type_output."/sort.bam -t ".$thread.
     " -l 0 -u ".$type_output."/alignment.bam");
    unlink_file($type_output."/alignment.bam");

    # markdup
    #system_call("sambamba markdup -l 9 -r --overflow-list-size=400000 --io-buffer-size=1024 -t ".$thread.
    #  " --tmpdir ".$type_output."/tmp ".$type_output."/sort.bam ".$type_output."/".$type.".bam");
    system_call("java -jar ".$picard." MarkDuplicates INPUT=".$type_output."/sort.bam OUTPUT=".$type_output."/".$type.".bam ".
        " CREATE_INDEX=true VALIDATION_STRINGENCY=STRICT REMOVE_SEQUENCING_DUPLICATES=true METRICS_FILE=".$type_output."/dupmark_metrics.txt");
    system_call("rm -f -r ".$type_output."/dupmark_metrics.txt");
    unlink_file($type_output."/sort.bam");

    # generate final bam and mpileup files
    system("rm -f -r ".$type_output."/tmp");
    system("rm -f -r ".$type_output."/sambamba_tmp");
    system("rm -f -r ".$type_output."/*bai");
}

sub unlink_file
{
  my $file=$_[0];
  unless (-f $file) {die $file." doesn't exist!\n";}
  unlink($file);
}

sub system_call
{
  my $command=$_[0];
  print "\n\n".$command."\n\n";
  system($command);
}

#perl /home2/twang6/software/cancer/somatic/cnv.pl \
#/project/BICF/shared/Kidney/Projects/ActiveSurveillance/Exome/RAW/9A_GATGAATC_ACBLLMANXX_L008_001.R1.fastq.gz \
#/project/BICF/shared/Kidney/Projects/ActiveSurveillance/Exome/RAW/9A_GATGAATC_ACBLLMANXX_L008_001.R2.fastq.gz \
#/project/BICF/shared/Kidney/Projects/ActiveSurveillance/Exome/RAW/9B_CTAAGGTC_ACBLLMANXX_L005_001.R1.fastq.gz \
#/project/BICF/shared/Kidney/Projects/ActiveSurveillance/Exome/RAW/9B_CTAAGGTC_ACBLLMANXX_L005_001.R2.fastq.gz \
#32 /project/shared/xiao_wang/data/hg38/hs38d1.fa \
#/project/bioinformatics/Xiao_lab/shared/neoantigen/data/somatic/exome_seq/005T/somatic_mutations_hg38.txt \
#/project/bioinformatics/Xiao_lab/shared/neoantigen/data/tmp2
#human /project/shared/xiao_wang/software/disambiguate_pipeline
