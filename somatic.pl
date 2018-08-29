#####################  human somatic mutation calling ########################
#
# The instructions must be followed exactly!!! 
# prerequisite in path: 
# Rscript, bwa (>=0.7.15), STAR (if using RNA-Seq data), sambamba, speedseq, varscan, samtools (>=1.6), shimmer
# annovar (database downloaded in default folder: refGene,ljb26_all,cosmic70,esp6500siv2_all,exac03,1000g2015aug)
# python (2), strelka (>=2.8.3, note: strelka is tuned to run exome-seq or RNA-seq), manta(>=1.4.0), java (1.8)
# perl (need Parallel::ForkManager), lofreq_star (>=2.1.3, for tumor-only calling)
# input format: 
# fastq files: (1) fastq1 and fastq2 of normal sample, fastq1 and fastq2 of tumor sample (must be .gz)
#                  default input is full path to the 4 fastq files for tumor and normal samples 
#              (2) if need to directly input bam files, use "bam path_to_bam.bam" in replace of the two corresponding fastq input files
#                  can be a mixture of fastq and bam input
#              (3) if RNA-Seq data are used, use "RNA:fastq1" or "RNA:bam" at the first or third slot
#              (4) if Agilent SureSelect (Deep exome sequencing) data, use "Deep:fastq1" at the first or third slot
#                  optional: run somatic_script/SurecallTrimmer.jar on the fastq files before running somatic.pl
#              (5) for tumor-only calling, put "NA NA" in the slots of the normal samples. Results will be written to *germline* files
# thread: number of threads to use. Recommended: 32
# build: human genome build, hg19 or hg38
# index: path (including file name) to the human reference genome
# java17: path (including the executable file name) to java 1.7 (needed only for mutect)
# output: the output folder, it will be deleted (if pre-existing) and re-created during analysis
# need at least 128GB of memory
# commonly appearing human genetics polymorphism will be filtered out
#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';
use File::Copy;
use Parallel::ForkManager;

my ($fastq1_normal,$fastq2_normal,$fastq1_tumor,$fastq2_tumor,$thread,$build,$index,$java17,$output)=@ARGV;
my ($line0,$line1,$line2,$read_name,$valid,$path,$picard,$mutect,$gatk,$resource,$dict,$locatit,$annovar_path,$rna,$star_index);
my ($resource_1000g,$resource_mills,$resource_dbsnp,$resource_cosmic,$normal_bam,$tumor_bam,$type);
my ($strelka_exome,$zcat,$bam2fastq,$pm,$ppm,$pid,$command);
my ($normal_output,$tumor_output,$mutect_tmp,$speed_tmp);

my $read_ct_max=75000000; # put a cap on how many reads to use, for sake of memory
my $read_ct=0;

my $manta="";
my @command_array=("fun_mutect();","fun_speed_var_shimmer();","fun_strelka();");
my %vcfs=("passed.somatic.indels.vcf"=>"strelka","passed.somatic.snvs.vcf"=>"strelka","mutect.vcf"=>"mutect",
  "somatic_diffs.readct.vcf"=>"shimmer","speedseq2.vcf"=>"speedseq",
  "varscan.indel.Somatic.hc.vcf"=>"varscan","varscan.snp.Somatic.hc.vcf"=>"varscan");
my %fastq=("tumor"=>[$fastq1_tumor,$fastq2_tumor],"normal"=>[$fastq1_normal,$fastq2_normal]);

##############################  path of files  ###############################

$resource=$index."_resource";
$dict=$index;
$dict=~s/fa$/dict/;
$star_index=$index;
$star_index=~s/\/[^\/]*?$/\/STAR/;

$path=abs_path($0);
$path=~s/somatic\.pl//;
$gatk=$path."/somatic_script/GenomeAnalysisTK.jar";
$mutect=$path."/somatic_script/mutect-1.1.7.jar";
$picard=$path."/somatic_script/picard.jar";
$bam2fastq=$path."/somatic_script/bam2fastq.pl";
$locatit=$path."/somatic_script/LocatIt_v4.0.1.jar";

$resource_1000g=$resource."/1000G_phase1.snps.high_confidence.".$build.".vcf";
$resource_mills=$resource."/Mills_and_1000G_gold_standard.indels.".$build.".vcf";
$resource_dbsnp=$resource."/dbsnp.".$build.".vcf";
$resource_cosmic=$resource."/CosmicCodingMuts.".$build.".vcf";

$strelka_exome="--exome";

########################  prepare the work folders  ##########################

$normal_output=$output."/normal";
$tumor_output=$output."/tumor";
$mutect_tmp=$output."/mutmp";
$speed_tmp=$output."/sptmp";

system("date");
system_call("rm -f -r ".$output);
system_call("mkdir ".$output);
system_call("mkdir ".$output."/intermediate");
if (!-w $output) {die "Error: directory ".$output." is not writable!\n";}
system_call("mkdir ".$normal_output);
system_call("mkdir ".$tumor_output);
system_call("mkdir ".$mutect_tmp);
system_call("mkdir ".$speed_tmp);

################################### alignment ###################################

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
if (((! -e $normal_bam) && ${$fastq{"normal"}}[0] ne "NA") || (! -e $tumor_bam))
{
  print "Error: At least one bam file doesn't exist!\n";
  exit;
}

################################# somatic calling ################################

if (${$fastq{"normal"}}[0] ne "NA")
{
  $ppm=Parallel::ForkManager->new(3);
  foreach $command (@command_array)
  {
    $pid = $ppm -> start and next;
    eval $command;
    $ppm->finish;  
  }
  $ppm->wait_all_children;
}else #tumor only calling
{
  system_call("lofreq call-parallel -s -f ".$index." -S ".$resource_dbsnp.".gz --call-indels ".
    " -l ".$index.".exon.bed --pp-threads ".$thread." -o ".$output."/lofreq.vcf ".$tumor_bam);
}

###############################  integrate vcf files  #############################

# left alignment
if (${$fastq{"normal"}}[0] ne "NA")
{ 
  $ppm=Parallel::ForkManager->new(7);
  foreach (keys %vcfs)
  {
    $pid = $ppm -> start and next;
    system_call("java -jar ".$gatk." -T LeftAlignAndTrimVariants -R ".$index." --variant ".$output."/".$_." -o ".
      $output."/left_".$_);
    system_call("cp ".$output."/left_".$_." ".$output."/intermediate/left_".$_);
    $ppm->finish;
  }
  $ppm->wait_all_children;
}

# filter vcfs
system_call("Rscript ".$path."/somatic_script/filter_vcf.R ".$output." ".${$fastq{"normal"}}[0]);

# annovar
open(AV,"which table_annovar.pl |") or die "Cannot call table_annovar.pl!\n";
$annovar_path=<AV>;
close(AV);
$annovar_path=~s/table_annovar.pl\n//;

if (${$fastq{"normal"}}[0] ne "NA")
{ 
  system_call("table_annovar.pl ".$output."/somatic_mutations.txt ".$annovar_path."/humandb/ -buildver ".$build." -out ".$output.
    "/somatic_mutations -remove -protocol refGene,ljb26_all,cosmic70,esp6500siv2_all,exac03,1000g2015aug_all".
    " -operation g,f,f,f,f,f -nastring .");
  system_call("Rscript ".$path."/somatic_script/filter_vcf2.R ".$output."/somatic_mutations.txt ".$output."/somatic_mutations.".
    $build."_multianno.txt ".$output."/somatic_mutations_".$build.".txt");
  unlink($output."/somatic_mutations.txt");
  unlink($output."/somatic_mutations.".$build."_multianno.txt");
}

system_call("table_annovar.pl ".$output."/germline_mutations.txt ".$annovar_path."/humandb/ -buildver ".$build." -out ".$output.
  "/germline_mutations -remove -protocol refGene,ljb26_all,cosmic70,esp6500siv2_all,exac03,1000g2015aug_all".
  " -operation g,f,f,f,f,f -nastring .");
system_call("Rscript ".$path."/somatic_script/filter_vcf2.R ".$output."/germline_mutations.txt ".$output."/germline_mutations.".
  $build."_multianno.txt ".$output."/germline_mutations_".$build.".txt");
unlink($output."/germline_mutations.txt");
unlink($output."/germline_mutations.".$build."_multianno.txt");

# clean up
system("date");
system("rm -f ".$output."/*out");
system("rm -f ".$output."/*idx");
system("rm -f ".$output."/varscan.*vcf");
system("rm -f ".$output."/mutect.vcf");
system("rm -f ".$output."/speedseq*");
system("rm -f ".$output."/som_counts*");
system("rm -f ".$output."/passed_*");
system("rm -f ".$output."/indel_counts*");
system("rm -f ".$output."/het_counts*");
system("rm -f ".$output."/somatic_diffs*");
system("rm -f ".$output."/somatic_indels.vs");
system("rm -f ".$output."/passed.somatic.*");
system("rm -rf ".$normal_output);
system("rm -rf ".$tumor_output);
system("rm -f ".$output."/indel_counts*");
system("rm -f ".$output."/het_counts*");
system("rm -f ".$output."/somatic_diffs*");
system("rm -f ".$output."/somatic_indels.vs");
system("rm -f ".$output."/passed.somatic.*");
system("rm -rf ".$normal_output);
system("rm -rf ".$tumor_output);
system("rm -rf ".$mutect_tmp);
system("rm -rf ".$speed_tmp);

################  functions  ######################

sub alignment{
    my $type = $_[0];
    if (${$fastq{$type}}[0] eq "NA") {return(0);} # tumor-only calling
    my $type_output = $output.'/'.$type;

    # check whether tumor sample is RNA-Seq
    $rna=0; # exome-seq
    if (${$fastq{$type}}[0]=~s/RNA\://) {$rna=1;} # rna-seq
    if (${$fastq{$type}}[0]=~s/Deep\://) {$rna=2;} # deep exome-seq

    if (${$fastq{$type}}[0] eq "bam")  # convert bam files 
    {
      system_call("perl ".$bam2fastq." ".${$fastq{$type}}[1]." ".$type_output." ".$thread);
    }else # or pre-process fastq.gz files
    {
      open(FQ_IN1,"zcat ".${$fastq{$type}}[0]." |");
      open(FQ_OUT1,">".$type_output."/fastq1.fastq");
      open(FQ_IN2,"zcat ".${$fastq{$type}}[1]." |");
      open(FQ_OUT2,">".$type_output."/fastq2.fastq");

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
        $line2=$read_name; # force read name to match between R1 and R2
        $line0=<FQ_IN2>;if (length($line0)<=5) {$valid=0;};$line2.=$line0;
        $line0=<FQ_IN2>;$line2.="+\n";
        $line2.=<FQ_IN2>;

        if ($read_ct++>$read_ct_max) {last;}

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
    }  
  
    ${$fastq{$type}}[0]=$type_output."/fastq1.fastq";
    ${$fastq{$type}}[1]=$type_output."/fastq2.fastq";
    $zcat="";
  
    # check again
    if (! -e ${$fastq{$type}}[0])
    { 
      print "Error: Fastq file doesn't exist!\n";
      exit;
    }

    # align
    if ($rna==0 || $rna==2)
    {
        system_call("bwa mem -v 1 -t ".$thread." -a -M ".$index." ".${$fastq{$type}}[0]." ".${$fastq{$type}}[1]." > ".$type_output."/alignment.sam");
    }else
    {
        system_call("mkdir ".$type_output."/tmp");
        system_call("mkdir ".$type_output."/tmp/pass2");
        system_call("STAR --genomeDir ".$star_index." --readFilesIn ".${$fastq{$type}}[0]." ".${$fastq{$type}}[1]." --runThreadN ".$thread.
        " ".$zcat." --outFileNamePrefix ".$type_output."/tmp/pass1");
        system_call("cd ".$type_output.";STAR --runMode genomeGenerate --genomeDir ".$type_output."/tmp/pass2 --genomeFastaFiles ".$index." --sjdbFileChrStartEnd ".
        $type_output."/tmp/pass1SJ.out.tab --sjdbOverhang 75 --runThreadN ".$thread);
        system_call("STAR --genomeDir ".$type_output."/tmp/pass2 --readFilesIn ".${$fastq{$type}}[0]." ".${$fastq{$type}}[1]." --runThreadN ".$thread.
        " ".$zcat." --outFileNamePrefix ".$type_output."/");
        system_call("mv ".$type_output."/Aligned.out.sam ".$type_output."/alignment.sam");
        unlink_file($type_output."/SJ.out.tab");
    }

    #add read group
    system_call("java -Djava.io.tmpdir=".$type_output."/tmp -jar ".$picard." AddOrReplaceReadGroups INPUT=".$type_output."/alignment.sam OUTPUT=".
        $type_output."/rgAdded.bam SORT_ORDER=coordinate RGID=".$type." RGLB=".$type." RGPL=illumina RGPU=".$type." RGSM=".$type.
        " CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT COMPRESSION_LEVEL=0");
    unlink_file($type_output."/alignment.sam");

    # mark duplicates
    if ($rna==0 || $rna==1)
    {
      system_call("sambamba markdup -l 0 -r --overflow-list-size=400000 --io-buffer-size=256 -t ".$thread.
        " --tmpdir ".$type_output."/tmp ".$type_output."/rgAdded.bam ".$type_output."/dupmark.bam");
    }else # SureSelect deep sequencing
    {
      system_call("samtools view -H ".$type_output."/rgAdded.bam > ".$type_output."/header.sam");
      system_call("java -Xmx200g -jar ".$locatit." -X ".$type_output."/tmp -PM:xm,Q:xq,q:nQ,r:nR,I:ni ".
        "-U -q 25 -m 1 -IB -OB -C -i -r -c 2500 -H ".$type_output."/header.sam ".
        "-o ".$type_output."/locatit.bam ".$type_output."/rgAdded.bam ".${$fastq{$type}}[0]." ".${$fastq{$type}}[1]);
      system_call("sambamba sort --memory-limit=64GB --tmpdir=".$type_output."/tmp -o ".$type_output."/dupmark.bam -t ".$thread.
        " -l 0 -u ".$type_output."/locatit.bam");
      system_call("sambamba index -t ".$thread." ".$type_output."/dupmark.bam");     
      unlink_file($type_output."/header.sam");
      unlink_file($type_output."/locatit.bam");
    }
   
    system_call("rm -f ".$type_output."/fastq1.fastq*");
    system_call("rm -f ".$type_output."/fastq2.fastq*");
    unlink_file($type_output."/rgAdded.bam");
    unlink_file($type_output."/rgAdded.bai");

    if ($rna==1)
    {
        system_call("java -jar ".$gatk." -T SplitNCigarReads -R ".$index." -I ".$type_output."/dupmark.bam -o ".$type_output."/split.bam -rf ReassignOneMappingQuality ".
        "-RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS");
        system_call("mv ".$type_output."/split.bam ".$type_output."/dupmark.bam");
        system_call("mv ".$type_output."/split.bai ".$type_output."/dupmark.bam.bai");
    }

    # indel realignment
    system_call("java -Djava.io.tmpdir=".$type_output."/tmp -jar ".$gatk." -T RealignerTargetCreator -R ".$index." --num_threads ".$thread." -known ".$resource_mills.
        " -known ".$resource_1000g." --bam_compression 0 -o ".$type_output."/".$type."_intervals.list -I ".$type_output."/dupmark.bam > ".$type_output."/index.out");

    system_call("java -Djava.io.tmpdir=".$type_output."/tmp -jar ".$gatk." -T IndelRealigner --bam_compression 0 ".
        " --filter_bases_not_stored --disable_auto_index_creation_and_locking_when_reading_rods -R ".$index." -known ".
        $resource_mills." -known ".$resource_1000g." -targetIntervals ".$type_output."/".$type."_intervals.list -I ".
        $type_output."/dupmark.bam -o ".$type_output."/realigned.bam >".$type_output."/".$type."_realign.out");

    unlink_file($type_output."/dupmark.bam");
    unlink_file($type_output."/dupmark.bam.bai");

    # base recalibration
    system_call("java -Djava.io.tmpdir=".$type_output."/tmp -jar ".$gatk." -T BaseRecalibrator -R ".$index." --bam_compression 0 ".
        " -knownSites ".$resource_dbsnp." -knownSites ".$resource_mills." -I ".$type_output."/realigned.bam -o ".$type_output."/".
        $type."_bqsr > ".$type_output."/table.out");
    system_call("java -Djava.io.tmpdir=".$type_output."/tmp -jar ".$gatk." -T PrintReads -rf NotPrimaryAlignment --bam_compression 0 -R ".
        $index." -I ".$type_output."/realigned.bam -BQSR ".$type_output."/".$type."_bqsr -o ".
        $type_output."/".$type.".bam > ".$type_output."/".$type."_recal.out");

    unlink_file($type_output."/realigned.bam");
    unlink_file($type_output."/realigned.bai");
    unlink_file($type_output."/".$type."_bqsr");
    unlink_file($type_output."/".$type."_intervals.list");

    # further process bam files
    system_call("sambamba index -t ".$thread." ".$type_output."/".$type.".bam");    		
    system_call("sambamba mpileup --tmpdir ".$type_output."/tmp -t ".$thread." --buffer-size=2048000000 ".		
    $type_output."/".$type.".bam --samtools \"-C 50 -f ".$index."\" > ".$type_output."/".$type.".mpileup");		
    if (-s $type_output."/".$type.".mpileup"<100000) {die "Error: mpileup has failed!\n";}		
}		

sub fun_mutect{
  system_call($java17." -Djava.io.tmpdir=".$mutect_tmp." -Xmx32g -jar ".$mutect." --analysis_type MuTect --reference_sequence ".$index.
    " --dbsnp ".$resource_dbsnp." --cosmic ".$resource_cosmic." --input_file:tumor ".$tumor_bam." --input_file:normal ".$normal_bam.
    " --vcf ".$output."/mutect.vcf --out ".$output."/mutect.out");
}

sub fun_speed_var_shimmer{
  # speedseq
  system_call("speedseq somatic -F 0.01 -q 10 -t ".$thread." -T ".$speed_tmp." -o ".$output."/speedseq ".$index." ".$normal_bam." ".$tumor_bam);
  system_call("gzip -d ".$output."/speedseq.vcf.gz");
  system_call("java -jar ".$picard." UpdateVcfSequenceDictionary I=".$output."/speedseq.vcf O=".$output."/speedseq1.vcf SEQUENCE_DICTIONARY=".$dict);
  system_call("java -jar ".$picard." SortVcf I=".$output."/speedseq1.vcf O=".$output."/speedseq2.vcf SEQUENCE_DICTIONARY=".$dict);

  # varscan
  system_call("VarScan somatic ".$normal_output."/normal.mpileup ".$tumor_output."/tumor.mpileup ".$output."/varscan --output-vcf 1");
  system_call("VarScan processSomatic ".$output."/varscan.indel.vcf --min-tumor-freq 0.01");
  system_call("VarScan processSomatic ".$output."/varscan.snp.vcf --min-tumor-freq 0.01");

  # shimmer 
  system_call("shimmer.pl --minqual 25 --ref ".$index." ".$normal_bam." ".$tumor_bam." --outdir ".$output);
  system_call("perl ".$path."/somatic_script/add_readct_shimmer.pl ".$output."/som_counts.bh.txt ".$output."/somatic_diffs.vcf ".$output.
    "/somatic_diffs.readct.vcf");
}

sub fun_strelka{
  system_call("rm -f -r ".$output."/manta");
  system_call("rm -f -r ".$output."/strelka");

  #run manta first
  system_call("configManta.py --normalBam ".$normal_bam." --tumorBam ".$tumor_bam.
    " --referenceFasta ".$index." --runDir ".$output."/manta ".$strelka_exome);
  system_call($output."/manta/runWorkflow.py -m local -j ".$thread);
  if (-e $output."/manta/results/variants/candidateSmallIndels.vcf.gz") # manta tends to fail for mate-pair sequencing 
    {$manta=" --indelCandidates ".$output."/manta/results/variants/candidateSmallIndels.vcf.gz";}

  # run strelka2
  system_call("configureStrelkaSomaticWorkflow.py --normalBam ".$normal_bam." --tumorBam ".$tumor_bam." --referenceFasta ".$index.
    " --runDir ".$output."/strelka ".$strelka_exome." ".$manta);
  system_call($output."/strelka/runWorkflow.py -m local -j ".$thread);
  system_call("gzip -d ".$output."/strelka/results/variants/*gz");
  system_call("mv ".$output."/strelka/results/variants/somatic.indels.vcf ".$output."/passed.somatic.indels.vcf"); # legacy codes
  system_call("mv ".$output."/strelka/results/variants/somatic.snvs.vcf ".$output."/passed.somatic.snvs.vcf"); # legacy codes
  
  system_call("rm -f -r ".$output."/strelka");
  system_call("rm -f -r ".$output."/manta"); 
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

##############  example code  #################

#perl /home2/twang6/software/cancer/somatic/somatic.pl \
#/archive/BICF/shared/Kidney/exome/RAW/SAM9259823_1_1.gne.fastq.gz \
#/archive/BICF/shared/Kidney/exome/RAW/SAM9259823_1_2.gne.fastq.gz \
#/archive/BICF/shared/Kidney/exome/RAW/SAM19944142_1_1.gne.fastq.gz \
#/archive/BICF/shared/Kidney/exome/RAW/SAM19944142_1_2.gne.fastq.gz \
#32 hg38 /home2/twang6/data/genomes/hg38/hs38d1.fa \
#/cm/shared/apps/java/oracle/jdk1.7.0_51/bin/java \
#~/iproject/test/
