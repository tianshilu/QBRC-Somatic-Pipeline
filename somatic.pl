#####################  human somatic mutation calling ########################
#
# The instructions must be followed exactly!!! 
# Need at least 32GB of memory for DNA jobs and 128GB for RNA jobs
# Commonly appearing SNVs in the population are filtered out from the somatic output
# if inputing bam files, the bam files are assumed to be paired-ended
#
# prerequisite in path: 
# Rscript, bwa (>=0.7.15), STAR (>=2.7.2), sambamba, speedseq, varscan, samtools (>=1.6), shimmer
# annovar (>=2019Oct24,refGene,ljb26_all,cosmic70,esp6500siv2_all,exac03,1000g2015aug downloaded in humandb and refGene downloaded in mousedb)
# python (2.7), strelka (>=2.8.3, note: strelka is tuned to run exome-seq or RNA-seq), manta(>=1.6.0), java (1.8)
# perl (need Parallel::ForkManager), lofreq_star (>=2.1.3), bowtie2 (>=2.3.4.3, for PDX mode)
# disambiguate (use conda env or install by the corresponding env.yml, conda env create -f environment.yml)
#
# input format: 
# fastq files: (1) fastq1 and fastq2 of normal sample, fastq1 and fastq2 of tumor sample (must be .gz)
#                  default input is full path to the 4 fastq files for tumor and normal samples 
#              (2) if need to directly input bam files, use "bam path_to_bam.bam" in replace of the two corresponding fastq input files
#                  can be a mixture of fastq and bam input
#              (3) if RNA-Seq data are used, use "RNA:fastq" or "RNA:bam" at the first or third slot
#              (4) if Agilent SureSelect (Deep exome sequencing) data, use "Deep:fastq1" at the first or third slot
#                  optional: run somatic_script/SurecallTrimmer.jar on the fastq files before running somatic.pl
#              (5) for tumor-only calling, put "NA NA" in the slots of the normal samples. Results will be written to *germline* files
#                  mutation callling is super sensitive and not specific in this case
#              (6) If only single end fastq data are available, put the fastq file(s) at the first and/or the third slots, 
#                  then put NA in the second and/or fourth slot
# thread: number of threads to use. Recommended: 32
# build: hg19 or hg38 or mm10
# index: path (including file name) to the human/mouse reference genome
# java17: path (including the executable file name) to java 1.7 (needed only for mutect). must be strictly followed
# output: the output folder, it will be deleted (if pre-existing) and re-created during analysis
# pdx: "PDX" or "human" or "mouse". if "PDX", the tumor (not normal) sample will be processed by disambiguate. But it can only handle paired-end sequencing reads
# keep_coverage: whether to keep per-base coverage information. Default is 0. Set to 1 to enable
# disambiguate: disambiguate path, /project/shared/xiao_wang/software/disambiguate_pipeline. or install your own version by env.yml
#
#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';
use File::Copy;
use Parallel::ForkManager;

# hidden paramters
my $lofreq=1; # 1 or 0. For tumor-only mode, use lofreq or strelka. For targeted panel sequencing, lofreq=1 is strongly suggested
my $skip_recal=0; # 1 or 0. Skip base recalibration or not
my $read_ct_max=50000000; # put a cap on how many reads to use, for sake of memory

my ($fastq1_normal,$fastq2_normal,$fastq1_tumor,$fastq2_tumor,$thread,$build,$index,
  $java17,$output,$pdx,$keep_coverage,$disambiguate)=@ARGV;
my ($line0,$line1,$line2,$read_name,$valid,$path,$picard,$mutect,$gatk,$resource,$dict,$locatit,$annovar_path,$rna,$star_index);
my ($resource_1000g,$resource_mills,$resource_dbsnp,$resource_cosmic,$normal_bam,$tumor_bam,$type);
my ($known,$strelka_exome,$zcat,$bam2fastq,$pm,$ppm,$pid,$command,$annovar_db,$annovar_protocol);
my ($normal_output,$tumor_output,$mutect_tmp,$speed_tmp);

my $manta="";
my @command_array=("fun_mutect();","fun_speed_var_shimmer();","fun_strelka_lofreq();");
my %vcfs=("passed.somatic.indels.vcf"=>"strelka","passed.somatic.snvs.vcf"=>"strelka","mutect.vcf"=>"mutect",
  "somatic_diffs.readct.vcf"=>"shimmer","speedseq2.vcf"=>"speedseq",
  "varscan.indel.Somatic.hc.vcf"=>"varscan","varscan.snp.Somatic.hc.vcf"=>"varscan",
  "lofreq_n.vcf"=>"lofreq","lofreq_t.vcf"=>"lofreq");
my %fastq=("tumor"=>[$fastq1_tumor,$fastq2_tumor],"normal"=>[$fastq1_normal,$fastq2_normal]);

##############################  path of files  ###############################

$resource=$index."_resource";
$dict=$index;
$dict=~s/fa$/dict/;$dict=~s/fasta$/dict/;
$star_index=$index;
$star_index=~s/\/[^\/]*?$/\/STAR/;

$path=abs_path($0);
$path=~s/somatic\.pl//;
$gatk=$path."/somatic_script/GenomeAnalysisTK.jar";
$mutect=$path."/somatic_script/mutect-1.1.7.jar";
$picard=$path."/somatic_script/picard.jar";
$bam2fastq=$path."/somatic_script/bam2fastq.pl";
$locatit=$path."/somatic_script/LocatIt_v4.0.1.jar";

$resource_1000g=$resource."/1000G_phase1.snps.high_confidence.".$build.".vcf"; # human only
$resource_mills=$resource."/Mills_and_1000G_gold_standard.indels.".$build.".vcf"; # human only
$resource_dbsnp=$resource."/dbsnp.".$build.".vcf"; # human+mouse
$resource_cosmic=$resource."/CosmicCodingMuts.".$build.".vcf"; # human only

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
  if ($lofreq==1) # two options
  {
    system_call("lofreq call -s --sig 1 --bonf 1 -C 5 -f ".$index." --call-indels ".
      " --verbose --use-orphan -o ".$output."/lofreq.vcf ".$tumor_bam);
  }else
  {
    system_call("configureStrelkaGermlineWorkflow.py --bam ".$tumor_bam." --referenceFasta ".$index.
      " --runDir ".$output."/strelka/ ".$strelka_exome); 
    system_call($output."/strelka/runWorkflow.py -m local -j ".$thread);
    system_call("cp ".$output."/strelka/results/variants/variants.vcf.gz ".$output);
    system_call("cp ".$output."/strelka/results/variants/variants.vcf.gz ".$output."/intermediate");
    system_call("gzip -d ".$output."/variants.vcf.gz");
    system_call("rm -f -r ".$output."/strelka");
  }
}

###############################  integrate vcf files  #############################

# left alignment
if (${$fastq{"normal"}}[0] ne "NA")
{ 
  $ppm=Parallel::ForkManager->new(9);
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

system_call("cp ".$output."/varscan.indel.Germline.vcf ".$output."/intermediate/");
system_call("cp ".$output."/varscan.snp.Germline.vcf ".$output."/intermediate/");
system_call("cp ".$output."/varscan.indel.LOH.vcf ".$output."/intermediate/");
system_call("cp ".$output."/varscan.snp.LOH.vcf ".$output."/intermediate/");

# filter vcfs
open(AV,"which table_annovar.pl |") or die "Cannot call table_annovar.pl!\n";
$annovar_path=<AV>;
close(AV);
$annovar_path=~s/table_annovar.pl\n//;

if ($pdx=~/mouse/)
{
  $annovar_db="/mousedb/";
  $annovar_protocol=" -protocol refGene -operation g ";
}else
{
  $annovar_db="/humandb/";
  $annovar_protocol=" -protocol refGene,ljb26_all,cosmic70,esp6500siv2_all,exac03,1000g2015aug_all -operation g,f,f,f,f,f ";
}

system_call("Rscript ".$path."/somatic_script/filter_vcf.R ".$output." ".${$fastq{"normal"}}[0].
  " ".$annovar_path." ".$annovar_db." ".$build);

# annovar
foreach $type (("germline","somatic"))
{
  system_call("table_annovar.pl ".$output."/".$type."_mutations.txt ".$annovar_path.$annovar_db." -buildver ".$build." -out ".$output.
    "/".$type."_mutations -remove ".$annovar_protocol." -nastring .");
  system_call("Rscript ".$path."/somatic_script/filter_vcf2.R ".$output."/".$type."_mutations.txt ".$output."/".$type."_mutations.".
    $build."_multianno.txt ".$output."/".$type."_mutations_".$build.".txt ".$type);
  unlink($output."/".$type."_mutations.txt");
  unlink($output."/".$type."_mutations.".$build."_multianno.txt");
}

# sv
system("rm -f ".$output."/*vcf");
#if (-e $output."/lumpy_sv.vcf.txt")
#{
#  system_call("grep -v IMPRECISE ".$output."/lumpy_sv.vcf.txt > ".$output."/lumpy_sv_".$build.".vcf");
#  system_call("rm -f ".$output."/lumpy_sv.vcf.txt");
#}

# clean up
system("date");
system("rm -f ".$output."/*out");
system("rm -f ".$output."/*idx");
system("rm -f ".$output."/speedseq*");
system("rm -f ".$output."/som_counts*");
system("rm -f ".$output."/passed_*");
system("rm -f ".$output."/indel_counts*");
system("rm -f ".$output."/het_counts*");
system("rm -f ".$output."/somatic_diffs*");
system("rm -f ".$output."/somatic_indels.vs");
system("rm -f ".$output."/passed.somatic.*");
if ($keep_coverage==1)
{
  system_call("cat ".$tumor_output."/tumor.mpileup | cut -f 1-4 > ".$output."/coverage.txt");
}
system("rm -rf ".$normal_output);
system("rm -rf ".$tumor_output);
system("rm -f ".$output."/indel_counts*");
system("rm -f ".$output."/het_counts*");
system("rm -f ".$output."/somatic_diffs*");
system("rm -f ".$output."/somatic_indels.vs");
system("rm -f ".$output."/passed.somatic.*");
system("rm -rf ".$mutect_tmp);
system("rm -rf ".$speed_tmp);

################  functions  ######################

sub alignment{
    my $type = $_[0];
    if (${$fastq{$type}}[0] eq "NA") {return(0);} # tumor-only calling
    my $type_output = $output.'/'.$type;
    my $n_line;

    # check whether tumor sample is RNA-Seq
    $rna=0; # exome-seq
    if (${$fastq{$type}}[0]=~s/RNA\://) {$rna=1;} # rna-seq
    if (${$fastq{$type}}[0]=~s/Deep\://) {$rna=2;} # deep exome-seq

    if (${$fastq{$type}}[0] eq "bam")  # convert bam files 
    {
      system_call("perl ".$bam2fastq." ".${$fastq{$type}}[1]." ".$type_output." ".$thread);
    }else # or pre-process fastq.gz files
    {
      # create R2 if single end sequencing data are provided
      if (${$fastq{$type}}[1] eq "NA")
      {
        system_call("perl ".$path."somatic_script/reverse_complement.pl ".${$fastq{$type}}[0]." ".
          $type_output."/R2.fastq.gz");
        ${$fastq{$type}}[1]=$type_output."/R2.fastq.gz";
      }

      # number of lines in a file
      open(NL,"zcat ".${$fastq{$type}}[0]." | wc -l |");
      $n_line=<NL>;
      $n_line=~s/\s.*//;
      close(NL);

      # other filtering
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

    # PDX
    if ($pdx=~/PDX/ && $type eq "tumor")
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
    if ($rna==0 || $rna==2)
    {
        system_call("bwa mem -v 1 -t ".$thread." -a -M ".$index." ".${$fastq{$type}}[0]." ".${$fastq{$type}}[1]." > ".$type_output."/alignment.sam");
    }else
    {
        system_call("mkdir ".$type_output."/tmp");
        system_call("mkdir ".$type_output."/tmp/pass2");
        system_call("STAR --genomeDir ".$star_index." --readFilesIn ".${$fastq{$type}}[0]." ".${$fastq{$type}}[1]." --runThreadN ".$thread.
          " ".$zcat." --outFileNamePrefix ".$type_output."/tmp/pass1 --outFilterMatchNminOverLread 0.4 --outFilterScoreMinOverLread 0.4");
        system_call("cd ".$type_output.";STAR --runMode genomeGenerate --genomeDir ".$type_output."/tmp/pass2 --genomeFastaFiles ".$index." --sjdbFileChrStartEnd ".
          $type_output."/tmp/pass1SJ.out.tab --sjdbOverhang 100 --runThreadN ".$thread);
        system_call("STAR --genomeDir ".$type_output."/tmp/pass2 --readFilesIn ".${$fastq{$type}}[0]." ".${$fastq{$type}}[1]." --runThreadN ".$thread.
          " ".$zcat." --outFileNamePrefix ".$type_output."/ --outFilterMatchNminOverLread 0.4 --outFilterScoreMinOverLread 0.4");
        system_call("mv ".$type_output."/Aligned.out.sam ".$type_output."/alignment.sam");
        unlink_file($type_output."/SJ.out.tab");
    }

    if ($rna!=2)
    {
      unlink_file(${$fastq{$type}}[0]);
      unlink_file(${$fastq{$type}}[1]);
    }

    #add read group
    system_call("java -Djava.io.tmpdir=".$type_output."/tmp -jar ".$picard." AddOrReplaceReadGroups INPUT=".$type_output."/alignment.sam OUTPUT=".
        $type_output."/rgAdded.bam SORT_ORDER=coordinate RGID=".$type." RGLB=".$type." RGPL=illumina RGPU=".$type." RGSM=".$type.
        " CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT");
    unlink_file($type_output."/alignment.sam");

    # mark duplicates
    if ($rna==0 || $rna==1)
    {
      #system_call("sambamba markdup -l 0 -r --overflow-list-size=400000 --io-buffer-size=1024 -t ".$thread.
      #  " --tmpdir ".$type_output."/tmp ".$type_output."/rgAdded.bam ".$type_output."/dupmark.bam");
      system_call("java -jar ".$picard." MarkDuplicates INPUT=".$type_output."/rgAdded.bam OUTPUT=".$type_output."/dupmark.bam ".
        " CREATE_INDEX=true VALIDATION_STRINGENCY=STRICT REMOVE_SEQUENCING_DUPLICATES=true METRICS_FILE=".$type_output."/dupmark_metrics.txt");
      system_call("rm -f -r ".$type_output."/dupmark_metrics.txt");
      system_call("mv ".$type_output."/dupmark.bai ".$type_output."/dupmark.bam.bai");
    }else # SureSelect deep sequencing
    {
      system_call("samtools view -H ".$type_output."/rgAdded.bam > ".$type_output."/header.sam");
      system_call("java -Xmx15g -jar ".$locatit." -X ".$type_output."/tmp -PM:xm,Q:xq,q:nQ,r:nR,I:ni ".
        "-U -N 800000 -q 25 -m 1 -IB -OB -C -i -r -c 2500 -H ".$type_output."/header.sam ".
        "-o ".$type_output."/locatit.bam ".$type_output."/rgAdded.bam ".${$fastq{$type}}[0]." ".${$fastq{$type}}[1]);
      system_call("sambamba sort --memory-limit=15GB --tmpdir=".$type_output."/tmp -o ".$type_output."/dupmark.bam -t ".$thread.
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
    if ($pdx=~/mouse/) {$known=" -known ".$resource_dbsnp." ";} else {$known=" -known ".$resource_mills." -known ".$resource_1000g." ";}
    system_call("java -Djava.io.tmpdir=".$type_output."/tmp -jar ".$gatk." -T RealignerTargetCreator -R ".$index." --num_threads ".$thread.
      $known." -o ".$type_output."/".$type."_intervals.list -I ".$type_output."/dupmark.bam > ".$type_output."/index.out");

    system_call("java -Djava.io.tmpdir=".$type_output."/tmp -jar ".$gatk." -T IndelRealigner ".
      " --filter_bases_not_stored --disable_auto_index_creation_and_locking_when_reading_rods -R ".$index.$known.
      " -targetIntervals ".$type_output."/".$type."_intervals.list -I ".
      $type_output."/dupmark.bam -o ".$type_output."/realigned.bam >".$type_output."/".$type."_realign.out");

    unlink_file($type_output."/dupmark.bam");
    unlink_file($type_output."/dupmark.bam.bai");

    # base recalibration
    if ($skip_recal==0)
    {
      if ($pdx=~/mouse/) {$known=" ";} else {$known=" -knownSites ".$resource_mills." ";}
      system_call("java -Djava.io.tmpdir=".$type_output."/tmp -jar ".$gatk." -T BaseRecalibrator -R ".$index." ".
        " -knownSites ".$resource_dbsnp.$known." -I ".$type_output."/realigned.bam -o ".$type_output."/".
        $type."_bqsr > ".$type_output."/table.out");
      system_call("java -Djava.io.tmpdir=".$type_output."/tmp -jar ".$gatk." -T PrintReads -rf NotPrimaryAlignment -R ".
        $index." -I ".$type_output."/realigned.bam -BQSR ".$type_output."/".$type."_bqsr -o ".
        $type_output."/".$type.".bam > ".$type_output."/".$type."_recal.out");

      unlink_file($type_output."/realigned.bam");
      unlink_file($type_output."/realigned.bai");
      unlink_file($type_output."/".$type."_bqsr");
      unlink_file($type_output."/".$type."_intervals.list");
    }else
    {
      system("mv ".$type_output."/realigned.bam ".$type_output."/".$type.".bam");
    }

    # further process bam files
    system_call("sambamba index -t ".$thread." ".$type_output."/".$type.".bam");    		
    #system_call("sambamba mpileup --tmpdir ".$type_output."/tmp -t ".$thread." --buffer-size=16000000000 ".		
    #  $type_output."/".$type.".bam --samtools \"-C 50 -f ".$index."\" > ".$type_output."/".$type.".mpileup");		
    system_call("samtools mpileup -d 7000 -I -f ".$index." ".$type_output."/".$type.".bam > ".$type_output."/".$type.".mpileup");
}		

sub fun_mutect{
  if ($pdx=~/mouse/) {$known=" ";} else {$known=" --cosmic ".$resource_cosmic." ";}
  system_call($java17." -Djava.io.tmpdir=".$mutect_tmp." -Xmx31g -jar ".$mutect." --analysis_type MuTect --reference_sequence ".$index.
    " --dbsnp ".$resource_dbsnp.$known." --input_file:tumor ".$tumor_bam." --input_file:normal ".$normal_bam.
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
 
  # lumpy
  #system_call("lumpyexpress -B ".$normal_bam.",".$tumor_bam." -R ".$index." -o ".$output."/lumpy_sv.vcf.txt -T ".$output."/normal/lumpy_temp");
}

sub fun_strelka_lofreq{
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

  # run lofreq
  system_call("lofreq call-parallel --pp-threads ".$thread." -s --sig 0.1 --bonf 1 -C 7 -f ".$index.
    " -S ".$resource_dbsnp." --call-indels -l ".$index.".exon.bed -o ".$output."/lofreq_t.vcf ".$tumor_bam);
  system_call("lofreq call-parallel --pp-threads ".$thread." -s --sig 1 --bonf 1 -C 7 -f ".$index.
    " -S ".$resource_dbsnp." --call-indels -l ".$index.".exon.bed -o ".$output."/lofreq_n.vcf ".$normal_bam); 
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
#/archive/BICF/shared/Kidney/exome/RAW/SAM9259827_1_1.gne.fastq.gz \
#/archive/BICF/shared/Kidney/exome/RAW/SAM9259827_1_2.gne.fastq.gz \
#/archive/BICF/shared/Kidney/exome/RAW/SAM19944142_1_1.gne.fastq.gz \
#/archive/BICF/shared/Kidney/exome/RAW/SAM19944142_1_2.gne.fastq.gz \
#32 hg38 /project/shared/xiao_wang/data/hg38/hs38d1.fa \
#/cm/shared/apps/java/oracle/jdk1.7.0_51/bin/java \
#/project/DPDS/Xiao_lab/shared/neoantigen/data/tmp/ human 0 \
#/project/shared/xiao_wang/software/disambiguate_pipeline
#
#perl /home2/twang6/software/cancer/somatic/somatic.pl \
#/archive/BICF/shared/Kidney/exome/RAW/LIB27355_SAM19944430_36695_R1.fastq.gz \
#/archive/BICF/shared/Kidney/exome/RAW/LIB27355_SAM19944430_36695_R2.fastq.gz \
#/archive/BICF/shared/Kidney/exome/RAW/LIB27324_SAM19944399_36728_R1.fastq.gz \
#/archive/BICF/shared/Kidney/exome/RAW/LIB27324_SAM19944399_36728_R2.fastq.gz \
#32 mm10 /project/shared/xiao_wang/data/mm10/mm10.fasta \
#/cm/shared/apps/java/oracle/jdk1.7.0_51/bin/java \
#/project/DPDS/Xiao_lab/shared/neoantigen/data/tmp mouse 0 \
#/project/shared/xiao_wang/software/disambiguate_pipeline
