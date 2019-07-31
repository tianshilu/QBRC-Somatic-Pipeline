# The QBRC somatic mutation calling pipeline
![preview](https://github.com/tianshilu/QBRC-Somatic-Pipeline/blob/master/somatic_flow.jpg)
## Introduction
The QBRC mutation calling pipeline is a flexible and comprehensive pipeline for mutation calling that has glued together a lot of commonly used software and data processing steps for mutation calling. The mutation calling software include: sambamba, speedseq, varscan, shimmer, strelka, manta, lofreq_tar. It identifies somatic and germline variants from whole exome sequencing (WXS), RNA sequencing and deep sequencing data. It can be used for human, PDX, and mouse data (fastq files or bam files as input). Please refer to https://qbrc.swmed.edu/labs/wanglab/index.php for more information. If you used our pipeline in your publication, please cite our paper ["Neoantigen Clonal Balance Predicts Response to Checkpoint Inhibitor"] (under review) and please refer to https://github.com/Somatic-pipeline/QBRC-Somatic-Pipeline/blob/master/Liscense.txt for liscence information.
## Running time
For a paird of 'fastq.gz'files of 200M, it takes around 2 hours to finish somatic mutation calling.
## Hardwares/Softwares Dependencies
64 bit linux operating system  
BWA (version >=0.7.15)  
STAR (required if applied for RNA sequencing data)  
sambamba  
speedseq  
varscan  
samtools (version >=1.6)  
shimmer  
annovar (database downloaded in default folder: refGene,ljb26_all,cosmic70,esp6500siv2_all,exac03,1000g2015aug)  
python2  
strelka (version >=2.8.3, note: strelka is tuned to run exome sequencing or RNA sequencing)  
manta (version >=1.4.0)  
java (version 1.8)  
perl (Parallel::ForkManager)  
lofreq_star (version >=2.1.3, for tumor-only calling)  
bowtie2 (version>= 2.3.4.3, for Patient Derived Xenograft models)

## Input files
Input can be fastq files or bam files or a mixture of fastq and bam files.
## Main procedures:
* Genome Alignment:  
Genome sequencing files are aligned to the human reference genome by BWA-MEM (Please contact Tianshi.Lu@UTSouthwestern.edu for genome reference files). Picard was used to add read group information and sambamba was used to mark PCR duplilcates. GATK toolkt was used to perform base quality score relcalibration adn local realignment around Indels.
* Variant Calling:  
MuTect, VarScan Shimmer, SpeedSeq, Manta, and Strelka2 were used to call SNPs and Indels. A mutation that was repeatedly called by any two of these softwares was retained.
* Mutation Annotation:  
Annovar was used to annotate SNPs, and Indels and protein sequence changes. Somatic mutations and germline mutations were annotated according to the mutation allele frequencies in the normal and tumor samples.
* Filter False Mutations:  
All SNPs and Indels were combined ony kept if there were at least 7 total( wild type and variant) reads in the normal sample and at least 3 variant reads in the tumor sample. Variants with allele frequency more than 2 times allele frequency of the according normal allele are kept. Variants with allele frequency less than 5% in background sample are kept.
## Guided Tutorial
## somatic.pl
The code for somatic and germline mutation calling for a pair of normal and tumor sequencing files.

### Usage
```
perl /Path/to/somatic.pl <normal_fastq1> <normal_fastq2/NA> <tumor_fastq1> <tumor_fastq2/NA> \ 
<thread> <build> <index> <java17> </Path/to/output> <pdx>
```
* fastq files:  
  * fastq1 and fastq2 of normal sample, fastq1 and fastq2 of tumor sample (must be .gz) default input is full path to the 4 fastq files for tumor and normal samples.  
  * If need to directly input bam files, use "bam path_to_bam.bam" in replace of the two corresponding fastq input files can be a mixture of fastq and bam input.  
  * If RNA-Seq data are used, use "RNA:fastq1" or "RNA:bam" at the first
  * If Agilent SureSelect (Deep exome sequencing) data, use "Deep:fastq1" at the first or third slot.  
 optional: run somatic_script/SurecallTrimmer.jar on the fastq files before running somatic.pl.  
  * For tumor-only calling, put "NA NA" in the slots of the normal samples. Results will be written to *germline* files
  * If only single end fastq data are available, put the fastq file(s) at the first and/or the third slots, then put NA in the second and/or fourth slot.  
* thread: number of threads to use. Recommended: 32  
* build: hg19 or hg38 or mm10  
* index: path (including file name) to the human/mouse reference genome  
* java17: path (including the executable file name) to java 1.7 (needed only for mutect)  
* output: the output folder, it will be deleted (if pre-existing) and re-created during analysis  
* pdx: "PDX" or "human" or "mouse"  

### Example:
```
perl ~/somatic/somatic.pl ~/seq/1799-01N.R1.fastq.gz ~/seq/1799-01N.R2.fastq.gz ~/seq/1799-01T.R1.fastq.gz ~/seq/1799-01T.R2.fastq.gz 32 hg38 ~/ref/hg38/hs38d1.fa /cm/shared/apps/java/oracle/jdk1.7.0_51/bin/java ~/somatic_result/1799-01/ human
```
### Note:
Input seuqencing files:  
(1) If input are fastq files, they must be 'gz' files. 'sequencing_file_1', 'sequencing_file_2' are path to fastq1 and fastq2 of normal sample; 'sequencing_file_3', 'sequencing_file_4' are path to fastq1 and fastq2 of tumor samples.  
(2) If input are bam files, use "bam /path/to/bam/files.bam" in replace of the tow corresponding fastq input files.  
(3) If input are RNA sequencing files, use "RNA:fastq1" or "RNA:bam" at the first or third slot.  
(4) If input are deep exome sequencing data, use "Deep:fastq1" at the first or third slot.  
(5) For tumor-only calling, put "NA NA" in the first two slots. Results will be written to germline output files.  
(6) Optional: run somatic_script/SurecallTrummer.jar on the fastq files before runnign somatic.pl for deep seuquencing files.  
(7) If only single end fastq data are available, put the fastq file(s) at the first and/or the third slots, then put NA in the second and/or fourth slot.  
* thread : number of threads to use.  
* build : genome build, hg19 or hg38.  
* index : path (including file names) to the reference genome fasta file of the reference bundle hg38 or hg19. (The pipeline will search for other files in that bundle folder automatically.)  
* java17 : path (including the executable file name) to java 1.7 (needed only for MuTect).  
* output : the output folder, it will be deleted (if pre-existing) adn re-created during analysis.  
* pdx : "PDX" or "human" if this is PDX sample, reads will be aligned to mouse genome first. And unmapped reads will be mapped to the human genome.  
(4) Example data to run the QBRC somatic mutation pipeline can be found at https://github.com/Somatic-pipeline/QBRC-Somatic-Pipeline/tree/master/example/example_dataset/sequencing. The output for the example data can be found at https://github.com/Somatic-pipeline/QBRC-Somatic-Pipeline/tree/master/example/example_dataset/example_output.
## job_somatic.pl
Slurm wrapper for somatic.pl for a batch of sampels and it is easy to change for other job scheduler system by revising this line of code: "system("sbatch ".$job)" and using proper demo job submission shell script.
### Command
```
perl /Directory/to/folder/of/code/job_somatic.pl 
design.txt 
example_file 
thread build index java17 n\
somatic_design.txt example (5 columns; columns seperated by tab):\
~/seq/1799-01N.R1.fastq.gz ~/seq/1799-01N.R2.fastq.gz ~/seq/1799-01T.R1.fastq.gz ~/seq/1799-01T.R2.fastq.gz ~/out/1799-01/ human \
~/seq/1799-02N.R1.fastq.gz ~/seq/1799-02N.R2.fastq.gz ~/seq/1799-02T.R1.fastq.gz ~/seq/1799-02T.R2.fastq.gz ~/out/1799-02/ human \
~/seq/1799-03N.R1.fastq.gz ~/seq/1799-03N.R2.fastq.gz ~/seq/1799-03T.R1.fastq.gz ~/seq/1799-03T.R2.fastq.gz ~/out/1799-03/ human 
```
### Command example
```
perl ~/somatic/job_somatic.pl somatic_design.txt ~/somatic/example/example.sh 32 hg38 ~/ref/hg38/hs38d1.fa /cm/shared/apps/java/oracle/jdk1.7.0_51/bin/java 2
```
### Note:
* design.txt: the batch job design file. It has 6 columns separated by '\t', the first four slots are fastq files or bam files for normal and tumor samples. The fifth is the output folder, and the last is "PDX" or "human". 
* example_file : the demo job submission shell script. A default one is in example/. 
* thread : number of threads to use. Recommended: 32 
* build : genome build, hg19 or hg38. 
* index : path (including file names) to the reference genome in the reference bundle. 
* java17 : path (including the executable file name) to java 1.7 (needed only for MuTect). 
* n : bundle $n somatic calling job into one submission.
## filter.R
Post-processing script for somatic mutations for a batch of sampels.
### Usage
```
Rscript filter.R  design.txt output build index VAF_cutoff filter
```
### Note:
* design.txt : tab-delimited file with three columns: sample_id, patient_id, output folder. 
* output : the output folder to place all filtering results. 
* build : the reference genome build, hg38, hg19 etc. 
* index : the path to the reference genome file in the reference bundle. 
* VAF_cutoff : the minimum VAF of the mutations in the tumor sample (recommended: 0.001-0.05). 
* filter : TRUE or FALSE. Whether to filter out extremely long genes in the list "TTN","KCNQ1OT1","MUC16","ANKRD20A9P","TSIX","SYNE1","ZBTB20","OBSCN", "SH3TC2","NEB","MUC19","MUC4","NEAT1","SYNE2","CCDC168","AAK1","HYDIN","RNF213","LOC100131257","FSIP2". These genes usually turn out ot have somatic muitations in any cohort of patients. Default is FALSE.\
## filter_design.txt example (3 columns; columns seperated by tab; header): 
sample_id patient_id folder 
1799-01 pat-01 ~/filter/1799-01/ 
1799-02 pat-02 ~/filter/1799-02/ 
1799-03 pat-03 ~/filter/1799-03/
### Command example: 
```
Rscript ~/somatic/filter.R filter_design.txt ~/filter/ hg38 ~/ref/hg38/hs38d1.fa 0.01 FALSE
```

## cnv.pl
Pipeline for somatic copy number variation calling and quality check for each sample
### Command
```
perl cnv.pl 
sequencing_file_1 
sequencing_file_2 
sequencing_file_3 
sequencing_file_4 
thread index somatic_mutation_result output
```
### Note:
* prerequisite in path: R; BWA; sambamba; perl (Parallel::ForkManager); samtools (version>=1.6); cnvkit; 
* fastqc Input seuqencing files: 
  * (1) If input are fastq files, they must be 'gz' files. 
  * 'sequencing_file_1', 'sequencing_file_2' are path to fastq1 and fastq2 of normal sample; 
  * 'sequencing_file_3', 'sequencing_file_4' are path to fastq1 and fastq2 of tumor samples. 
  * (2) If input are bam files, use "bam /path/to/bam/files.bam" in replace of the tow corresponding fastq input files. 
* thread : number of threads to use. Recommended: 32 
* index : the path to the reference genome file in the reference bundle. 
* somatic_mutation_result": somatic mutation calling output file. THis is for adjusting CNV by somatic mutation VAF. Set to 1 to turn off this adjustment. 
* output :the output folder. it will be deleted (if pre-existing) adn re-created during analysis. 
The CNV calling needs at least 128GB of memory.
### Example
```
perl ~/somatic/cnv.pl ~/seq/1799-01N.R1.fastq.gz ~/seq/1799-01N.R2.fastq.gz ~/seq/1799-01T.R1.fastq.gz ~/seq/1799-01T.R2.fastq.gz 32 ~/ref/hg38/hs38d1.fa ~/somatic_result/1799-01/somatic_mutation_hg38.txt ~/cnv_result/1799-01
```
## job_cnv.pl:
Slurm wrapper for cnv.pl for a batch of samples and it is easy to change for other job scheduler system by revising this line of code: "system("sbatch ".$job)" and using proper demo job submission shell script.
### Command
```
perl job_cnv.pl 
design.txt 
example.sh thread index n\
cnv_design.txt example (6 columns; columns seperated by tab): 
~/seq/1799-01N.R1.fastq.gz ~/seq/1799-01N.R2.fastq.gz ~/seq/1799-01T.R1.fastq.gz ~/seq/1799-01T.R2.fastq.gz ~/somatic_result/1799-01/somatic_mutations_hg38.txt ~/cnv_result/1799-01/ 
~/seq/1799-02N.R1.fastq.gz ~/seq/1799-02N.R2.fastq.gz ~/seq/1799-02T.R1.fastq.gz ~/seq/1799-02T.R2.fastq.gz ~/somatic_result/1799-02/somatic_mutations_hg38.txt ~/cnv_result/1799-02/ 
~/seq/1799-03N.R1.fastq.gz ~/seq/1799-03N.R2.fastq.gz ~/seq/1799-03T.R1.fastq.gz ~/seq/1799-03T.R2.fastq.gz ~/somatic_result/1799-03/somatic_mutations_hg38.txt ~/cnv_result/1799-03/
```
### Command example
```
perl ~/somatic/job_cnv.pl cnv_design.txt ~/somatic/example/example.sh 32 ~/ref/hg38/hs38d1.fa 2
```
## summarize_cnv.R
Summarizing script for CNV and quality check callings for a batch of samples.
### Command
```
Rscript summarize_cnv.R design.txt output index
cnv_sum_design.txt example (2 columns; columns seperated by tab; header): 
sample_id folder 
1799-01 ~/cnv_result/1799-01 
1799-02 ~/cnv_result/1799-02 
1799-03 ~/cnv_result/1799-03
```
### Command example: 
```
Rscript ~/somatic/summarize_cnv.R cnv_sum_design.txt ~/cnv_sum/ ~/ref/hg38/
```
