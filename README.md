![QBRC logo](https://github.com/Somatic-pipeline/Somatic-pipeline/blob/master/QBRC.jpg)

# Somatic mutation calling pipeline
## Introduction
The somatic mutation calling pipeline of Wang lab identifies somatic and germline variants within whole exome sequencing (WXS), RNA sequencing and deep sequencing data. Mutations are identified by comparing allele frequencies in normal and tumor samples, annotating each mutation, and aggregating mutations from multiple cases into ine project file. The pipeline can be applied to fastq files or bam files for tumor tissue, normal tissue and patient derived xenografts.

Please refer to our paper for more detail of somatic mutation calling pipeline:["Neoantigen Clonal Balance Predicts Response to Checkpoint Inhibitor"](url peding)

# Dependencies
BWA (version >=0.7.15); STAR (required if applied for RNA sequencing data); sambamba; speedseq; varscan, samtools (version >=1.6); shimmer; annovar (database downloaded in default folder: refGene,ljb26_all,cosmic70,esp6500siv2_all,exac03,1000g2015aug); python (version 2); strelka (version >=2.8.3, note: strelka is tuned to run exome sequencing or RNA sequencing); manta (version >=1.4.0); java (version 1.8); perl (Parallel::ForkManager); lofreq_star (version >=2.1.3, for tumor-only calling); bowtie2 (version>= 2.3.4.3, for Patient Derived Xenograft models)

# Input files
Input can be fastq files or bam files. Input can also be a mixture of fastq adn bam files.

# Main procedures:
## Genome Alignment
Genome sequencing files are aligned to the human or mouse reference genome by BWA-MEM (Please contact Tianshi.Lu@UTSouthwestern.edu for genome reference files).
## Alignment Co-Cleaning
Picard was used to add read group information and sambamba was used to mark PCR duplilcates. GATK toolkt was used to perform base quality score relcalibration adn local realignment around Indels.
## Variant Calling
MuTect, VarScan Shimmer, SpeedSeq, Manta, and Strelka2 were used to call SNPs and Indels. A mutation that was repeatedly called by any two of these softwares was retained.
## Mutation Annotation
Annovar was used to annotate SNPs, and Indels and protein sequence changes. Somatic mutations and germline mutations were annotated according to the mutation allele frequencies in the normal and tumor samples.
## Filter False Mutations
All SNPs and Indels were combined ony kept if there were at least 7 total( wild type adn variant) reads in the normal sample and at least 3 variant reads in the tumor sample.

# Guided Tutorial
## somatic.pl
The code for somatic and germline mutation calling for a pair of normal and tumor sequencing files.
### Command
perl /Directory/to/folder/of/code/somatic.pl \
sequencing_file_1 \
sequencing_file_2 \
sequencing_file_3 \
sequencing_file_4 \
thread build index java17 /Directory/to/output pdx 

Example: \
perl ~/somatic/somatic.pl ~/seq/1799-01N.R1.fastq.gz ~/seq/1799-01N.R2.fastq.gz ~/seq/1799-01T.R1.fastq.gz ~/seq/1799-01T.R2.fastq.gz 32 hg38 ~/ref/hg38/hs38d1.fa /cm/shared/apps/java/oracle/jdk1.7.0_51/bin/java ~/somatic_result/1799-01/ human

#### Note: 
Input seuqencing files: (1) If input are fastq files, they must be 'gz' files. 'sequencing_file_1', 'sequencing_file_2' are path to fastq1 and fastq2 of normal sample; 'sequencing_file_3', 'sequencing_file_4' are path to fastq1 and fastq2 of tumor samples. (2) If input are bam files, use "bam /path/to/bam/files.bam" in replace of the tow corresponding fastq input files. (3) If input are RNA sequencing files, use "RNA:fastq1" or "RNA:bam" at the first or third slot. (4) If input are deep exome sequencing data, use "Deep:fastq1" at the first or third slot. (5) For tumor-only calling, put "NA NA" in the first two slots. Results will be written to *germline* output files. (6) Optional: run somatic_script/SurecallTrummer.jar on the fastq files before runnign somatic.pl for deep seuquencing files. \
"thread": number of threads to use. Recommended: 32 \
"build": genome build, hg19 or hg38 or mm10. \
"index": path (including file names) to the reference genome. \
"java17": path (including the executable file name) to java 1.7 (needed only for MuTect). \
"ouput": the output folder, it will be deleted (if pre-existing) adn re-created during analysis. \
"pdx": "PDX" or "human" if this is PDX sample, reads will be aligned to mouse genome first. And unmapped reads will be mapped to the human genome. 

## job_somatic.pl
BioHPC submission wrapper for somatic.pl for a batch of sampels.
### Command
perl /Directory/to/folder/of/code/job_somatic.pl \
design.txt \
example_file \
thread build index java17 n

somatic_design.txt example:
~/seq/1799-01N.R1.fastq.gz ~/seq/1799-01N.R2.fastq.gz ~/seq/1799-01T.R1.fastq.gz ~/seq/1799-01T.R2.fastq.gz ~/out/1799-01/ \
~/seq/1799-02N.R1.fastq.gz ~/seq/1799-02N.R2.fastq.gz ~/seq/1799-02T.R1.fastq.gz ~/seq/1799-02T.R2.fastq.gz ~/out/1799-02/ \
~/seq/1799-03N.R1.fastq.gz ~/seq/1799-03N.R2.fastq.gz ~/seq/1799-03T.R1.fastq.gz ~/seq/1799-03T.R2.fastq.gz ~/out/1799-03/ 

Command example:  \
perl ~/somatic/job_somatic.pl somatic_design.txt ~/somatic/example/example.sh 32 hg38 ~/ref/hg38/hs38d1.fa /cm/shared/apps/java/oracle/jdk1.7.0_51/bin/java 2

#### Note:
"design.txt" is the batch job design file. It has 6 columns separated by '\t', the first four slots are fastq files or bam files for normal and tumor samples. The fifth is the output folder, and the last is "PDX" or "". \
"example_file" is the demo job submission shell script. A default one is in example/. \
"thread": number of threads to use. Recommended: 32 \
"build": genome build, hg19 or hg38 or mm10. \
"index": path (including file names) to the reference genome. \
"java17": path (including the executable file name) to java 1.7 (needed only for MuTect). \
"n": bundle $n somatic calling job into one submission. 

## filter.R
Post-processing script for somatic mutations for a batch of sampels. \
### Command
Rscript filter.R \
design.txt \
output build index VAF_cutoff filter 
#### Note:
"design.txt": tab-delimited file with three columns: sample_id, patient_id, output folder. \
"output": the output folder to place all filtering results. \
"build": the reference genome build, hg38, hg19 etc. \
"index": the path to the reference genome file. \
"VAF_cutoff": the minimum VAF of the mutations in the tumor sample (recommended: 0.001-0.05). \
"filter": TRUE or FALSE. Whether to filter out extremely long genes in the list "TTN","KCNQ1OT1","MUC16","ANKRD20A9P","TSIX","SYNE1","ZBTB20","OBSCN", "SH3TC2","NEB","MUC19","MUC4","NEAT1","SYNE2","CCDC168","AAK1","HYDIN","RNF213","LOC100131257","FSIP2". These genes usually turn out ot have somatic muitations in any cohort of patients. Default is FALSE.

filter_design.txt example: \
1799-01 pat-01 ~/filter/1799-01/ \
1799-02 pat-02 ~/filter/1799-02/ \
1799-03 pat-03 ~/filter/1799-03/ 

Command example: \
Rscript ~/somatic/filter.R filter_design.txt ~/filter/ hg38 ~/ref/hg38/hs38d1.fa 0.01 FALSE 

## cnv.pl
Pipeline for somatic copy number variation calling and quality check for each sample 
### Command 
perl cnv.pl \
sequencing_file_1 \
sequencing_file_2 \
sequencing_file_3 \
sequencing_file_4 \
thread index somatic_mutation_result output 
##### Note:
prerequisite in path: R; BWA; sambamba; perl (Parallel::ForkManager); samtools (version>=1.6); cnvkit; fastqc
Input seuqencing files: (1) If input are fastq files, they must be 'gz' files. 'sequencing_file_1', 'sequencing_file_2' are path to fastq1 and fastq2 of normal sample; 'sequencing_file_3', 'sequencing_file_4' are path to fastq1 and fastq2 of tumor samples. (2) If input are bam files, use "bam /path/to/bam/files.bam" in replace of the tow corresponding fastq input files. \
"thread": number of threads to use. Recommended: 32 \
"index": the path to the reference genome file. \
"somatic_mutation_result": somatic mutation calling output file. THis is for adjusting CNV by somatic mutation VAF. Set to 1 to turn off this adjustment. \
"output":the output folder. it will be deleted (if pre-existing) adn re-created during analysis. \
The CNV calling needs at least 128GB of memory. 

Example: \
perl ~/somatic/cnv.pl ~/seq/1799-01N.R1.fastq.gz ~/seq/1799-01N.R2.fastq.gz ~/seq/1799-01T.R1.fastq.gz ~/seq/1799-01T.R2.fastq.gz 32 ~/ref/hg38/hs38d1.fa ~/somatic_result/1799-01/somatic_mutation_hg38.txt ~/cnv_result/1799-01

## job_cnv.pl: 
BioHPC submission wrapper for cnv.pl for a batch of samples.
### Command
perl job_cnv.pl \
design.txt \
example.sh thread index n

cnv_design.txt example: \
~/seq/1799-01N.R1.fastq.gz ~/seq/1799-01N.R2.fastq.gz ~/seq/1799-01T.R1.fastq.gz ~/seq/1799-01T.R2.fastq.gz ~/somatic_result/1799-01/somatic_mutations_hg38.txt ~/cnv_result/1799-01/ \
~/seq/1799-02N.R1.fastq.gz ~/seq/1799-02N.R2.fastq.gz ~/seq/1799-02T.R1.fastq.gz ~/seq/1799-02T.R2.fastq.gz ~/somatic_result/1799-02/somatic_mutations_hg38.txt ~/cnv_result/1799-02/ \
~/seq/1799-03N.R1.fastq.gz ~/seq/1799-03N.R2.fastq.gz ~/seq/1799-03T.R1.fastq.gz ~/seq/1799-03T.R2.fastq.gz ~/somatic_result/1799-03/somatic_mutations_hg38.txt ~/cnv_result/1799-03/ 

Command example: \
perl ~/somatic/job_cnv.pl cnv_design.txt ~/somatic/example/example.sh 32 ~/ref/hg38/hs38d1.fa 2 

## summarize_cnv.R
Summarizing script for CNV and quality check callings for a batch of samples.
### Command
Rscript summarize_cnv.R design.txt output index

cnv_sum_design.txt example: \
sample_id folder \
1799-01 ~/cnv_result/1799-01 \
1799-02 ~/cnv_result/1799-02 \
1799-03 ~/cnv_result/1799-03 

Command example: \
Rscript ~/somatic/summarize_cnv.R cnv_sum_design.txt ~/cnv_sum/ ~/ref/hg38/
