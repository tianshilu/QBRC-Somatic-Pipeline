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
Genome sequencing files are aligned to the human or mouse reference genome by BWA-MEM. 
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
#### Note: 
Input seuqencing files: (1) If input are fastq files, they must be 'gz' files. 'sequencing_file_1', 'sequencing_file_2' are path to fastq1 and fastq2 of normal sample; 'sequencing_file_3', 'sequencing_file_4' are path to fastq1 and fastq2 of tumor samples. (2) If input are bam files, use "bam /path/to/bam/files.bam" in replace of the tow corresponding fastq input files. (3) If input are RNA sequencing files, use "RNA:fastq1" or "RNA:bam" at the first or third slot. (4) If input are deep exome sequencing data, use "Deep:fastq1" at the first or third slot. (5) For tumor-only calling, put "NA NA" in the first two slots. Results will be written to *germline* output files. (6) Optional: run somatic_script/SurecallTrummer.jar on the fastq files before runnign somatic.pl for deep seuquencing files.
"thread": number of threads to use. Recommended: 32
"build": genome build, hg19 or hg38 or mm10.
"index": path (including file names) to the reference genome.
"java17": path (including the executable file name) to java 1.7 (needed only for MuTect).
"ouput": the output folder, it will be deleted (if pre-existing) adn re-created during analysis.
"pdx": "PDX" or "human" if this is PDX sample, reads will be aligned to mouse genome first. And unmapped reads will be mapped to the human genome.

## job_somatic.pl
BioHPC submission wrapper for somatic.pl for a batch of sampels.
### Command
perl /Directory/to/folder/of/code/job_somatic.pl \
design.txt \
example_file \
thread build index java17 n
#### Note:
"design.txt" is the batch job design file. It has 6 columns separated by '\t', the first four slots are fastq files or bam files for normal and tumor samples. The fifth is the output folder, and the last is "PDX" or "".
"example_file" is the demo job submission shell script. A default one is in example/.
"thread": number of threads to use. Recommended: 32
"build": genome build, hg19 or hg38 or mm10.
"index": path (including file names) to the reference genome.
"java17": path (including the executable file name) to java 1.7 (needed only for MuTect).
"n": bundle $n somatic calling job into one submission.


filter.R: post-processing script for somatic mutations (for a batch of samples)
cnv.pl: pipeline for somatic CNV calling and quality check (for each sample)
job_cnv.pl: biohpc submission wrapper for cnv.pl (for a batch of samples)
summarize_cnv.R: summarizing script for CNV and quality check callings (for a batch of samples)
