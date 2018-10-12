![QBRC logo](https://github.com/Somatic-pipeline/Somatic-pipeline/blob/master/QBRC.jpg)

# Somatic mutation calling pipeline
## Introduction
The somatic mutation calling pipeline of Wang lab can be applied to call somatic and germline mutations for exome-seq, RNA-seq and ultra deep sequencing. 
Input can be fastq files or bam files.

somatic.pl: pipeline for somatic and germline mutation calling (for each sample)
job_somatic.pl: biohpc submission wrapper for somatic.pl (for a batch of samples)
filter.R: post-processing script for somatic mutations (for a batch of samples)

cnv.pl: pipeline for somatic CNV calling and quality check (for each sample)
job_cnv.pl: biohpc submission wrapper for cnv.pl (for a batch of samples)
summarize_cnv.R: summarizing script for CNV and quality check callings (for a batch of samples)
