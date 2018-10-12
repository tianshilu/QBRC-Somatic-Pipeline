![QBRC logo](https://github.com/Somatic-pipeline/Somatic-pipeline/blob/master/QBRC.jpg)

# Somatic mutation calling pipeline
## Introduction
The somatic mutation calling pipeline of Wang lab identifies somatic and germline variants within whole exome sequencing (WXS), RNA sequencing and ultra deep sequencing data. Mutations are identified by comparing allele frequencies in normal and tumor samples, annotating each mutation, and aggregating mutations from multiple cases into ine project file.

Please refer to our paper for more detail of somatic mutation calling pipeline:["Neoantigen Clonal Balance Predicts Response to Checkpoint Inhibitor"](url peding)

# Dependencies:
BWA (version >=0.7.15); STAR (required if applied for RNA sequencing data); sambamba; speedseq; varscan, samtools (version >=1.6); shimmer; annovar (database downloaded in default folder: refGene,ljb26_all,cosmic70,esp6500siv2_all,exac03,1000g2015aug); python (version=2); strelka (version >=2.8.3, note: strelka is tuned to run exome sequencing or RNA sequencing); manta (>=
# Main procedures:

Input can be fastq files or bam files.

somatic.pl: pipeline for somatic and germline mutation calling (for each sample)
job_somatic.pl: biohpc submission wrapper for somatic.pl (for a batch of samples)
filter.R: post-processing script for somatic mutations (for a batch of samples)

cnv.pl: pipeline for somatic CNV calling and quality check (for each sample)
job_cnv.pl: biohpc submission wrapper for cnv.pl (for a batch of samples)
summarize_cnv.R: summarizing script for CNV and quality check callings (for a batch of samples)
