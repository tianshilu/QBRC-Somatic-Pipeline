# Somatic-pipeline
Somatic mutation calling pipeline
#To use somatic.pl, please use following instructions:
prerequisite in path: 
Rscript, bwa (>=0.7.15), STAR (if using RNA-Seq data), sambamba, speedseq, varscan, samtools (>=1.6), shimmer
annovar (database downloaded in default folder: refGene,ljb26_all,cosmic70,esp6500siv2_all,exac03,1000g2015aug)
python (2), strelka (>=2.8.3, note: strelka is tuned to run exome-seq or RNA-seq), manta(>=1.2.0), java (1.8)
perl (need Parallel::ForkManager)
input format: 
fastq files: (1) fastq1 and fastq2 of normal sample, fastq1 and fastq2 of tumor sample
                  default input is full path to the 4 fastq files for tumor and normal samples 
              (2) if need to directly input bam files, use "bam path_to_bam.bam" in replace of the two corresponding fastq input files
                  can be a mixture of fastq and bam input
              (3) if RNA-Seq data are used, use "RNA:fastq1" or "RNA:bam" at the first or third slot
              (4) if Agilent SureSelect (Deep exome sequencing) data, use "Deep:fastq1" at the first or third slot
                  optional: run somatic_script/SurecallTrimmer.jar on the fastq files before running somatic.pl
 thread: number of threads to use. Recommended: 32
 build: human genome build, hg19 or hg38
 index: path (including file name) to the human reference genome
 java17: path (including the executable file name) to java 1.7 (needed only for mutect)
 output: the output folder, it will be deleted (if pre-existing) and re-created during analysis
 need at least 128GB of memory
#To use job_somatic.pl, please use following instructions:
input format:
jobs: the batch job design file, it has 5 columns separated by \t, the first four are fastq files or bam files (normal+tumor), 
      the last one is the output folder. Commented lines ("#" at the front) are skipped
 example: the demo job submission shell script. A default one is in this folder
thread,build,index,java17: follow those in somatic.pl
n: bundle $n somatic calling jobs into one submission
