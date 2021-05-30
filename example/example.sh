#!/bin/bash
#SBATCH --job-name=c
#SBATCH --partition=256GB
#SBATCH --nodes=1
#SBATCH --time=60-00:00:00
#SBATCH --output=./sbatch_output_%j
#SBATCH --error=./sbatch_error_%j
source ~/.bash_profile
##########JOBSTART###########################
perl job_somatic.pl /project/job_somatic.txt example.sh 32 hg38 /project/data/hg38/hs38d1.fa /cm/shared/apps/java/oracle/jdk1.7.0_51/bin/java 3
