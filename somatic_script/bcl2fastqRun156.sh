#!/bin/bash

#SBATCH --job-name=bcl2fastq 

#SBATCH --partition=super

#SBATCH --nodes=1

#SBATCH --time=2-10:00:00

#SBATCH --output=bcl2fastq.%j.out 

#SBATCH --error=bcl2fastq.%j.err


module load bcl2fastq/2.17.1.14  


bcl2fastq --runfolder-dir /project/CRI/Zhu_lab/shared/NGS/Raw/Run156_HZ_11-24-17/171124_NS500717_0197_AHW23YBGX3 --output-dir /project/shared/zhu_bicf/Run156/fastq --sample-sheet /project/shared/zhu_bicf/Run156/SampleSheetXTHS.csv --use-bases-mask Y148,I8,Y10,Y148 --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0 
