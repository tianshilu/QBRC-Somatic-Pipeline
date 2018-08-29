#!/bin/bash
#SBATCH --job-name=example                                           # job name
#SBATCH --partition=256GB                                         # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                                    # number of nodes requested by user
#SBATCH --ntasks=32                                                  # number of total tasks
#SBATCH --time=10-00:00:00                                            # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=./sbatch_output_%j                                 # redirect both standard output and erro output to the same file
#SBATCH --error=./sbatch_error_%j
source /home2/twang6/.bash_profile

###########JOBSTART########################

perl /home2/twang6/software/immune/neoantigen/detect_neoantigen_aa.pl \
/home2/twang6/software/immune/mhc_i \
/home2/twang6/software/immune/mhc_ii \
~/iproject/test/example.fa \
"B_HLA-B*44:02_HLA-B*44:03" 2 \
~/iproject/test/neoantigen_aa

