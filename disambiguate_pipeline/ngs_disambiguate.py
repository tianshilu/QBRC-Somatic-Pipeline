# -*- coding: utf-8 -*-
'''
A fastq2fastq script aimed to separate mixed human and mouse sequences from a 
hybrid sample such as PDX. 

Requirements
--------
This script need the standard QBRC somatic pipeline environment, plus 
an additional dependency: `ngs-disambiguate` from conda. 
Install it with 
```
conda install ngs-disambiguate
```

Usage
--------
usage: ngs_disambiguate.py [-h] [-o [OUTPUT_FOLDER]] [-a [ALIGNER]]
                           [-r [REFERENCE]] [-i [INTERMEDIATE_PATH]]
                           [-b [BAM2FASTQ]] [-k] [-n [NUMBER_THREADS]]
                           fastq1 fastq2

Wrapper function for running disambiguate. Need to load ngs-disambiguate conda
env first.

positional arguments:
  fastq1                Fastq read 1
  fastq2                Fastq read 2

optional arguments:
  -h, --help            show this help message and exit
  -o [OUTPUT_FOLDER], --output_folder [OUTPUT_FOLDER]
                        disambiguated output (default: ./disambiguated)
  -a [ALIGNER], --aligner [ALIGNER]
                        The aligner used to generate these reads. Some
                        aligners set different tags. Select from {bwa,star}.
                        (default: star)
  -r [REFERENCE], --reference [REFERENCE]
                        Path to genome reference. Skip this if using default
                        refs. if entered it should be in this format
                        [human_ref|mouse_ref] (default: default)
  -i [INTERMEDIATE_PATH], --intermediate_path [INTERMEDIATE_PATH]
                        Path to store intermediate files (default: ./tmp)
  -b [BAM2FASTQ], --bam2fastq [BAM2FASTQ]
                        Path to bam2fastq script. (default: /project/shared/xi
                        ao_wang/software/somatic/somatic_script/bam2fastq.pl)
  -k, --keep_intermedaite
                        Whether to keep intermediate files. (default: False)
  -n [NUMBER_THREADS], --number_threads [NUMBER_THREADS]
                        Number of threads to use. (default: 32)

Output
--------
[OUTPUT_FOLDER]
    fastq1.fastq.human
    fastq2.fastq.human
    fastq1.fastq.mouse
    fastq2.fastq.mouse

Example Usage
--------
conda activate [ngs_disambiguate_env]
python ngs_disambiguate.py example_data/hs1m_mm1m_R1.fastq.gz example_data/hs1m_mm1m_R2.fastq.gz
'''

import os
import argparse

if __name__ == "__main__":
    # Default references
    ref_dict = {
        'human':{
            'star':'~/work_personal/ref/hg38/STAR',
            'bwa':'/home2/twang6/data/genomes/hg38/hs38d1.fa',
        },
        'mouse':{
            'star':'~/work_personal/ref/mm10/star',
            'bwa':'/home2/twang6/data/genomes/mm10/mm10.fasta',
        }
    }

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Wrapper function for running disambiguate. \
            Need to load ngs-disambiguate conda env first.'
    )
    parser.add_argument(
        "fastq1", type=str, help="Fastq read 1")
    parser.add_argument(
        "fastq2", type=str, help="Fastq read 2")
    parser.add_argument(
        "-o", "--output_folder", type=str, nargs='?',default='./disambiguated', 
        help="disambiguated output")
    parser.add_argument(
        "-a", "--aligner", type=str, nargs='?',default='star',
        help="The aligner used to generate these reads. Some aligners set different tags. \
            Select from {bwa,star}.")
    parser.add_argument(
        "-r", "--reference", type=str, nargs='?',default='default',
        help="Path to genome reference. Skip this if using default refs. \
            if entered it should be in this format [human_ref|mouse_ref]")
    parser.add_argument(
        "-i", "--intermediate_path", type=str, nargs='?', default='./tmp',
        help="Path to store intermediate files")
    parser.add_argument(
        "-b", "--bam2fastq", type=str, nargs='?',
        default='/project/shared/xiao_wang/software/somatic/somatic_script/bam2fastq.pl',
        help="Path to bam2fastq script.")
    parser.add_argument(
        "-k", "--keep_intermedaite", default=False,action='store_true',
        help="Whether to keep intermediate files.")
    parser.add_argument(
        "-n", "--number_threads", type=int, nargs='?', default=32,
        help="Number of threads to use.")
    # Parse args
    args = parser.parse_args()
    fastq1 = args.fastq1
    fastq2 = args.fastq2
    output = args.output_folder
    aligner = args.aligner
    int_path = args.intermediate_path
    ref = args.reference
    keep = args.keep_intermedaite
    bam2fastq = args.bam2fastq
    nt = args.number_threads
    if ref == 'default':
        ref_path_h = ref_dict['human'][aligner]
        ref_path_m = ref_dict['mouse'][aligner]
    else:
        ref_path_h, ref_path_m = ref.split('|')

    # make intermediate file path.
    alignment_human_path = '/'.join([int_path,'alignment_human/'])
    alignment_mouse_path = '/'.join([int_path,'alignment_mouse/'])
    if not os.path.exists(int_path):
        os.mkdir(int_path)
    for f in [alignment_human_path, alignment_mouse_path]:
        if not os.path.exists(f):
            os.mkdir(f)
    
    # Alignment step
    if aligner == 'star':     
        os.system(
            "STAR --runThreadN {} --genomeDir {} --readFilesIn {} {} \
            --outFileNamePrefix {} --outSAMtype BAM Unsorted \
            --outSAMattributes NM".format(
                nt, ref_path_h, fastq1, fastq2, alignment_human_path
            )
        ) # human
        os.system(
            "STAR --runThreadN {} --genomeDir {} --readFilesIn {} {} \
            --outFileNamePrefix {} --outSAMtype BAM Unsorted \
            --outSAMattributes NM".format(
                nt, ref_path_m, fastq1, fastq2, alignment_mouse_path
            )
        ) # mouse
        bam_h = alignment_human_path + 'Aligned.out.bam'
        bam_m = alignment_mouse_path + 'Aligned.out.bam'
    else:
        samfile_path_h = alignment_human_path + "alignment.sam"
        samfile_path_m = alignment_mouse_path + "alignment.sam"
        bam_h = alignment_human_path + "alignment.bam"
        bam_m = alignment_mouse_path + "alignment.bam"
        os.system(
            "bwa mem -v 1 -t {} -a -M {} {} {} > {}".format(
                nt, ref_path_h, fastq1, fastq2, samfile_path_h
            )
        ) # human

        os.system(
            "bwa mem -v 1 -t {} -a -M {} {} {} > {}".format(
                nt, ref_path_m, fastq1, fastq2, samfile_path_m
            )
        ) # mouse

        os.system(
            "sambamba view -f bam -h -v -S -l 0 -t {} {} -o {}".format(
                nt,samfile_path_h, bam_h
            )
        ) # human bam
        
        os.system(
            "sambamba view -f bam -h -v -S -l 0 -t {} {} -o {}".format(
                nt,samfile_path_m, bam_m
            )
        ) # mouse bam

    # Disambiguate step
    os.system(
        "ngs_disambiguate -a {} -s pdx -o {} {} {}".format(
            aligner, int_path, bam_h, bam_m)
    )
    disambiguated_hbam = os.path.join(int_path, 'pdx.disambiguatedSpeciesA.bam')
    disambiguated_mbam = os.path.join(int_path, 'pdx.disambiguatedSpeciesB.bam')

    # Convert to fastq
    # human
    os.system(
        "perl {} {} {} {}".format(
            bam2fastq, disambiguated_hbam, output, nt)
    )
    for i in range(1,3):
        fastq_name = os.path.join(output,'fastq{}.fastq'.format(i))
        os.rename(fastq_name, fastq_name + '.human')
    # mouse
    os.system(
        "perl {} {} {} {}".format(
            bam2fastq, disambiguated_mbam, output, nt)
    )
    for i in range(1,3):
        fastq_name = os.path.join(output,'fastq{}.fastq'.format(i))
        os.rename(fastq_name, fastq_name + '.mouse')
    
    # deleting intermadiate files
    if not keep:
        os.system('rm -f -r {}/'.format(int_path))
