### inferring clonal evolution from somatic mutation and CNV results
#
# Require: R, python (2.7.x, not 3), PyClone, The clonevol and packcircles R packages
#
# design: design file. The samples in it will be used to infer one clonal evolutionary tree/graph. So these samples should be from one patient
#         three columns separated by \t, sample id, path to CNV result folder, and path to the mutation result folder (no headers)
# output: output folder. This folder will be deleted and recreated each time
# build: reference genome build

###########  prepare  ###################

#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';

my ($design,$output,$build)=@ARGV;
my ($tumor_contents,$tsv,$path,@items);

$path=abs_path($0);
$path=~s/evolution\.pl//;

###########  run PyClone  ###########

system("rm -f -r ".$output);
system("mkdir ".$output);
system("Rscript ".$path."/somatic_script/prepare_pyclone.R ".$design." ".$output." ".$build);

open(FILE,$output."/pyclone/parameters.txt");
$tsv=<FILE>;
$tsv=~s/\n//;
$tumor_contents=<FILE>;
$tumor_contents=~s/\n//;
close(FILE);

system("cd ".$output."/pyclone;PyClone run_analysis_pipeline --in_files ".$tsv." --working_dir ".$output."/pyclone ".
  "--prior total_copy_number --seed 1 --num_iters 5000 --burnin 1000 --thin 5 --tumour_contents ".$tumor_contents);

###########  run clonevol ##########

system("Rscript ".$path."/somatic_script/clonevol.R ".$output);

###########  example  ############
#
# perl evolution.pl \
# /home2/twang6/software/cancer/somatic/example/evolution/design.txt \
# /home2/twang6/software/cancer/somatic/example/evolution/results \
# hg38
