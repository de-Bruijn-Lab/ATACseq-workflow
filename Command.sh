#!/bin/sh

#$ -cwd
#$ -q batchq
#$ -M username
#$ -m eas

module load homer
module load python
module load bedtools
python /t1-data/user/jharman/Scripts/ATACseq/ATAC_Pipeline.py -genome mm9 -s ATAC_Intersect_Samples.txt -sDB ATAC_DiffBind_Samples.csv 
