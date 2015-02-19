#!/bin/sh
#
#
# preserve shell environment vars:
#$ -V
#
# run a job-array consisting of following range of task ids
#$ -t 1-22
#
# write output and err files to following dirs
#$ -o /mnt/lustre/home/gmcvicker/sge/out
#$ -e /mnt/lustre/home/gmcvicker/sge/err
#
# jobs need 12GB of mem:
#$ -l h_vmem=12g

# get chromosome for this SGE task id
#CHROM_ARRAY=({1..22} X Y) 
CHROM_ARRAY=({1..22}) 
i=`expr $SGE_TASK_ID-1`
CHROM=`echo ${CHROM_ARRAY[$i]}`

# POP=YRI
POP=CEU

echo "chrom$CHROM, host:$HOSTNAME" >&2

OUTPUT_PATH=/KG/gmcvicker/IMPUTE/hapmap_haplotypes/phaseII/r22/hg18/haplotypes_chr$CHROM.$POP.hg18.txt.gz

python $HOME/proj/impute/python/make_phased_hapmap_genotype_file.py $POP chr$CHROM | gzip > $OUTPUT_PATH
