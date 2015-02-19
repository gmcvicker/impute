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
# jobs need 2 GB of mem:
#$ -l h_vmem=2g


CHROM_ARRAY=({1..22}) 
i=`expr $SGE_TASK_ID-1`
CHROM=`echo ${CHROM_ARRAY[$i]}`

POP=CEU

echo "chrom$CHROM, host:$HOSTNAME" >&2

SCRIPT=$HOME/proj/script/perl/IMPUTE/convert_known_haps_to_vcf_format.pl

DIR=/KG/gmcvicker/IMPUTE/

perl $SCRIPT \
    -hap $DIR/output/$POP/hg18/chr$CHROM.hg18.impute2_haps.gz \
    -vcf $DIR/output/$POP/hg18/vcf/chr$CHROM.hg18.impute2_haps.vcf.gz \
    -allele_p $DIR/output/$POP/hg18/chr$CHROM.hg18.impute2_allele_probs.gz \
    -chr $CHROM \
    -samp_tab $DIR/samples/hapmap_phase2_YRI_parents.sample_table



