#!/bin/sh
#
#
# preserve shell environment vars:
#$ -V
#
# run a job-array consisting of following range of task ids
#$ -t 1-571
#
# write output and err files to following dirs
#$ -o /mnt/lustre/home/gmcvicker/sge/out
#$ -e /mnt/lustre/home/gmcvicker/sge/err
#
# jobs need 4GB of mem:
#$ -l h_vmem=4g,bigio=1

TASKFILE=/KG/gmcvicker/IMPUTE/tasklist.hg19.txt
# POP=YRI
POP=CEU

LINE=`head -$SGE_TASK_ID $TASKFILE | tail -1`

CHR=`echo $LINE | awk '{print $2}'`
START=`echo $LINE | awk '{print $3}'`
END=`echo $LINE | awk '{print $4}'`
 
echo "$CHR:$START-$END, hostname=$HOSTNAME" >&2

mkdir -p /KG/gmcvicker/IMPUTE/output/$POP/hg19/$CHR

# IMPUTE2=/data/share/impute_v2.2.2_x86_64_dynamic/impute2
IMPUTE2=/mnt/lustre/home/gmcvicker/impute2/impute_v2.2.2_x86_64_dynamic/impute2


GENETIC_MAP=/mnt/lustre/home/bhowie/genetic_maps_b37/genetic_map_${CHR}_combined_b37.txt

PREPHASE_HAPS=/KG/gmcvicker/IMPUTE/hapmap_haplotypes/phaseII/r22/hg19/haplotypes_${CHR}.$POP.hg19.txt.gz

REFERENCE_HAPS=/home/bhowie/reference_panels/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_${CHR}_impute.hap.gz

LEGEND=/home/bhowie/reference_panels/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_${CHR}_impute.legend.gz

OUT=/KG/gmcvicker/IMPUTE/output/$POP/hg19/$CHR/${CHR}.$START.$END.impute2

$IMPUTE2 \
    -use_prephased_g \
    -m $GENETIC_MAP \
    -h $REFERENCE_HAPS \
    -l $LEGEND \
    -int $START $END \
    -known_haps_g $PREPHASE_HAPS \
    -Ne 20000 \
    -filt_rules_l 'afr.maf<0.004' \
    -phase \
    -o_gz \
    -o $OUT


if [ "$?" -ne "0" ]; then
    echo "IMPUTE2 failed" 1>&2
    exit 1
fi

echo "done" >&2