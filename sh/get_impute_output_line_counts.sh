#!/bin/sh
#


TASKFILE=/KG/gmcvicker/IMPUTE/tasklist.txt

N_TASK=`wc -l $TASKFILE | awk '{print $1}'`

for i in $(seq 1 $N_TASK)
do
    echo $i >&2

    LINE=`head -$i $TASKFILE | tail -1`

    CHR=`echo $LINE | awk '{print $2}'`
    START=`echo $LINE | awk '{print $3}'`
    END=`echo $LINE | awk '{print $4}'`
    
    PREFIX=/KG/gmcvicker/IMPUTE/output/hg19/${CHR}/${CHR}.$START.$END.impute2

    IMPUTE_FILE=$PREFIX.gz
    HAPS_FILE=${PREFIX}_haps.gz
    ALLELE_PROBS_FILE=${PREFIX}_allele_probs.gz
    
    echo $IMPUTE_FILE >&2
    echo $HAPS_FILE >&2
    echo $ALLELE_PROBS_FILE >&2

    IMPUTE_LINES=`gunzip -c $IMPUTE_FILE | wc -l | awk '{print $1}'`
    HAPS_LINES=`gunzip -c $HAPS_FILE | wc -l | awk '{print $1}'`
    ALLELE_PROBS_LINES=`gunzip -c $ALLELE_PROBS_FILE | wc -l | awk '{print $1}'`

    echo "$IMPUTE_LINES $HAPS_LINES $ALLELE_PROBS_LINES" >&2
done

