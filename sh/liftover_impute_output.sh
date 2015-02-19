#!/bin/sh
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
# jobs need 4GB of mem:
#$ -l h_vmem=4g


#CHROM_ARRAY=({1..22} X Y) 
CHROM_ARRAY=({1..22}) 

# get chromosome for this SGE task id
i=`expr $SGE_TASK_ID-1`
CHROM=`echo ${CHROM_ARRAY[$i]}`

echo "chr$CHROM, pop:$POP, hostname:$HOSTNAME" >&2

###
# CHANGE THIS to the pop used (YRI/CEU)
#
POP=YRI

###
# CHANGE THIS to change the type of file that is
# lifted over:
#      impute2, impute2_haps, impute2_allele_probs.gz
# are lifted over
# TYPE=impute2
# TYPE=impute2_haps
TYPE=impute2_allele_probs


INPUT_DIR=/KG/gmcvicker/IMPUTE/output/$POP/hg19
INPUT_FILE=chr$CHROM.hg19.$TYPE.gz

OUT_DIR=/KG/gmcvicker/IMPUTE/output/$POP/hg18

mkdir -p $OUT_DIR

BED_FILE_IN=`echo $INPUT_FILE | sed s/.gz/.hg18.bed/`
BED_FILE_OUT=`echo $INPUT_FILE | sed s/.gz/.hg19.bed/`
BED_FILE_UNMAPPED=`echo $INPUT_FILE | sed s/.gz/.unmapped.bed/`
OUT_FILE=$OUT_DIR/chr${CHROM}.hg18.$TYPE.gz

TEMP=/KG/gmcvicker/tmp

echo "making BED file" >&2
# the columns of the bed file that we use are:
# 1) CHROM
# 2) START (0-based)
# 3) END
# 4) NAME (put refsnp name in here)
# 5) SCORE (put 0.0 here)
# 6) STRAND (put + here)
# 
# The remaining columns are just the lines from the input file
# (the data we actually want to be lifted over)
# We include the last three so that strand can be indicated
# is mapped to - of hg18 assembly, because all SNPs need to be
# reported wrt + strand.
gunzip -c $INPUT_DIR/$INPUT_FILE |  awk 'BEGIN{ORS=""} {print "chr"'$CHROM',$3-1,$3,$1,"0.0","+ "; gsub(" ", ",", $0); print $0"\n"}' > $TEMP/$BED_FILE_IN

if [ "$?" -ne "0" ]; then
    echo "making BED file failed" 1>&2
    exit 1
fi

echo "lifting over BED file"  >&2
# map intermediate BED file to hg18 with liftover
/home/rpique/bin/x86_64/liftOver -bedPlus=6 $TEMP/$BED_FILE_IN $HOME/data/ucsc/hg19/liftover/hg19ToHg18.over.chain.gz  $TEMP/$BED_FILE_OUT $TEMP/$BED_FILE_UNMAPPED

if [ "$?" -ne "0" ]; then
    echo "liftover failed" 1>&2
    exit 1
fi


echo "parsing and filtering BED output"  >&2
# Sort and filter SNPs that mapped to different strand or chromosome
# on hg19. Also discard SNPs with ordering that is inconsistent between
# hg18 and hg19. Convert BED file back to haplotype format.
python $HOME/proj/script/python/10_IND/IMPUTE/filter_bed_intermediate.py chr$CHROM $TEMP/$BED_FILE_OUT | gzip > $OUT_FILE

if [ "$?" -ne "0" ]; then
    echo "filtering of BED file failed" 1>&2
    exit 1
fi


# cleanup temporary files
rm $TEMP/$BED_FILE_IN
rm $TEMP/$BED_FILE_UNMAPPED
rm $TEMP/$BED_FILE_OUT


echo "done" >&2
