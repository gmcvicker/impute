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
# jobs need 2GB of mem:
#$ -l h_vmem=2g


#CHROM_ARRAY=({1..22} X Y) 
CHROM_ARRAY=({1..22}) 

# get chromosome for this SGE task id
i=`expr $SGE_TASK_ID-1`
CHROM=`echo ${CHROM_ARRAY[$i]}`

# POP=YRI
POP=CEU

echo "chr$CHROM, hostname:$HOSTNAME" >&2

# PANEL_DIR=/data/share/hapGEplab/Data/Archives/SNP/HapMap/r28_7_Oct_2010_nb
# PANEL_DIR=/data/share/HapMap/genotypes/r28/2010-08_phaseII+III/forward
PANEL_DIR=/KG/gmcvicker/IMPUTE/hapmap_haplotypes/phaseII/r22/hg18

OUT_DIR=/KG/gmcvicker/IMPUTE/hapmap_haplotypes/phaseII/r22/hg19

HAP_FILE=haplotypes_chr${CHROM}.$POP.hg18.txt.gz

BED_FILE_IN=`echo $HAP_FILE | sed s/.gz/.hg18.bed/`
BED_FILE_OUT=`echo $HAP_FILE | sed s/.gz/.hg19.bed/`
BED_FILE_UNMAPPED=`echo $HAP_FILE | sed s/.gz/.unmapped.bed/`
OUT_FILE=$OUT_DIR/haplotypes_chr${CHROM}.$POP.hg19.txt.gz


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
gunzip -c $PANEL_DIR/$HAP_FILE |  awk 'BEGIN{ORS=""} {print "chr"'$CHROM',$3-1,$3,$1,"0.0","+ "; gsub(" ", ",", $0); print $0"\n"}' > $TEMP/$BED_FILE_IN

if [ "$?" -ne "0" ]; then
    echo "making BED file failed" 1>&2
    exit 1
fi

echo "lifting over BED file"  >&2
# map intermediate BED file to hg18 with liftover
/home/rpique/bin/x86_64/liftOver -bedPlus=6 $TEMP/$BED_FILE_IN $HOME/data/ucsc/hg18/liftover/hg18ToHg19.over.chain.gz  $TEMP/$BED_FILE_OUT $TEMP/$BED_FILE_UNMAPPED

if [ "$?" -ne "0" ]; then
    echo "liftover failed" 1>&2
    exit 1
fi


echo "parsing and filtering BED output"  >&2
# Sort and filter SNPs that mapped to different strand or chromosome
# on hg19. Also discard SNPs with ordering that is inconsistent between
# hg18 and hg19. Convert BED file back to haplotype format.
python $HOME/proj/impute/python/filter_bed_intermediate.py chr$CHROM $TEMP/$BED_FILE_OUT | gzip > $OUT_FILE

if [ "$?" -ne "0" ]; then
    echo "filtering of BED file failed" 1>&2
    exit 1
fi


# cleanup temporary files
rm $TEMP/$BED_FILE_IN
rm $TEMP/$BED_FILE_UNMAPPED
rm $TEMP/$BED_FILE_OUT


echo "done" >&2
