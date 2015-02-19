

RS_ID=rs9616970
# RS_ID=rs12260013

FILES=`ls /home/bhowie/reference_panels/ALL_1000G_phase1integrated_v3_impute/*.legend.gz`

for FILE in $FILES
do
	echo $FILE >&2
	gunzip -c $FILE | grep "^$RS_ID"
done


FILES=`ls /data/share/10_IND/IMPUTE/hg18/*.impute2.gz`

for FILE in $FILES
do
	echo $FILE >&2
	gunzip -c $FILE | grep "$RS_ID"
done
