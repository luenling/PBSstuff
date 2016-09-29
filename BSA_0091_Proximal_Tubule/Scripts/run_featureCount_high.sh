# 500 nt highest
BAMFILES=~/BSA_0091_Proximal_Tubule/mm10/featCounts/all_RN_sorted.bam
GTF=~/BSA_0091_Proximal_Tubule/mm10/mm10_genome/ucsc_mm10_ensembl_exons_500nt_most_highly.gtf
OUTFILE=~/BSA_0091_Proximal_Tubule/mm10/featCounts/mm10.pseudo.500nt.high.featurecount.counts
LOGFILE=${OUTFILE}.log
featureCounts -T 14 -pB -Q 10 -g gene_id -s 0 -F GTF -a $GTF -o $OUTFILE $(cat $BAMFILES) 2> $LOGFILE 

EX=$?
if [ $EX -ne 0 ]
then
        echo "creating $OUTFILE died at `date` with exit status $EX"  >>$LOGFILE
else
        echo $OUTFILE "finished on `date`" >> $LOGFILE
        grep Success $LOGFILE | awk -v FS="[(%]" 'BEGIN{max=0; min=100000} {sum+=$2; count+=1; if ($2>max) max=$2; if ($2 < min) min = $2; } END{print sum/count, max, min}' >> $LOGFILE
fi


# 1000 nt highest
BAMFILES=~/BSA_0091_Proximal_Tubule/mm10/featCounts/all_RN_sorted.bam
GTF=~/BSA_0091_Proximal_Tubule/mm10/mm10_genome/ucsc_mm10_ensembl_exons_1000nt_most_highly.gtf
OUTFILE=~/BSA_0091_Proximal_Tubule/mm10/featCounts/mm10.pseudo.1000nt.high.featurecount.counts
LOGFILE=${OUTFILE}.log
featureCounts -T 14 -pB -Q 10 -g gene_id -s 0 -F GTF -a $GTF -o $OUTFILE $(cat $BAMFILES) 2> $LOGFILE 

EX=$?
if [ $EX -ne 0 ]
then
        echo "creating $OUTFILE died at `date` with exit status $EX"  >>$LOGFILE
else
        echo $OUTFILE "finished on `date`" >> $LOGFILE
        grep Success $LOGFILE | awk -v FS="[(%]" 'BEGIN{max=0; min=100000} {sum+=$2; count+=1; if ($2>max) max=$2; if ($2 < min) min = $2; } END{print sum/count, max, min}' >> $LOGFILE
fi
