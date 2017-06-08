#!/bin/bash
#----------
# author: Lukas Endler
# date: 20.9.2015 at 16:46
# takes a bam file calls variants with varscan2
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses_short.fasta
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
GATK=/Volumes/Temp/Lukas/LCMV_project/Tools/GenomeAnalysisTK-3.4-46.jar
SAMTOOLS=/usr/local/bin/samtools
VARSCAN=$BASEDIR/Tools/varscan/VarScan.v2.4.0.jar
if [[ $1 =~ "mpileup" ]] ; then
    FN=`basename $1 .mpileup`
else
   FN=`basename $1 .list`
fi

LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log

echo starting mpileup and varscan at `date` >> $LOGFILE

if  [[ $1 =~ "mpileup" ]] ; then
    echo cat $1 \| java -Xmx8g -jar $VARSCAN mpileup2snp --min-coverage 50 --strand-filter 1 \
	 --output-vcf 1 --min-var-freq 0.001 --vcf-sample-list $2 \>  ${FN}_b_varscan.vcf >> $LOGFILE
    cat $1 | java -Xmx8g -jar $VARSCAN mpileup2snp --min-coverage 50 --strand-filter 1 \
	 --output-vcf 1 --min-var-freq 0.001 --vcf-sample-list $2 >  ${FN}_b_varscan.vcf 2>> $ERRORLOG
else
    echo samtools mpileup -q 30 -Q 30 -d 10000000 -f $REFGENOME -b $1 \| \
    	 java -Xmx8g -jar $VARSCAN mpileup2snp --min-coverage 50 --strand-filter 1 \
    	 --output-vcf 1  --min-var-freq 0.001 --vcf-sample-list  $2 \>  ${FN}_varscan.vcf >> $LOGFILE
    samtools mpileup -q 30 -Q 30 -d 10000000 -f $REFGENOME -b $1 | \
    	java -Xmx8g -jar $VARSCAN mpileup2snp --min-coverage 50 --strand-filter 1 \
    	     --output-vcf 1 --min-var-freq 0.001 --vcf-sample-list $2 >  ${FN}_varscan.vcf 2>> $ERRORLOG
fi
ES=$?
echo finished varscan at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES


exit 0

java -Xmx10g -jar  ~/LCMV_project/Tools/popoolation2_1201/mpileup2sync.jar --input *.mpileup --output q30_Q30.sync
python  ~/LCMV_project//Scripts/get_allele_freqs.py -i q30_Q30.sync --all 
for i in 0.01 0.05 0.1; do
python  ~/LCMV_project/Scripts/add_afs_filter_varscan.py -i  ${FN}_varscan.vcf -d 50 -m $i -a *.af > ${FN}_varscan_afs_${i}.vcf
bash  ~/LCMV_project/Scripts/call_snpeff.sh ${FN}_varscan_afs_${i}.vcf
done



# #samtools mpileup -B -q 30 -Q 30 -d 10000000 -f $REFGENOME -b $1 | java -Xmx8g -jar $VARSCAN mpileup2snp --min-coverage 50 --strand-filter 1 --output-vcf 1  --vcf-sample-list $2 >  ${FN}_varscan.vcf
# echo starting mpileup on file at `date` >> $LOGFILE
# echo samtools mpileup -B -q 30 -Q 30 -d 10000000 -f $REFGENOME -b $1 >> $LOGFILE
# samtools mpileup -B -q 30 -Q 30 -d 10000000 -f $REFGENOME -b $1 > ${FN}.mpileup 2>> $ERRORLOG &
# # ES=$?
# # echo finished mpileup at `date` with exit state $ES >> $LOGFILE
# # [ $ES -eq 0 ] || exit $ES
# echo cat ${FN}.mpileup \| java -Xmx8g -jar $VARSCAN mpileup2snp --min-coverage 50 --strand-filter 1 --output-vcf 1 --vcf-sample-list $2 \> ${FN}_varscan.vcf >> $LOGFILE
# cat ${FN}.mpileup | java -Xmx8g -jar $VARSCAN mpileup2snp --min-coverage 50 --strand-filter 1 --output-vcf 1  --vcf-sample-list $2 >  ${FN}_varscan.vcf
# ES=$?
# echo finished mpileup at `date` with exit state $ES >> $LOGFILE
# [ $ES -eq 0 ] || exit $ES
