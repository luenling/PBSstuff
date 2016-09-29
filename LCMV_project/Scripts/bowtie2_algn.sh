#!/bin/bash
#----------
# author: Lukas Endler
# date: 06.10.2015 at 15:346
# takes 2 fastq files and creates a bowtie2 aligment 
#--------------

BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses.bt2
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
FN=`basename $1 .gz`
FN=`basename $FN _1.fq`
FN=`basename $FN _1.fastq`
# just to replace the whole bloody name tag
FN=${FN#*H325VBBXX_5_}
FN=${FN%_trimm*}

LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log

echo  start mapping of $1 at `date` >> $LOGFILE
echo bowtie2 -x $REFGENOME --rg-id $FN --rg "SM:$FN" -p 10 --met-stderr -N 0 -L 16 -1 $1 -2 $2 >> $LOGFILE
bowtie2 -x $REFGENOME --rg-id $FN --rg "SM:$FN" -p 10 --met-stderr -N 0 -L 16 -1 $1 -2 $2\
	2>> $ERRORLOG | samtools view -Sbh - > ${FN}.bam
ES=$?
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
echo start sorting at `date`  >> $LOGFILE
echo java  -Xmx6g -jar $PICARD SortSam I=${FN}.bam O=${FN}_sorted.bam SO=coordinate  >> $LOGFILE
java  -Xmx6g -jar $PICARD SortSam I=${FN}.bam O=${FN}_sorted.bam SO=coordinate  >> $LOGFILE 2>> $ERRORLOG
ES=$?
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
samtools index ${FN}_sorted.bam
echo flagstat >> $LOGFILE
samtools flagstat ${FN}_sorted.bam >> $LOGFILE
echo idxstats >> $LOGFILE
samtools idxstats ${FN}_sorted.bam >> $LOGFILE

