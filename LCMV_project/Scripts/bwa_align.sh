#!/bin/bash
#----------
# author: Lukas Endler
# date: 06.10.2015 at 15:346
# takes 2 fastq files and creates a bowtie2 aligment 
#--------------

BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFBWA=$BASEDIR/References/viruses_bw6.2.fasta
REFGENOME=$BASEDIR/References/viruses.fasta
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
BWA=${BASEDIR}/Tools/bwa-0.6.2/bwa

FN=`basename $1 .gz`
FN=`basename $FN _1.fq`
FN=`basename $FN _1.fastq`
# just to replace the whole bloody name tag
FN=${FN#*H325VBBXX_5_}
FN=${FN%_trimm*}

LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log

SAMPLE=$FN
LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log

RG='@RG\tID:'$FN'\tSM:'$FN
echo "start bwa aln  at" `date` >> $LOGFILE
echo $BWA aln -n 10 -o 1 -e 15 -l 500 -k 1 -t 10 $REFBWA $1 \> ${FN}_1.sai >> $LOGFILE
echo $BWA aln -n 10 -o 1 -e 15 -l 500 -k 1 -t 10 $REFBWA $2 \> ${FN}_2.sai >> $LOGFILE

 $BWA aln  -n 10 -o 1 -e 15 -l 500 -k 1 -t 10 $REFBWA $1 > ${FN}_1.sai 2>> $ERRORLOG &
 $BWA aln  -n 10 -o 1 -e 15 -l 500 -k 1 -t 10 $REFBWA $1 > ${FN}_2.sai 2>> $ERRORLOG
wait
ES=$?
echo finished bwa aln mapping at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
echo "start bwa sampe  at" `date` >> $LOGFILE
echo  $BWA sampe -r $RG $REFBWA ${FN}_1.sai ${FN}_2.sai $1 $2 \|\
     samtools view -Sbh - \> ${FN}.bam >> $LOGFILE
 $BWA sampe -r $RG $REFBWA ${FN}_1.sai ${FN}_2.sai $1 $2 2>> $ERRORLOG |\
    samtools view -Sbh - > ${FN}.bam
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
