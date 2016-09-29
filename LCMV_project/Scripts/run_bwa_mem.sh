#!/bin/bash
#----------
# author: Lukas Endler
# date: 17.9.2015 at 16:34
# takes two fastq files runs bwa with 12 threads and outputs a bam file called like the fastq prefix
# call with bash command blub_1.fq blub_2.fq >> logfile.log 2>> log.error.log
# also sort the file afterwards
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses.fasta
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
BWA=/usr/local/Cellar/bwa/0.7.12/bin/bwa
FN=`basename $1 .gz`
FN=`basename $FN _1.fq`
FN=`basename $FN _1.fastq`
FN=`basename $FN _1P.fq`
# just to replace the whole bloody name tag
FN=${FN/_H325VBBXX_5#/_smp_}
# or just the neary as bloody part
FN=${FN/_H325VBBXX_5/}
SAMPLE=$FN
LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log

#RG='@RG\tID:'$FN'\tSM:'$FN
echo "start bwa mem  at" `date` >> $LOGFILE
echo bwa mem -R "@RG\tID:${FN}\tSM:${FN}" -k 15 -r 1.25 -M -t 10 $REFGENOME $1 $2 \| samtools view -Shb - \> $FN".bam"  >> $LOGFILE
bwa mem -R "@RG\tID:${FN}\tSM:${FN}"  -k 15 -r 1.25 -M -t 10 $REFGENOME $1 $2 | samtools view -Shb - > $FN".bam" 2>> $ERRORLOG
ES=$?
echo finished bwa mem mapping at `date` with exit state $ES >> $LOGFILE
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


