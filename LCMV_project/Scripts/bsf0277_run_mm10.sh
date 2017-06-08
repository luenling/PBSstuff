#!/bin/bash
#----------
# author: Lukas Endler
# date: 15.12.2015 at 11:07
# Time-stamp: <2017-01-25 14:32:34 lukasendler>
# takes two fastq files runs bwa with 12 threads and outputs a bam file called like the fastq prefix
# call with bash command blub_1.fq blub_2.fq >> logfile.log 2>> log.error.log
# also sort the file afterwards
#--------------


BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses_Mus_musculus.GRCm38.fa.gz
PICARD=/usr/local/Cellar/picard-tools/2.5.0/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
BWA=/usr/local//Cellar/bwa/0.7.15/bin/bwa

FN=`basename $1 .gz`
FN=`basename $FN _1.fq`
FN=`basename $FN _1.fastq`
# just to replace the whole bloody name tag
R2=`echo $1 | sed 's/_1\.f/_2.f/'`
FN=${FN/*_HFNMYBBXX_5\#/}
SAMPLE=${FN/_S19[1-9]/}
LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log

RG='@RG\tID:'$FN'\tLB:'$SAMPLE'\tSM:'$SAMPLE'\tPL:illumina\tPU:unitx'
echo "start bwa mem  at" `date` >> $LOGFILE
echo $BWA mem -R $RG -k 17 -r 1.25 -M -t 17 $REFGENOME $1 $R2  2\>\> $ERRORLOG  \| $SAMTOOLS view -Shb - \| $SAMTOOLS sort -T ${FN}_temp - \> $FN"_sorted.bam"  >> $LOGFILE

$BWA mem -R $RG  -k 17 -r 1.25 -M -t 17 $REFGENOME $1 $R2 2>> $ERRORLOG | $SAMTOOLS view -Shb - |  $SAMTOOLS sort -T ${FN}_temp - > $FN"_sorted.bam"
ES=$?
echo finished bwa mem mapping at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
#echo start sorting at `date`  >> $LOGFILE
#echo java  -Xmx6g -jar $PICARD SortSam I=${FN}.bam O=${FN}_sorted.bam SO=coordinate  >> $LOGFILE
#java  -Xmx6g -jar $PICARD SortSam I=${FN}.bam O=${FN}_sorted.bam SO=coordinate  >> $LOGFILE 2>> $ERRORLOG
#ES=$?
#echo finished at `date` with exit state $ES >> $LOGFILE
#[ $ES -eq 0 ] || exit $ES
{
samtools index ${FN}_sorted.bam
echo flagstat >> $LOGFILE
samtools flagstat ${FN}_sorted.bam >> $LOGFILE
echo idxstats >> $LOGFILE
samtools idxstats ${FN}_sorted.bam >> $LOGFILE
# extract only the viral genomes
samtools view -bh -f 2 -F 256 ${FN}_sorted.bam 'gi|86440167|gb|DQ361066.1|' 'gi|116563461|gb|DQ361065.2|' > ${FN}_sorted_viral.bam
samtools index ${FN}_sorted_viral.bam
echo flagstat >> $LOGFILE
samtools flagstat ${FN}_sorted_viral.bam >> $LOGFILE
echo idxstats >> $LOGFILE
samtools idxstats ${FN}_sorted_viral.bam >> $LOGFILE
} &




