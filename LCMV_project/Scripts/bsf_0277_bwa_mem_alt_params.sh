#!/bin/bash
#----------
# author: Lukas Endler
# date: 17.9.2015 at 16:34
# takes a bam file and runs bwa mem with high gap penalty
# call with bash command bla.bam >> logfile.log 2>> log.error.log
# also sort the file afterwards
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses_short.fasta
PICARD=/usr/local/Cellar/picard-tools/2.5.0/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
BWA=/usr/local//Cellar/bwa/0.7.15/bin/bwa
FN=`basename $1 .gz`
FN=`basename $FN .bam`
FN=`basename $FN _1.fq`
FN=`basename $FN _1.fastq`
FN=`basename $FN _1P.fq`
# just to replace the whole bloody name tag
SAMPLE=${FN/_S19[1-9]*/}
SAMPLE=${Sample/_sorted*/}
FN=${FN/_sorted*/}
FN=${FN}_hgap
LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log

#RG='@RG\tID:'$FN'\tSM:'$SAMPLE
echo "start bwa mem  at" `date` >> $LOGFILE
#echo $BWA mem -R "@RG\tID:${FN}\tSM:${SAMPLE}" -B 20 -A 3 -O 30 -E 3 -k 15 -r 1.25 -M -t 10 $REFGENOME $1 $2 \| samtools view
echo picard SamToFastq Inter=True F=/dev/stdout I=$1 \| $BWA mem -R "@RG\tID:${FN}\tSM:${SAMPLE}" -B 4 -A 1 -O 10 -E 1 -k 15 -r 1.25 -M -t 10 -p $REFGENOME - \| samtools view -Shb -F 256 - \| $SAMTOOLS sort -T ${FN}_temp - \> $FN"_sorted.bam"   >> $LOGFILE
picard SamToFastq Inter=True F=/dev/stdout I=$1 | $BWA mem -R "@RG\tID:${FN}\tSM:${SAMPLE}" -B 4 -A 1 -O 10 -E 1 -k 15 -r 1.25 -M -t 10 -p $REFGENOME - | samtools view -Shb  -F 256 - | $SAMTOOLS sort -T ${FN}_temp - > $FN"_sorted.bam" 2>> $ERRORLOG
ES=$?
echo finished bwa mem mapping at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES

$SAMTOOLS index ${FN}_sorted.bam
echo flagstat ${FN}_sorted.bam >> $LOGFILE
samtools flagstat ${FN}_sorted.bam >> $LOGFILE
echo idxstats ${FN}_sorted.bam >> $LOGFILE
samtools idxstats ${FN}_sorted.bam >> $LOGFILE

sam-stats ${FN}_sorted.bam > ${FN}_sorted.stats





