#!/bin/bash
#----------
# author: Lukas Endler
# created: 2016-12-19 10:34:11
# Time-stamp: <2017-01-25 14:39:03 lukasendler>
# takes a bam file, reheaders it to only the two viral segments in short format 
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses_short.fasta
PICARD=/usr/local/Cellar/picard-tools/2.5.0/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
BWA=/usr/local//Cellar/bwa/0.7.15/bin/bwa

FN=`basename $1 .bam`

LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log

echo reheader $FN at `date` >> $LOGFILE
echo $SAMTOOLS reheader \<\($SAMTOOLS view -H $1 \| grep -Ev \''^@SQ.*SN:[^g]'\' \| sed \''s/gi|86440167|gb|DQ361066.1|/L/; s/gi|116563461|gb|DQ361065.2|/S/'\'\) $1 \> ${FN}_rh.bam >> $LOGFILE
$SAMTOOLS reheader <($SAMTOOLS view -H $1 | grep -Ev '^@SQ.*SN:[^g]' | sed 's/gi|86440167|gb|DQ361066.1|/L/; s/gi|116563461|gb|DQ361065.2|/S/') $1 > ${FN}_rh.bam 2>> $ERRORLOG
ES=$?
echo finished reheader at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES

$SAMTOOLS index ${FN}_rh.bam
echo flagstat ${FN}_rh.bam >> $LOGFILE
$SAMTOOLS flagstat ${FN}_rh.bam >> $LOGFILE
echo idxstats ${FN}_rh.bam >> $LOGFILE
$SAMTOOLS idxstats ${FN}_rh.bam >> $LOGFILE
sam-stats ${FN}_rh.bam > ${FN}_rh.stats

#echo creating coverages at `date` >> $LOGFILE
#if [ ! -d Coverages ]; then mkdir Coverages ; fi




