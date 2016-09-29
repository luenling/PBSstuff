#!/bin/bash
#----------
# author: Lukas Endler
# date: 17.9.2015 at 16:34
# takes two fastq files runs bwa with 12 threads and outputs a bam file called like the fastq prefix
# call with bash command blub_1.fq blub_2.fq >> logfile.log 2>> log.error.log
# also sort the file afterwards
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOMEM=$BASEDIR/References/viruses_short.dat
REFGENOME=$BASEDIR/References/viruses_short.fasta
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
MOSAIKB=$BASEDIR/Tools/MOSAIK/bin/MosaikBuild
MOSAIKA=$BASEDIR/Tools/MOSAIK/bin/MosaikAligner

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


echo $MOSAIKB -q $1 -q2 $2 -out ${FN}_reads.mkb -sam ${FN} -id ${FN} -mfl 300 -ln ${FN} -st illumina >> $LOGFILE
$MOSAIKB -q $1 -q2 $2 -out ${FN}_reads.mkb -sam ${FN} -id ${FN} -mfl 300 -st illumina 2>> ${FN}.err.log
echo start aligning with mosaik at `date`  >> $LOGFILE
echo $MOSAIKA -ia $REFGENOMEM -gop 50 -hgop 40 -gep 15 -mmp 0.07 -minp 0.9 -p 10 -in ${FN}_reads.mkb -out ${FN}_mosaik -annpe $BASEDIR/Tools/MOSAIK/src/networkFile/2.1.78.pe.ann -annse $BASEDIR/Tools/MOSAIK/src/networkFile/2.1.78.se.ann >> $LOGFILE
$MOSAIKA -ia $REFGENOMEM -gop 50 -hgop 40 -gep 15 -mmp 0.07 -minp 0.9 -p 10 -in ${FN}_reads.mkb -out ${FN}_mosaik -annpe $BASEDIR/Tools/MOSAIK/src/networkFile/2.1.78.pe.ann -annse $BASEDIR/Tools/MOSAIK/src/networkFile/2.1.78.se.ann 2>> ${FN}.err.log
ES=$?
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES

ID=${FN/_nop*/}

echo picard AddOrReplaceReadGroups I=${FN}_mosaik.bam O=${FN}_mosaik_sorted.bam  VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate RGID=${ID}_noprime RGLB=$ID RGPL=illumina RGPU=unitx RGSM=$ID >> $LOGFILE
picard AddOrReplaceReadGroups I=${FN}_mosaik.bam O=${FN}_mosaik_sorted.bam  VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate RGID=${ID}_noprime RGLB=$ID RGPL=illumina RGPU=unitx RGSM=$ID >> $LOGFILE 2>> ${FN}.err.log
ES=$?
echo finished adding RG  at `date` with exit state $ES >> $LOGFILE
$SAMTOOLS index ${FN}_mosaik_sorted.bam
echo flagstat ${FN}_mosaik_sorted.bam >> $LOGFILE
samtools flagstat ${FN}_mosaik_sorted.bam >> $LOGFILE
echo idxstats ${FN}_mosaik_sorted.bam >> $LOGFILE
samtools idxstats ${FN}_mosaik_sorted.bam >> $LOGFILE

sam-stats ${FN}_mosaik_sorted.bam > ${FN}_mosaik_sorted.stats




