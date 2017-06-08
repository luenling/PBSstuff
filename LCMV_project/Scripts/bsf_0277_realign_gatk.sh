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
GATK=/Volumes/Temp/Lukas/LCMV_project/Tools/GenomeAnalysisTK-3.4-46.jar



FN=`basename $1 .bam`

LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log

echo "start Realtargetcreator at" `date` >> $LOGFILE

echo java -Xmx4g -jar $GATK -T RealignerTargetCreator -dt none -R $REFGENOME -I $1 -o ${FN}.intervals   >> $LOGFILE
java -Xmx4g -jar $GATK -T RealignerTargetCreator -dt none -R $REFGENOME -I $1 -o ${FN}.intervals >> $LOGFILE 2>> $ERRORLOG

ES=$?
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES

echo "start Realigner  at" `date` >> $LOGFILE

echo java -Xmx4g  -jar $GATK -T IndelRealigner -dt none  -R $REFGENOME -I $1 -model USE_SW -targetIntervals ${FN}.intervals -o ${FN}_real.bam >> $LOGFILE



java -Xmx4g  -jar $GATK -T IndelRealigner -dt none  -R $REFGENOME -I $1 -model USE_SW -targetIntervals ${FN}.intervals -o ${FN}_real.bam >> $LOGFILE 2>> $ERRORLOG
ES=$?
echo finished  at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
																										       
