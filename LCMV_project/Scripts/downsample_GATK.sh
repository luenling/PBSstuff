#!/bin/bash
#----------
# author: Lukas Endler
# date: 20.9.2015 at 12:23
# takes a bam file and downsamples it to a certain max coverage (default 1000)
# removes unmapped reads, improper pairs (on different chromosomes)
# bash command.sh {bamfile} {max coverage}
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses.fasta
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
GATK=/Volumes/Temp/Lukas/LCMV_project/Tools/GenomeAnalysisTK-3.4-46.jar
SAMTOOLS=/usr/local/bin/samtools
BWA=/usr/local/Cellar/bwa/0.7.12/bin/bwa
FN=`basename $1 .bam`
MAXCOV=10000
if [ $2 ] ; then
    MAXCOV=$2
fi
LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log

RFS="-rf BadMate -rf MappingQuality -mmq 30 -rf MateSameStrand -rf NotPrimaryAlignment -rf UnmappedRead"

echo start GATK for downsampling and filtering at `date` >> $LOGFILE
echo java -Xmx4g -jar $GATK -T PrintReads -R $REFGENOME -S lenient -I $1 -o ${FN}_ds${MAXCOV}.bam \
     -dcov $MAXCOV $RFS >> $LOGFILE
java -Xmx4g -jar $GATK -T PrintReads -R $REFGENOME -S lenient -I $1 -o ${FN}_ds${MAXCOV}.bam \
     -dcov $MAXCOV $RFS >> $LOGFILE 2>> $ERRORLOG
ES=$?
echo finished downsampling and filtering at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
# echo start depthofcoverage for file  ${FN}_ds${MAXCOV}.bam  at `date` >> $LOGFILE
# if [ ! -d ${FN} ] ; then mkdir $FN; fi 
# echo java -Xmx4g -jar $GATK -T DepthOfCoverage -R $REFGENOME -S lenient -I ${FN}_ds${MAXCOV}.bam -o ${FN}/${FN}_ds${MAXCOV} >> $LOGFILE
# java -Xmx4g -jar $GATK -T DepthOfCoverage -R $REFGENOME -S lenient -I $1 -o ${FN}/${FN}_ds${MAXCOV} >> $LOGFILE 2>> $ERRORLOG &
# ES=$?
# echo finished depth of coverage at `date` with exit state $ES >> $LOGFILE
# [ $ES -eq 0 ] || exit $ES




