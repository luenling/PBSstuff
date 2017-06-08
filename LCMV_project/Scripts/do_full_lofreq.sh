#!/bin/bash
#----------
# author: Lukas Endler
# Time-stamp: <2017-01-23 12:54:05 lukasendler>
# date: 20.9.2015 at 12:23
# takes a file with list of bam files and vcf from samtools for depths and calls variants with lofreq2
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
SCRIPTS=${BASEDIR}/Scripts
REFGENOME=$BASEDIR/References/viruses_short.fasta
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
GATK=/Volumes/Temp/Lukas/LCMV_project/Tools/GenomeAnalysisTK-3.4-46.jar
SAMTOOLS=/usr/local/bin/samtools
LOFREQ=$BASEDIR/Tools/lofreq_star-2.1.2/bin/lofreq
FN=`basename $1 .list`
FN=${FN/_noprime*/}
FN=${FN/_sorted*/}
LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log

# do call lofreq
while read file; do
    bash $SCRIPTS/call_lofreq.sh $file
done < $1

bash $SCRIPTS/combine_lofreq_vars.sh

bash $SCRIPTS/call_lofreq_withbed.sh ./all_samp_bcf.bed $2



