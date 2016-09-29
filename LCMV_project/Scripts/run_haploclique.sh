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
HAPLOCLIQUE=$BASEDIR/Tools/haploclique/bin/haploclique-assembly
export SAF=$BASEDIR/Tools/

FN=`basename $1 .bam`
# just to replace the whole bloody name tag
if [ ! -d $FN ]; then
    mkdir $FN
fi
cd $FN
LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log
echo "start Haploclique for " $1 " at" `date` >> $LOGFILE
echo $HAPLOCLIQUE -r $REFGENOME -i $1 >> $LOGFILE
 $HAPLOCLIQUE -r $REFGENOME -i $1 2>> $ERRORLOG >> $LOGFILE
ES=$?
echo finished haploclique at `date` with exit state $ES >> $LOGFILE


