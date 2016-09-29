#!/bin/bash
#----------
# author: Lukas Endler
# date: 17.9.2015 at 16:34
# takes two fastq files runs bwa with 12 threads and outputs a bam file called like the fastq prefix
# call with bash command blub_1.fq blub_2.fq >> logfile.log 2>> log.error.log
# also sort the file afterwards
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses_short.fasta
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
VPHASER2=$BASEDIR/Tools/VPhaser-2-02112013/bin/variant_caller
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$BASEDIR/Tools/VPhaser-2-02112013/bamtools-2.3.0/lib
FN=`basename $1 .bam`
# just to replace the whole bloody name tag
#FN=${FN/*_smp/smp}
FN=${FN/_noprime*/}
#if [ ! -d $FN ]; then
#    mkdir $FN
#fi
#cd $FN
LOGFILE=${FN}.vphaser2.log
ERRORLOG=${FN}.vphaser2.err.log
OMP_NUM_THREADS=8

echo "start VPHASER2 for " $1 " at" `date` >> $LOGFILE
if [ ! -d VPHASER2/$FN ]; then
    echo mkdir >> $LOGFILE
    mkdir -p VPHASER2/$FN
fi
echo $VPHASER2 -i $1 -o VPHASER2/$FN >> $LOGFILE
$VPHASER2 -i $1 -o VPHASER2/$FN 2>> $ERRORLOG >> $LOGFILE
ES=$?
echo finished vphaser2 at `date` with exit state $ES >> $LOGFILE


