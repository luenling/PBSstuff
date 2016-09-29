#!/bin/bash
#----------
# author: Lukas Endler
# date: 17.9.2015 at 16:34
# takes a bam file and uses quasirecomb 
# call with bash command blub_1.fq blub_2.fq >> logfile.log 2>> log.error.log
# also sort the file afterwards
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses.fasta
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
QUASIRECOMB=$BASEDIR/Tools/QuasiRecomb.jar

FN=`basename $1 .bam`
# just to replace the whole bloody name tag
if [ ! -d $FN ]; then
    mkdir $FN
fi
LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log
echo "start Quasirecomb for " $1 " at" `date` >> $LOGFILE
echo java -XX:+UseParallelGC -XX:+UseNUMA -XX:NewRatio=9 -Xms2G -Xmx10G -jar $QUASIRECOMB -o $FN -noRecomb -quality -maxDel 2 -maxPercDel 0.1  -i $1 >> $LOGFILE
java -XX:+UseParallelGC -XX:+UseNUMA -XX:NewRatio=9 -Xms2G -Xmx10G -jar $QUASIRECOMB -o $FN -noRecomb -quality -maxDel 2 -maxPercDel 0.1  -i $1 2>> $ERRORLOG >> $LOGFILE
ES=$?
echo finished quasirecomb at `date` with exit state $ES >> $LOGFILE
