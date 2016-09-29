#!/bin/bash
#-----------------------------------------------------------------------------------------------------------------------------------------------

BN=`basename $1 _1.fq.gz`
BN=${BN/\#/_smp}
echo at `date` started clipping of $BN >> $BN.log
echo trimmomatic PE -threads 10 -phred33 -basein $1 -baseout ${BN}.fq.gz HEADCROP:20 >> ${BN}.log
trimmomatic PE -threads 10 -phred33 -basein $1 -baseout ${BN}.fq.gz HEADCROP:20 >> ${BN}.log 2>&1
E=$?
echo finished at `date` with exit state $E >> ${BN}.log


