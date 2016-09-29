#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob

REFLAT=/Volumes/Temp/Lukas/BSA_0091_Proximal_Tubule/check_qual/mm10_refFlat.txt
PICARD=/usr/local/share/java/picard.jar
BASE=`basename $1 .bam`
BASE=${BASE#rnaseq_tophat_}
BASE=${BASE%_accep*}
echo java -Xmx1g -jar $PICARD CollectRnaSeqMetrics STRAND=NONE REF_FLAT=$REFLAT I=$1 CHART=${BASE}_chart.pdf O=${BASE}_out.txt 2\>\&1 \>\> ${BASE}.log
java -Xmx1g -jar $PICARD CollectRnaSeqMetrics STRAND=NONE REF_FLAT=$REFLAT I=$1 CHART=${BASE}_chart.pdf O=${BASE}_out.txt &> ${BASE}.log

