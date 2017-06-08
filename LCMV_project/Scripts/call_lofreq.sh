#!/bin/bash
#----------
# author: Lukas Endler
# Time-stamp: <2017-01-23 11:54:48 lukasendler>
# date: 20.9.2015 at 12:23
# takes a bam file calls variants with lofreq2
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses_short.fasta
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
GATK=/Volumes/Temp/Lukas/LCMV_project/Tools/GenomeAnalysisTK-3.4-46.jar
SAMTOOLS=/usr/local/bin/samtools
LOFREQ=$BASEDIR/Tools/lofreq_star-2.1.2/bin/lofreq
FN=`basename $1 .bam`
FN=${FN/_noprime*/}
FN=${FN/_sorted*/}
LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log
INF=$1
if [[ $1 != *IDSQ.bam ]]; then
    echo at `date`  >> $LOGFILE
    echo  $LOFREQ viterbi -f $REFGENOME $1 \| samtools sort -T ${FN}_temp - \> ${FN}_real_viterbi.bam >> $LOGFILE
    $LOFREQ viterbi -f $REFGENOME $1 | samtools sort -T ${FN}_temp - > ${FN}_real_viterbi.bam
    echo $LOFREQ indelqual --dindel -f $REFGENOME -o ${FN}_real_viterbi_IDQS.bam ${FN}_real_viterbi.bam  >> $LOGFILE
    $LOFREQ indelqual --dindel -f $REFGENOME -o ${FN}_real_viterbi_IDQS.bam ${FN}_real_viterbi.bam
    rm -f ${FN}_real_viterbi.bam
    samtools index ${FN}_real_viterbi_IDQS.bam
    INF=${FN}_real_viterbi_IDQS.bam
fi

echo "start lofreq2 at" `date`  >> $LOGFILE 
echo $LOFREQ call -f $REFGENOME --verbose -o ${FN}_lofreq.vcf -q 20 -Q 20 -m 20 -C 75 -a 0.05 --call-indels $INF >> $LOGFILE
$LOFREQ call -f $REFGENOME --verbose -o ${FN}_lofreq.vcf -q 20 -Q 20 -m 20 -C 75 -a 0.05 --call-indels $INF 2>> $ERRORLOG >> $LOGFILE
ES=$?
echo finished lofreq at `date` with exit state $ES >> $LOGFILE

