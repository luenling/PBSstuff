#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
BASEDIR=/Volumes/Temp/Lukas/LCMV_project/
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
TRIMJAR=$BASEDIR/Scripts/TrimFastq.jar
#REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
INFILE=$1
OUT=`basename $INFILE .bam`
LOGFILE=$OUT.log
ERRORFILE=$OUT.err.log
TRIMLOG=$OUT.trim.log
OUTDIR=''
if [ $2 ]; then 
    OUTDIR=$2'/'
    if [ ! -d $2 ]; then mkdir $2 ; fi
fi  


echo converting bam to fastq at `date` >> ${LOGFILE}
echo java -Xmx4g -jar ${PICARD} SamToFastq VALIDATION_STRINGENCY=SILENT INPUT=$INFILE FASTQ=${OUTDIR}${OUT}_1.fastq SECOND_END_FASTQ=${OUTDIR}${OUT}_2.fastq \>\> ${LOGFILE} 2\> ${ERRORFILE} >> ${LOGFILE}
java -Xmx4g -jar ${PICARD} SamToFastq VALIDATION_STRINGENCY=SILENT INPUT=$INFILE FASTQ=${OUTDIR}${OUT}_1.fastq SECOND_END_FASTQ=${OUTDIR}${OUT}_2.fastq  >> ${LOGFILE} 2>> ${ERRORFILE}
ES=$?
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
echo gzipping fastq ${OUTDIR}${OUT}_1.fastq at `date` >> ${LOGFILE}
gzip -fq ${OUTDIR}${OUT}_1.fastq  >> ${LOGFILE} 2>> ${ERRORFILE}
ES=$?
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
echo gzipping fastq ${OUTDIR}${OUT}_2.fastq at `date` >> ${LOGFILE}
gzip -fq ${OUTDIR}${OUT}_2.fastq >> ${LOGFILE} 2>> ${ERRORFILE}
ES=$?
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
#echo starting trimming of the read files at `date` >> $LOGFILE
#echo java -jar $TRIMJAR --input1 ${OUTDIR}${OUT}_1.fastq.gz --input2 ${OUTDIR}${OUT}_2.fastq.gz --min-length 80 --multi-thread --no-5p-trim --output ${OUT}_trimmed --quality-threshold 20 \>\> $TRIMLOG 2\>\> $ERRORFILE >> $LOGFILE
#java -jar $TRIMJAR --input1 ${OUTDIR}${OUT}_1.fastq.gz --input2 ${OUTDIR}${OUT}_2.fastq.gz --min-length 80 --multi-thread --no-5p-trim --output ${OUT}_trimmed --quality-threshold 20 >> $TRIMLOG 2>> $ERRORFILE
#ES=$?
#echo finished at `date` with exit state $ES >> $LOGFILE
#[ $ES -eq 0 ] || exit $ES
#echo removing  ${OUTDIR}${OUT}_1.fastq.gz and  ${OUTDIR}${OUT}_1.fastq.gz >>  $TRIMLOG
#rm -f ${OUTDIR}${OUT}_[12].fastq.gz




