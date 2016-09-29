#!/bin/bash
#----------
# author: Lukas Endler
# date: 17.9.2015 at 16:34
# takes two fastq files runs bwa with 12 threads and outputs a bam file called like the fastq prefix
# call with bash command blub_1.fq blub_2.fq >> logfile.log 2>> log.error.log
# also sort the file afterwards
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses_Mus_musculus.GRCm38.fa.gz
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
BWA=/usr/local/Cellar/bwa/0.7.12/bin/bwa
FN=`basename $1 .gz`
FN=`basename $FN _1.fq`
FN=`basename $FN _1.fastq`
# just to replace the whole bloody name tag
R2=`echo $1 | sed 's/_1\.f/_2.f/'`
RM=`echo $1 | sed 's/_1\.f/.f/'`
FN=${FN/*_H325VBBXX_5_/}
SAMPLE=$FN
LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log

RG='@RG\tID:'$FN'_PE\tSM:'$FN
echo "start PE bwa mem  at" `date` >> $LOGFILE
echo bwa mem -R $RG -k 17 -r 1.25 -M -t 17 $REFGENOME $1 $R2 \| samtools view -Shb - \> $FN"_PE.bam"  >> $LOGFILE
bwa mem -R $RG  -k 17 -r 1.25 -M -t 17 $REFGENOME $1 $R2 2>> $ERRORLOG | samtools view -Shb - > $FN"_PE.bam"
ES=$?
echo finished bwa mem mapping at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
echo start sorting at `date`  >> $LOGFILE
echo java  -Xmx6g -jar $PICARD SortSam  VALIDATION_STRINGENCY=SILENT I=${FN}_PE.bam O=${FN}_sorted_PE.bam SO=coordinate  >> $LOGFILE
java  -Xmx6g -jar $PICARD SortSam  VALIDATION_STRINGENCY=SILENT I=${FN}_PE.bam O=${FN}_sorted_PE.bam SO=coordinate >> $LOGFILE 2>> $ERRORLOG
ES=$?
rm -f  $FN"_PE.bam"
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
samtools index ${FN}_sorted_PE.bam

{
echo flagstat ${FN}_sorted_PE.bam >> $LOGFILE
samtools flagstat ${FN}_sorted_PE.bam >> $LOGFILE
echo idxstats  ${FN}_sorted_PE.bam >> $LOGFILE
samtools idxstats  ${FN}_sorted_PE.bam >> $LOGFILE
} &


RG='@RG\tID:'$FN'_SE\tSM:'$FN
echo "start SE bwa mem  at" `date` >> $LOGFILE
echo bwa mem -R $RG -k 17 -r 1.25 -M -t 17 $REFGENOME $RM \| samtools view -Shb - \> $FN"_SE.bam"  >> $LOGFILE

bwa mem -R $RG  -k 17 -r 1.25 -M -t 17 $REFGENOME $RM 2>> $ERRORLOG | samtools view -Shb - > $FN"_SE.bam"
ES=$?
echo finished bwa mem mapping at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
echo start sorting at `date`  >> $LOGFILE
echo java  -Xmx6g -jar $PICARD SortSam  VALIDATION_STRINGENCY=SILENT I=${FN}_SE.bam O=${FN}_sorted_SE.bam SO=coordinate  >> $LOGFILE
java  -Xmx6g -jar $PICARD SortSam  VALIDATION_STRINGENCY=SILENT I=${FN}_SE.bam O=${FN}_sorted_SE.bam SO=coordinate >> $LOGFILE 2>> $ERRORLOG
ES=$?
rm -f  $FN"_SE.bam"
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
samtools index ${FN}_sorted_SE.bam
{
echo flagstat ${FN}_sorted_SE.bam >> $LOGFILE
samtools flagstat ${FN}_sorted_SE.bam >> $LOGFILE
echo idxstats  ${FN}_sorted_SE.bam >> $LOGFILE
samtools idxstats  ${FN}_sorted_SE.bam >> $LOGFILE
  } &


# extract only the viral genomes
samtools view -bh -F 256 ${FN}_sorted_SE.bam 'gi|86440167|gb|DQ361066.1|' 'gi|116563461|gb|DQ361065.2|' > ${FN}_sorted_viral_SE.bam
samtools reheader <(samtools view -H ${FN}_sorted_viral_SE.bam | grep -v 'SN:[^g]') ${FN}_sorted_viral_SE.bam > temp.bam
mv -f temp.bam ${FN}_sorted_viral_SE.bam
samtools view -bh -f 2 -F 256 ${FN}_sorted_PE.bam 'gi|86440167|gb|DQ361066.1|' 'gi|116563461|gb|DQ361065.2|' > ${FN}_sorted_viral_PE.bam
samtools reheader <(samtools view -H ${FN}_sorted_viral_PE.bam | grep -v 'SN:[^g]') ${FN}_sorted_viral_PE.bam > temp.bam
mv -f temp.bam ${FN}_sorted_viral_PE.bam
echo java  -Xmx6g -jar $PICARD MergeSamFiles VALIDATION_STRINGENCY=SILENT I=${FN}_sorted_viral_PE.bam I=${FN}_sorted_viral_SE.bam SO=coordinate O=${FN}_sorted_viral_PSE.bam >> $LOGFILE
java  -Xmx6g -jar $PICARD  MergeSamFiles VALIDATION_STRINGENCY=SILENT I=${FN}_sorted_viral_PE.bam I=${FN}_sorted_viral_SE.bam SO=coordinate O=${FN}_sorted_viral_PSE.bam >> $LOGFILE 2>> $ERRORLOG
ES=$?
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
rm -f ${FN}_sorted_viral_SE.bam ${FN}_sorted_viral_PE.bam
{
samtools index ${FN}_sorted_viral_PSE.bam
echo flagstat ${FN}_sorted_viral_PSE.bam >> $LOGFILE
samtools flagstat ${FN}_sorted_viral_PSE.bam >> $LOGFILE
echo idxstats ${FN}_sorted_viral_PSE.bam >> $LOGFILE
samtools idxstats  ${FN}_sorted_viral_PSE.bam >> $LOGFILE

} &
