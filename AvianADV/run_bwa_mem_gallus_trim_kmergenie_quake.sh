#!/bin/bash
#----------
# author: Lukas Endler
# Time-stamp: <2016-03-15 14:10:21 lukasendler>
# date: 03.03.2016 at 20:31
# takes two fastq files runs bwa with 12 threads and outputs a bam file called like the fastq prefix against gallus gallus
# only keeps the non mapped reads and puts them into fastq files
# trims them
# runs kmergenie
# runs quake for kmer correction
# call with bash command blub_1.fq blub_2.fq >> logfile.log 2>> log.error.log
#--------------
BASEDIR=/Volumes/Temp/Anna.Schachner
REFGENOME=$BASEDIR/References/Gallus_gallus.Galgal4.dna.toplevel.fa.gz
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
BBDUK=/Volumes/Temp/Lukas/LCMV_project/Tools/bbmap/bbduk.sh
BBMERGE=/Volumes/Temp/Lukas/LCMV_project/Tools/bbmap/bbmerge.sh
BWA=/usr/local/Cellar/bwa/0.7.12/bin/bwa
BFC=$BASEDIR/Tools/bfc/bfc
ADAPTS=/Volumes/Temp/Lukas/LCMV_project/Tools/bbmap/resources/adapters.fa
KMERGENIE=/Volumes/Temp/Lukas/Tools/kmergenie-1.6741/kmergenie

FN=`basename $1 .gz`
FN=`basename $FN _1.fq`
FN=`basename $FN _1.fastq`
# just to replace the whole bloody name tag
R2=`echo $1 | sed 's/_1\.f/_2.f/'`
FN=${FN/*_H325VBBXX_5_/}
SAMPLE=$FN
LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log

RG='@RG\tID:'$FN'\tSM:'$FN
echo "start bwa mem  at" `date` >> $LOGFILE
echo $BWA mem -R $RG -M -t 12 $REFGENOME $1 $R2 \| samtools view -Shb - \> $FN".bam"  >> $LOGFILE

$BWA mem -R $RG -M -t 12 $REFGENOME $1 $R2 2>> $ERRORLOG | samtools view -Shb - > $FN".bam"
ES=$?
echo finished bwa mem mapping at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
echo "remove gallus gallus at " `date` >> $LOGFILE
echo $SAMTOOLS view -bf 13  -F 256 $FN".bam" \|  $SAMTOOLS bam2fq - \| gzip -c - \> ${FN}_nogal_12.fq.gz >> $LOGFILE
$SAMTOOLS view -b -f 13  -F 256 $FN".bam" |  $SAMTOOLS bam2fq - | gzip -c - > ${FN}_nogal_12.fq.gz
echo flagstat >> $LOGFILE
samtools flagstat ${FN}.bam >> $LOGFILE &

echo trimming and removing adaptors at `date` >> ${LOGFILE}
echo $BBDUK -Xmx4g in1=${FN}_nogal_12.fq.gz out1=${FN}_trim_1.fq out2=${FN}_trim_2.fq qtrim=r minlen=100 ktrim=r ref=\"$ADAPTS\" k=20 mink=15 trimq=15 overwrite=t tbo=t tpe=t lhist=${FN}_rlen.hist stats=${FN}_trim.stats \>\> ${LOGFILE} 2\> ${ERRORLOG} >> ${LOGFILE}
$BBDUK -Xmx4g in1=${FN}_nogal_12.fq.gz out1=${FN}_trim_1.fq out2=${FN}_trim_2.fq qtrim=r minlen=100 ktrim=r ref=\"$ADAPTS\" k=20 mink=15 trimq=15 overwrite=t tbo=t tpe=t lhist=${FN}_rlen.hist stats=${FN}_trim.stats  >> ${LOGFILE} 2> ${ERRORLOG}
ES=$?
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES

fastq-stats ${FN}_trim_2.fq > ${FN}_trim_2.fq.stats &
fastq-stats ${FN}_trim_1.fq > ${FN}_trim_1.fq.stats &

#$KMERGENIE ${FN}_trim_1.fq.gz -k 60 -l 20 -s 10 -t 12 -o ${FN}_kmergenie/${FN}_kmergenie  >> ${LOGFILE} 2> ${ERRORLOG}

# downsample reads to coverage ~ 200 assuming mean length of 200 and 50K genome:
SEED=$RANDOM
echo subsampling at `date`  >> ${LOGFILE}
echo fastq-sample -s $SEED -o ${FN}_trim_sub_1 -n 50000 ${FN}_trim_1.fq >> ${LOGFILE}
fastq-sample -s $SEED -o ${FN}_trim_sub_1 -n 50000 ${FN}_trim_1.fq
echo fastq-sample -s $SEED -o ${FN}_trim_sub_2 -n 100000 ${FN}_trim_2.fq >> ${LOGFILE}
fastq-sample -s $SEED -o ${FN}_trim_sub_2 -n 50000 ${FN}_trim_2.fq

# create read file for quake
echo echo \-e \"${FN}_trim_sub_1.fastq\t${FN}_trim_sub_2.fastq\" \> reads.txt  >> $LOGFILE
echo -e "${FN}_trim_sub_1.fastq\t${FN}_trim_sub_2.fastq" > reads.txt
echo python ~/Tools/Quake/bin/quake.py -k 17 --ratio=100 -p 5 -q 33 -f reads.txt >> $LOGFILE
python /Volumes/Temp/Lukas/Tools/Quake/bin/quake.py -k 17 --ratio=100 -p 5 -q 33 -f reads.txt
if [ ! -e  ${FN}_trim_sub_1.cor.fastq ]
    then
    echo /Volumes/Temp/Lukas/Tools/Quake/bin/correct -f reads.txt -k 16 -m reads.txt.qcts c-p 10 -q 33 -b reads.dbm -c 2 >> $LOGFILE
    /Volumes/Temp/Lukas/Tools/Quake/bin/correct -f reads.txt -k 17 -m reads.txt.qcts c-p 10 -q 33 -b reads.dbm -c 2
fi
rm -f reads.dbm
#quorum -s 1G -q 33 -t 10 -p ${FN}_trim_sub_quorum ${FN}_trim_sub_1.fastq ${FN}_trim_sub_2.fastq
SPADES=/Volumes/Temp/Anna.Schachner/Tools/Spades/SPAdes-3.7.0/spades.py
echo $SPADES -k 21,33,55,77,99,127 --careful --only-assembler -1 ${FN}_trim_sub_1.cor.fastq -2 ${FN}_trim_sub_2.cor.fastq -o ${FN}_Spades \>\> ${FN}_spades.log 2\>\> ${FN}_spades.error.log  >> $LOGFILE
$SPADES -k 21,33,55,77,99,127 --careful --only-assembler -1 ${FN}_trim_sub_1.cor.fastq -2 ${FN}_trim_sub_2.cor.fastq -o ${FN}_Spades >> ${FN}_spades.log 2>> ${FN}_spades.error.log

abyss-fac ${FN}_Spades/*.fasta >> $LOGFILE

# map reads against contigs
echo $BWA index ${FN}_Spades/scaffolds.fasta >> $LOGFILE
$BWA index ${FN}_Spades/scaffolds.fasta
echo $BWA mem -t 10 ${FN}_Spades/scaffolds.fasta ${FN}_trim_1.fq  ${FN}_trim_2.fq \| samtools view -Shb - \| $SAMTOOLS sort -T ${FN}_to_quake - \>\> ${FN}_nogal_to_quake_sorted.bam >> $LOGFILE
$BWA mem -t 10 ${FN}_Spades/scaffolds.fasta ${FN}_trim_1.fq  ${FN}_trim_2.fq | samtools view -Shb - | $SAMTOOLS sort -T ${FN}_to_quake - >> ${FN}_nogal_to_quake_sorted.bam
$SAMTOOLS index ${FN}_nogal_to_quake_sorted.bam
echo $SAMTOOLS flagstat  ${FN}_nogal_to_quake_sorted.bam >> $LOGFILE
$SAMTOOLS flagstat  ${FN}_nogal_to_quake_sorted.bam >> $LOGFILE
echo $SAMTOOLS idxstats  ${FN}_nogal_to_quake_sorted.bam >> $LOGFILE
$SAMTOOLS idxstats  ${FN}_nogal_to_quake_sorted.bam >> $LOGFILE


exit

# merge fastq to interleaved
echo fastqutils merge -slash ${FN}_trim_sub_1.fastq ${FN}_trim_sub_2.fastq \> ${FN}_trim_sub_12.fastq  >> ${LOGFILE}
fastqutils merge -slash ${FN}_trim_sub_1.fastq ${FN}_trim_sub_2.fastq > ${FN}_trim_sub_12.fq

# error correction
echo $BFC -w 25 -s 50k -t 10 -c 10 ${FN}_trim_sub_12.fq \| awk '{print $1}' \> ${FN}_trim_sub_ec_12.fq  >> ${LOGFILE}
$BFC -w 15 -s 50k -t 10 -c 6  ${FN}_trim_sub_12.fq | awk '{print $1}' > ${FN}_trim_sub_ec_12.fq
fastqutils unmerge ${FN}_trim_sub_ec_12.fq ${FN}_trim_sub_ec

$SPADES -k 21,33,55,77,99,127 --careful --only-assembler -1 ${FN}_trim_sub_ec.1.fastq -2 ${FN}_trim_sub_ec.2.fastq -o ${FN}_Spades >> ${FN}_spades.log 2>> ${FN}_spades.error.log


# join reads using bbmerge
if [ ! -d Merged ]; then mkdir Merged ; fi
echo merging overlapping reads at `date` >> ${LOGFILE}
echo $BBMERGE -Xmx4g  in=${FN}_trim_sub_ec_12.fq out=Merged/${FN}_trim_sub_ec_m.fq.gz outu1=Merged/${FN}_trim_sub_ec_u_1.fq.gz outu2=Merged/${FN}_trim_sub_ec_u_2.fq.gz strict=t ihist=Merged/${FN}_trim_sub_ec_m.hist 2\> Merged/${FN}_trim_sub_ec_m.log >>  $LOGFILE
$BBMERGE -Xmx4g  in=${FN}_trim_sub_ec_12.fq out=Merged/${FN}_trim_sub_ec_m.fq outu1=Merged/${FN}_trim_sub_ec_u_1.fq outu2=Merged/${FN}_trim_sub_ec_u_2.fq strict=t ihist=Merged/${FN}_trim_sub_ec_m.hist 2> Merged/${FN}_trim_sub_ec_m.log
ES=$?
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES

# joining using abyss
abyss-pe name=${FN} k=40 lib="pe1" pe1="../Merged/${FN}_trim_sub_ec_u_1.fq ../Merged/${FN}_trim_sub_ec_u_2.fq" se="../Merged/${FN}_trim_sub_ec_m.fq"

export k
for k in 20 25 30 35 40 45 50 55 60 64; do
    mkdir k$k
    abyss-pe -C k$k name=${FN} k=$k j=10 q=15 lib="pe1" pe1="../../Merged/${FN}_trim_sub_ec_u_1.fq ../../Merged/${FN}_trim_sub_ec_u_2.fq" se="../../Merged/${FN}_trim_sub_ec_m.fq" >> k_${FN}.log
done
abyss-fac k*/${FN}-contigs.fa

#quorum -s 1G -q 33 -t 10 -p ${FN}_trim_sub_quorum ${FN}_trim_1_sub.fq ${FN}_trim_2_sub.fq
#flash --interleaved-input sample_id_35587_trim_ec_12.fq.gz -t 10 -z -d Flash --max-overlap 230 -p 33 -o sample_id_35587_trim_ec_flash -x 0.1
#flash -t 10 -d Flash_quake --max-overlap 300 -p 33 -o ${FN}_trim_flash_quake -x 0.1 ${FN}_trim_sub_1.cor.fq ${FN}_trim_sub_2.cor.fq >> flash.log
#flash -t 10 -d Flash_noec --max-overlap 300 -p 33 -o ${FN}_trim_flash -x 0.1 ${FN}_trim_sub_1.fastq ${FN}_trim_sub_2.fastq >> flash.log
echo echo \-e \"${FN}_trim_sub_1.fq\t${FN}_trim_sub_2.fq\" \> reads.fq  >> $LOGFILE

echo -e "${FN}_trim_sub_1.fastq\t${FN}_trim_sub_2.fastq" > reads.fq
echo python ~/Tools/Quake/bin/quake.py -k 17 -p 10 -q 33 -f reads.fq >> $LOGFILE
python ~/Tools/Quake/bin/quake.py -k 17 -p 5 -q 33 -f reads.fq
#quorum -s 1G -q 33 -t 10 -p ${FN}_trim_sub_quorum ${FN}_trim_sub_1.fastq ${FN}_trim_sub_2.fastq
SPADES=/Volumes/Temp/Anna.Schachner/Tools/Spades/SPAdes-3.7.0/spades.py
$SPADES -k 21,33,55,77,99,127 --careful --only-assembler -1 ${FN}_trim_sub_1.cor.fq -2 ${FN}_trim_sub_2.cor.fq -o Spade_Quake

 abyss-pe name=${FN} k=40 j=10 q=15 lib="pe1" pe1="../../${FN}_trim_sub_1.cor.fq ../../${FN}_trim_sub_2.cor.fq" >> k_${FN}.log

flash -t 10 -d Flash_${FN} --max-overlap 300 -p 33 -o ${FN}_trim_flash -x 0.1 ${FN}_trim_1.fq ${FN}_trim_2.fq >> ${FN}_flash.log
 fastq-sample -s $SEED -o  Flash_${FN}/${FN}_trim_flash_sub -n 10000 Flash_${FN}/${FN}_trim_flash.extendedFrags.fastq
SPADES=/Volumes/Temp/Anna.Schachner/Tools/Spades/SPAdes-3.7.0/spades.py
$SPADES -k 21,33,55,77,99,127 --careful --only-assembler --pe1-1 ${FN}_trim_sub_1.cor.fastq --pe1-2 ${FN}_trim_sub_2.cor.fastq --s1 Flash_${FN}/${FN}_trim_flash_sub.fastq  -o Spade_Quake_flash

$SPADES -k 21,33,55,77,99,111 --sc --careful --only-assembler --sc --pe1-1 ${FN}_trim_sub_1.cor.fastq --pe1-2 ${FN}_trim_sub_2.cor.fastq -  -o Spade_Quake_sc
 echo ../Flash_sample_id_35580/sample_id_35580_trim_flash_sub.fastq > flashreads ; LINKS -f contigs.fasta -d 550 -s flashreads -k 17 -e 0.5 -a 0.75
