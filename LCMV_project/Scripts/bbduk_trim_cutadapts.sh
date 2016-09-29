#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
BASEDIR=/Volumes/Temp/Lukas/LCMV_project/
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
BBDUK=$BASEDIR/Tools/bbmap/bbduk.sh
BBMERGE=$BASEDIR/Tools/bbmap/bbmerge.sh
ADAPTS=$BASEDIR/Tools/bbmap/resources/adapters.fa
PRIMERS=$BASEDIR/References/primers.fna
#REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
INF1=$1
INF2=`echo $INF1 | sed 's/_1\.f/_2.f/' `
OUT=`basename $INF1 .gz`
OUT=`basename $OUT .fastq`
OUT=`basename $OUT .fq`
OUT=`basename $OUT _1`
OUT=${OUT/BSF_0201_000000000-AMA24_1#/ND_}

LOGFILE=$OUT.log
ERRORFILE=$OUT.err.log
TRIMLOG=$OUT.trim.log
OUTDIR=''
if [ $2 ]; then
    OUTDIR=$2'/'
    if [ ! -d $2 ]; then mkdir $3 ; fi
fi

echo removing primers at `date` >> ${LOGFILE}
echo $BBDUK -Xmx4g in1=$INF1 in2=$INF2 out1=${OUT}_np_1.fq.gz out2=${OUT}_np_2.fq.gz  ref=\"$PRIMERS\" minlen=55 ktrim=l k=20 mink=12 overwrite=t restrictleft=22  \>\> ${LOGFILE} 2\>\> ${ERRORFILE} >> ${LOGFILE}
$BBDUK -Xmx4g in1=$INF1 in2=$INF2 out1=${OUT}_np_1.fq.gz out2=${OUT}_np_2.fq.gz  ref=\"$PRIMERS\" minlen=75 ktrim=l k=20 mink=12 overwrite=t restrictleft=22  >> ${LOGFILE} 2>> ${ERRORFILE}
ES=$?
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES

echo trimming and removing adaptors at `date` >> ${LOGFILE}
echo $BBDUK -Xmx4g in1=${OUT}_np_1.fq.gz in2=${OUT}_np_2.fq.gz out1=${OUT}_npa_1.fq.gz out2=${OUT}_npa_2.fq.gz qtrim=r minlen=75 ktrim=r ref=\"$ADAPTS\" k=20 mink=15 trimq=15 overwrite=t tbo=t tpe=t \>\> ${LOGFILE} 2\> ${ERRORFILE} >> ${LOGFILE}
$BBDUK -Xmx4g in1=${OUT}_np_1.fq.gz in2=${OUT}_np_2.fq.gz out1=${OUT}_npa_1.fq.gz out2=${OUT}_npa_2.fq.gz qtrim=r minlen=55 ktrim=r ref=\"$ADAPTS\" k=20 mink=15 trimq=15 overwrite=t tbo=t tpe=t >> ${LOGFILE} 2> ${ERRORFILE}
ES=$?
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
if [ ! -d fastqc ]; then mkdir fastqc ; fi
fastqc -q -o fastqc ${OUT}_npa_1.fq.gz ${OUT}_npa_2.fq.gz ${OUT}_np_1.fq.gz ${OUT}_np_2.fq.gz &
if [ ! -d Merged ]; then mkdir Merged ; fi
echo merging overlapping reads at `date` >> ${LOGFILE}
echo $BBMERGE -Xmx4g  in1=${OUT}_npa_1.fq.gz in2=${OUT}_npa_2.fq.gz out=Merged/${OUT}_m.fq.gz outu1=Merged/${OUT}_m_1.fq.gz outu2=Merged/${OUT}_m_2.fq.gz strict=t qtrim2=r trimq=10 ouq=t ihist=Merged/${OUT}_m.hist 2\> Merged/${OUT}_m.log >>  $LOGFILE
$BBMERGE -Xmx4g  in1=${OUT}_npa_1.fq.gz in2=${OUT}_npa_2.fq.gz out=Merged/${OUT}_m.fq.gz outu1=Merged/${OUT}_m_1.fq.gz outu2=Merged/${OUT}_m_2.fq.gz strict=t qtrim2=r trimq=10 ouq=t ihist=Merged/${OUT}_m.hist 2> Merged/${OUT}_m.log
ES=$?
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES





