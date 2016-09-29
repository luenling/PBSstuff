#!/bin/bash
#----------
# author: Lukas Endler
# date: 20.9.2015 at 16:46
# takes a bam file calls variants with freebayes
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses.fasta
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
GATK=/Volumes/Temp/Lukas/LCMV_project/Tools/GenomeAnalysisTK-3.4-46.jar
SAMTOOLS=/usr/local/bin/samtools
FB=$BASEDIR/Tools/freebayes/bin/freebayes
FN=`basename $1 .list`
MAXCOV=1000000
if [ $2 ] ; then
    MAXCOV=$2
fi
LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log
echo "start freebayes at" `date`
echo $FB --fasta-reference $REFGENOME --hwe-priors-off --pooled-continuous --allele-balance-priors-off \
     --haplotype-length 3 --use-best-n-alleles 3 --bam-list $1 --min-coverage 75\
        --min-alternate-count 3  --min-alternate-fraction 0.001 \
        --min-mapping-quality 20 --min-base-quality 30  >> $LOGFILE
$FB --fasta-reference $REFGENOME --hwe-priors-off --pooled-continuous --allele-balance-priors-off \
     --haplotype-length 3 --use-best-n-alleles 3 --bam-list $1 --min-coverage 75\
        --min-alternate-count 3  --min-alternate-fraction 0.001 \
        --min-mapping-quality 20 --min-base-quality 30  > ${FN}_freebayes.vcf 2>> $ERRORLOG
        echo "finished with freebayes at"`date`
