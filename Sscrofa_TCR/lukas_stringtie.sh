#!/bin/bash
#----------
# author: Lukas Endler
# authored: 2016-05-02 12:01:03 
# Time-stamp: <2016-07-11 15:19:30 lukasendler>
# command lines for cufflinks
#--------------

DATA=/Volumes/Temp/Hammer/data/
OUTDIR=/Volumes/Temp/Hammer/Stringtie/
STRINGTIE=/Volumes/Temp/Hammer/Tools/stringtie-1.2.3.OSX_x86_64/stringtie


function add_XS_to_bam () {
    samtools view -h ../../data/Lung_standardMapQual.Aligned.sorted.bam chr1:102394113-107658852 | awk -v OFS="\t" ' /^@/ {print} !/^@/ {if ($2 AND 64) print $0,"XS:A:-"; if  ($2 AND 128) print $0,"XS:A:-"   } ' | samtools view -Sb - > Lung_standardMapQual.Aligned.sorted.XS.chr1_102394113_107658852.bam

}

$STRINGTIE -p 10  -f 0.01 -G $GTF Lung_standardMapQual.Aligned.sorted.XS.chr1_102394113_107658852.bam


	GTF=/Volumes/Temp/Hammer/annotations/Sus_scrofa.Sscrofa10.2.84.with.chr.BAM=$1
    echo  ${CUFFDIR}/cuffcompare  -o ${OUTDIR}/cuffcompare -r gtf -CR `echo ${GTFS[@]}` 2\> ${OUTDIR}/cuffcompare.log  >> $LOGFILE
   ${CUFFDIR}/cuffcompare  -o ${OUTDIR}/cuffcompare -r gtf -CR `echo ${GTFS[@]}` 2> ${OUTDIR}/cuffcompare.log
    ES=$?
    echo finished at `date` with exit state $ES >> $LOGFILE
    [ $ES -eq 0 ] || exit $ES
    echo finished at `date` with exit state $ES

}

function run_cuffmerge () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME on $BN at "`date`
    echo "-----------------------------------------------------------------"
    echo "starting cuffmerge at:"`date`  >> $LOGFILE
    for i in ${GTFs[@]}; do echo $i >>  ${OUTDIR}/gtfs_to_merge.txt; done
    echo  ${CUFFDIR}/cuffmerge -p 12 -o ${OUTDIR}/cuffmerge -g $GTF -s $REFGENOME ${OUTDIR}/gtfs_to_merge.txt 2\> ${OUTDIR}/cuffmerge.log  >> $LOGFILE
     ${CUFFDIR}/cuffmerge -p 12 -o ${OUTDIR}/cuffmerge -g $GTF -s $REFGENOME ${OUTDIR}/gtfs_to_merge.txt 2>> ${OUTDIR}/cuffmerge.log
    ES=$?
    echo finished at `date` with exit state $ES >> $LOGFILE
    [ $ES -eq 0 ] || exit $ES
    echo finished at `date` with exit state $ES

}

# set BN array
BNS=()
# set GTFS array
GTFS=()
for BAM in ${DATA}/*_standardMapQual.Aligned.sorted.bam
do
    BN=`basename $BAM .bam`
    BNS+=($BN)
    #LOGFILE=${OUTDIR}/${BN}".log"
    #check_bam_index
    #run_cufflinks
    GTFS+=(${OUTDIR}/${BN}_cufflinks/transcripts.gtf)
done

run_cuffcompare
run_cuffmerge
