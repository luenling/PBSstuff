#!/bin/bash
#----------
# author: Lukas Endler
# authored: 2016-05-02 12:01:03 
# Time-stamp: <2016-07-11 14:35:50 lukasendler>
# command lines
# usage: command DATADIR BAMSUFFIX
# default: DATADIR: /Volumes/Temp/Hammer/data/ BAMSUFFIX: _standardMapQual.Aligned.sorted.bam
#--------------

DATA=/Volumes/Temp/Hammer/data/
SCRIPTS=/Volumes/Temp/Hammer/scripts/
CUFFDIR=/Volumes/Temp/Lukas/Tools/cufflinks-2.2.1/src
GTF=/Volumes/Temp/Hammer/annotations/Sus_scrofa.Sscrofa10.2.84.with.chr.no_gene.gtf
GTF_DIR=/Volumes/Temp/Hammer/annotations/
REFGENOME=/Volumes/Temp/Hammer/annotations/susScr3.fa
SUFFIX=_standardMapQual.Aligned.sorted.bam
# BAM=$1
# BN=`basename $1 .bam`
# LOGFILE=${BN}".log"

if [ -d $1 ]
then
    DATA=$1
fi


OUTDIR=$DATA/Cufflinks/
mkdir -p $OUTDIR

if [ $2 ]
then
    SUFFIX=$2
fi

LOGFILE=${OUTDIR}/"logfile.log"

function add_chr_correctly () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME on $BN at "`date`
    echo "-----------------------------------------------------------------"
    sed -E '/^[0-9XYM]/ s/^/chr/g; s/chrMT/chrM/g; s/(^[A-Z0-9]*)\.([1234])/\1-\2/g ' $GTF_DIR/Sus_scrofa.Sscrofa10.2.84.gtf > $GTF_DIR/Sus_scrofa.Sscrofa10.2.84.chr_added_correctly.gtf
    awk '$3 != "gene" {print}' $GTF_DIR/Sus_scrofa.Sscrofa10.2.84.chr_added_correctly.gtf > $GTF_DIR/Sus_scrofa.Sscrofa10.2.84.with.chr.no_gene.gtf
}


function run_cufflinks () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME on $BN at "`date`
    echo "-----------------------------------------------------------------"
    echo "starting cufflinks at:"`date`  >> $LOGFILE
    echo  ${CUFFDIR}/cufflinks -q -o ${OUTDIR}/${BN}_cufflinks -p 12  --min-isoform-fraction 0.01 --small-anchor-fraction 0.05 -g $GTF --library-type fr-firststrand $BAM 2\> ${OUTDIR}/${BN}_cufflinks.log  >> $LOGFILE
    ${CUFFDIR}/cufflinks -q -o ${OUTDIR}/${BN}_cufflinks -p 12  --min-isoform-fraction 0.01 --small-anchor-fraction 0.05 -g $GTF --library-type fr-firststrand $BAM 2> ${OUTDIR}/${BN}_cufflinks.log
    ES=$?
    echo finished at `date` with exit state $ES >> $LOGFILE
    [ $ES -eq 0 ] || exit $ES
    echo finished at `date` with exit state $ES

}

function check_bam_index () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME on $BAM at "`date`
    echo "-----------------------------------------------------------------"
    if [ ! -e $BAM".bai" ]
    then
	echo "indexing $BAM at `date`" >> $LOGFILE
	echo samtools index $BAM  >> $LOGFILE
	samtools index $BAM
	ES=$?
	echo finished at `date` with exit state $ES >> $LOGFILE
	[ $ES -eq 0 ] || exit $ES
    fi    
}

function run_cuffcompare () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME "`date`
    echo "-----------------------------------------------------------------"
    echo "starting cufflinks at:"`date`  >> $LOGFILE
    echo  ${CUFFDIR}/cuffcompare  -o ${OUTDIR}/cuffcompare -r $GTF -CR `echo ${GTFS[@]}` 2\> ${OUTDIR}/cuffcompare.log  >> $LOGFILE
   ${CUFFDIR}/cuffcompare  -o ${OUTDIR}/cuffcompare  -r $GTF -CR `echo ${GTFS[@]}` 2> ${OUTDIR}/cuffcompare.log
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
    rm -f ${OUTDIR}/gtfs_to_merge.txt
    for i in ${GTFS_CLEANED[@]}; do echo $i >>  ${OUTDIR}/gtfs_to_merge.txt; done
    # need to create folder and then change to it first as option -o gives error with cuffmerge
    mkdir ${OUTDIR}/cuffmerge
    CUR=`pwd`
    cd  ${OUTDIR}/cuffmerge
    echo  ${CUFFDIR}/cuffmerge -p 16 -g $GTF -s $REFGENOME ${OUTDIR}/gtfs_to_merge.txt 2\> ${OUTDIR}/cuffmerge.log  >> $LOGFILE
    ${CUFFDIR}/cuffmerge -p 16 -g $GTF -s $REFGENOME ${OUTDIR}/gtfs_to_merge.txt 2>> ${OUTDIR}/cuffmerge.log
    ES=$?
    cd $CUR
    echo finished at `date` with exit state $ES >> $LOGFILE
    [ $ES -eq 0 ] || exit $ES
    echo finished at `date` with exit state $ES

}

function run_cuffquant () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME on $BN at "`date`
    echo "-----------------------------------------------------------------"
    echo "starting cuffnorm at:"`date`  >> $LOGFILE
    # need to create folder and then change to it first as option -o gives error with cuffmerge
    mkdir -p ${OUTDIR}/cuffquant
    CUR=`pwd`
    cd  ${OUTDIR}/cuffquant
    echo  ${CUFFDIR}/cuffnorm -p 16 --library-type fr-firststrand -q  ${OUTDIR}/cuffmerge/merged_asm/merged.gtf  `ls $DATA/*$SUFFIX` 2\>\> ${OUTDIR}/cuffnorm.log  >> $LOGFILE
    ${CUFFDIR}/cuffnorm -p 16 --library-type fr-firststrand -q  ${OUTDIR}/cuffmerge/merged_asm/merged.gtf  `ls $DATA/*$SUFFIX` 2> ${OUTDIR}/cuffnorm.log
    ES=$?
    cd $CUR
    echo finished at `date` with exit state $ES >> $LOGFILE
    [ $ES -eq 0 ] || exit $ES
    echo finished at `date` with exit state $ES

}

function add_info_to_cuffnorm () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME on $BN at "`date`
    echo "-----------------------------------------------------------------"
    echo "starting add_info_to_cuffnorm at:"`date`  >> $LOGFILE
    # adding the info from the attr_table files to the cuffnorm output and change the sample names
    CUR=`pwd`
    cd ${OUTDIR}/cuffquant
    sub1=""; sub2="";
    while read line ; do
	read -r -a array <<< "$line";
	sub1=${sub1}"\\t"${array[0]};
	sn=`basename ${array[1]} _Aligned.sorted.bam`;
	sub2=${sub2}"\\t"${sn};
    done < <(tail -8 samples.table)
    for i in cds genes isoforms tss_groups ;
    do
	join -t $'\t' -1 1 -2 1 ${i}.count_table ${i}.attr_table > ${i}.count_table.ext ;
	join -t $'\t' -1 1 -2 1 ${i}.fpkm_table ${i}.attr_table > ${i}.fpkm_table.ext;
	gsed -i.bak '1 s/'$sub1'/'$sub2'/' ${i}.count_table.ext 
	gsed -i.bak '1 s/'$sub1'/'$sub2'/' ${i}.fpkm_table.ext 	
    done
    unset sub1
    unset sub2
    cd $CUR
    echo finished at `date` >> $LOGFILE

}





# index genome
if [ ! -e $REFGENOME".fai" ]
    then
	echo "indexing $REFGENOME at `date`" >> $LOGFILE
	echo samtools faidx  $REFGENOME  >> $LOGFILE
	samtools faidx  $REFGENOME
	ES=$?
	echo finished at `date` with exit state $ES >> $LOGFILE
	[ $ES -eq 0 ] || exit $ES
    fi    


# set BN array
BNS=()
# set GTFS array
GTFS=()
GTFS_CLEANED=()
for BAM in ${DATA}/*$SUFFIX
do
    BN=`basename $BAM .bam`
    BNS+=($BN)
    LOGFILE=${OUTDIR}/${BN}".log"
    check_bam_index
    #run_cufflinks
    GTFS+=(${OUTDIR}/${BN}_cufflinks/transcripts.gtf)
    if [ ! -e  ${OUTDIR}/${BN}_cufflinks/${BN}_transcripts.gtf  ] ||  \
	       [ ${OUTDIR}/${BN}_cufflinks/${BN}_transcripts.gtf -ot ${OUTDIR}/${BN}_cufflinks/transcripts.gtf  ]
    then
	echo  "cleaning gtf at at `date`" >> $LOGFILE
	echo python $SCRIPTS/clean_gtf.py -g ${OUTDIR}/${BN}_cufflinks/transcripts.gtf -f  $REFGENOME".fai" \> ${OUTDIR}/${BN}_cufflinks/${BN}_transcripts.gtf >> $LOGFILE
	python $SCRIPTS/clean_gtf.py -g ${OUTDIR}/${BN}_cufflinks/transcripts.gtf -f  $REFGENOME".fai" > ${OUTDIR}/${BN}_cufflinks/${BN}_transcripts.gtf
    fi
    GTFS_CLEANED+=(${OUTDIR}/${BN}_cufflinks/${BN}_transcripts.gtf)
done



#run_cuffcompare
#run_cuffmerge
run_cuffquant
