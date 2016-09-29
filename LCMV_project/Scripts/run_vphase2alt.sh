#!/bin/bash
#----------
# author: Lukas Endler
# date: 17.9.2015 at 16:34
# takes two fastq files runs bwa with 12 threads and outputs a bam file called like the fastq prefix
# call with bash command blub_1.fq blub_2.fq >> logfile.log 2>> log.error.log
# also sort the file afterwards
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses_short.fasta
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
VPHASER2=$BASEDIR/Tools/viral-ngs/intrahost.py
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$BASEDIR/Tools/VPhaser-2-02112013/bamtools-2.3.0/lib

while read FN; do
	FNB=`basename $FN .bam`
	SAMP=${FNB/_nopr*/}
	SAMPLES=${SAMPLES}${SAMP}" "
	ALGS=${ALGS}$FN" "
	ISNVS=${ISNVS}${FNB}.tab" "
	LOGFILE=${FNB}.vphaser2.log
	ERRORLOG=${FNB}.vphaser2.err.log
	OMP_NUM_THREADS=8
	#echo "start VPHASER2 for " $FN " at" `date` >> $LOGFILE
	#echo $VPHASER2 vphaser_one_sample   --vphaserNumThreads 10 --minReadsEach 2 $FN $REFGENOME ${FNB}.tab >> $LOGFILE
	#$VPHASER2 vphaser_one_sample  --vphaserNumThreads 10 --minReadsEach 2 $FN $REFGENOME ${FNB}.tab 2>> $ERRORLOG >> $LOGFILE
	#ES=$?
	#echo finished vphaser2 at `date` with exit state $ES >> $LOGFILE
done <$1

echo " merge to vcf " $1 " at" `date` >> $LOGFILE
echo $VPHASER2 merge_to_vcf --samples $SAMPLES --isnvs $ISNVS --alignments $ALGS $REFGENOME all_out.vcf >> all.log
$VPHASER2 merge_to_vcf --samples "$SAMPLES" --isnvs "$ISNVS" --alignments "$ALGS" $REFGENOME all_out.vcf 2>> all.err.log >> all.log

exit 0
if [ ! -d VPHASER2/$FN ]; then
    echo mkdir >> $LOGFILE
    mkdir -p VPHASER2/$FN
fi

/Volumes/Temp/Lukas/LCMV_project/Tools/viral-ngs/intrahost.py merge_to_vcf  /Volumes/Temp/Lukas/LCMV_project/References/viruses_short.fasta all_out.vcf --samples "ND_0,ND_ARM_infected_BHK21_S11759" --isnvs "ND_0_noprime_sorted_viral_sh_max_cov_10000_short_RG.tab,ND_ARM_infected_BHK21_S11759_noprime_sorted_viral_sh_max_cov_10000_short_RG.tab"  --alignments "../ND_0_noprime_sorted_viral_sh_max_cov_10000_short_RG.bam,../ND_ARM_infected_BHK21_S11759_noprime_sorted_viral_sh_max_cov_10000_short_RG.bam"

