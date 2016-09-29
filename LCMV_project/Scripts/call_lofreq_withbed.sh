#!/bin/bash
#----------
# author: Lukas Endler
# Time-stamp: <2016-03-15 15:58:01 lukasendler>
# date: 20.9.2015 at 12:23
# takes a bam file and a bed file and calls variants with lofreq2 for specific loci
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses_short.fasta
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
GATK=/Volumes/Temp/Lukas/LCMV_project/Tools/GenomeAnalysisTK-3.4-46.jar
SAMTOOLS=/usr/local/bin/samtools
LOFREQ=$BASEDIR/Tools/lofreq_star-2.1.2/bin/lofreq

for FFN in *real_viterbi_IDQS.bam; do
	FN=`basename $FFN .bam`
	FN=${FN/_real*/}
	LOGFILE=${FN}.log
	ERRORLOG=${FN}.err.log
	echo "start lofreq2 at" `date` >> $LOGFILE
	echo $LOFREQ call -f $REFGENOME --verbose -o ${FN}_lofreq_bed.vcf --bed $1 -q 20 -Q 20 -m 20 -C 75 -a 0.05 --call-indels $FFN >> $LOGFILE
	$LOFREQ call -f $REFGENOME --verbose -o ${FN}_lofreq_bed.vcf --bed $1 -q 20 -Q 20 -m 20 -C 75 -a 0.05 --call-indels $FFN 2>> $ERRORLOG >> $LOGFILE
	ES=$?
	echo finished lofreq at `date` with exit state $ES >> $LOGFILE
done

for i in *_lofreq_bed.vcf; do
   	SMP=${i%_lofreq*}
	awk -v OFS="\t" '/^\#\#[^I]/ {print} /^\#\#INFO/ {print $0; sub("INFO","FORMAT",$0); print $0; } /\#CH/ {print $0,"FORMAT","'$SMP'"} !/^\#/ {form=$8;  gsub(/=[^A-Z]+/,":",form); gsub(/;/,":",form); sub(/:$/,"",form); sub(/INDEL:/,"",form);  samp=$8; gsub(/[A-Z4]+=/,"",samp); gsub(/;/,":",samp); sub(/INDEL:/,"",samp); print $0,form,samp}' $i | bgzip -c > ${SMP}_samp.lofreq.bed.vcf.gz
	tabix -p vcf ${SMP}_samp.lofreq.bed.vcf.gz  
done

bcftools merge -O v *_samp.lofreq.bed.vcf.gz > all_samp_bed.vcf

