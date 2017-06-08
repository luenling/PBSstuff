#!/bin/bash
#----------
# author: Lukas Endler
# Time-stamp: <2017-01-25 15:41:21 lukasendler>
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
	awk -v OFS="\t" '/^\#\#[^I]/ {print} /^\#\#INFO/ {sub("AF,Number=1","AF,Number=A",$0); print $0; sub("INFO","FORMAT",$0); print $0; print "##FORMAT=<ID=PQ,Number=1,Type=Integer,Description=\"Phred-scaled variant call P value\">" } /\#CH/ {print $0,"FORMAT","'$SMP'"} !/^\#/ {form=$8;  gsub(/=[^A-Z]+/,":",form); gsub(/;/,":",form); sub(/:$/,"",form); sub(/INDEL:/,"",form);  samp=$8; gsub(/[A-Z4]+=/,"",samp); gsub(/;/,":",samp); sub(/INDEL:/,"",samp); print $0,form":PQ",samp":"$6}' $i | sed 's/_S19[0-9]*_npa//g; s/_S19[0-9]*//g'  | bgzip -c > ${SMP}_samp.lofreq.bed.vcf.gz
	tabix -p vcf ${SMP}_samp.lofreq.bed.vcf.gz  
done

bcftools merge -i "DP:sum,DP4:sum,AF:max,SB:max"  -m none -O v *_samp.lofreq.bed.vcf.gz > all_samp_bed.vcf

DPS=$2

if [[ -e $DPS ]];
then
    mv all_samp_bed.vcf all_samp_bed_tmp.vcf
    bgzip all_samp_bed_tmp.vcf
    tabix -p vcf all_samp_bed_tmp.vcf.gz
    ~/Tools/bcftools-1.3.1/bcftools annotate -a $DPS -c "+FORMAT/DP" all_samp_bed_tmp.vcf.gz > all_samp_bed.vcf
fi


bcftools norm -f $REFGENOME -m+any  all_samp_bed.vcf >  lofreq2_all_samp_bed_norm.vcf



for i in 0.1 0.05 0.001 ;
do
    LFVCF=lofreq2_all_samp_bed_norm_${i}.vcf
    bcftools view -i "AF>$i" all_samp_bed.vcf | bcftools norm -f $REFGENOME -m+any - > $LFVCF
    FN=`basename $LFVCF .vcf`
    LOGFILE=${FN}.log
    ERRORLOG=${FN}.err.log


    echo startingrunning snpeff at `date` >> $LOGFILE
    echo  snpeff lcmv -no-intergenic -no "INTRAGENIC" -no-downstream -no-upstream -stats ${FN}_snpeff_log.html $LFVCF \> ${FN}_snpeff.vcf >> $LOGFILE
    snpeff lcmv -no-intergenic  -no "INTRAGENIC" -no-downstream -no-upstream -stats ${FN}_snpeff_log.html $LFVCF > ${FN}_snpeff.vcf

    ES=$?
    echo finished snpeff at `date` with exit state $ES >> $LOGFILE
    [ $ES -eq 0 ] || exit $ES

    java -jar $BASEDIR/Tools/snpEff/SnpSift.jar extractFields ${FN}_snpeff.vcf CHROM POS REF "ANN[0].GENE" "ANN[0].ALLELE" "ANN[0].EFFECT" "ANN[0].AA" "ANN[1].ALLELE" "ANN[1].EFFECT" "ANN[1].AA" > ${FN}_snpeff.tab
    bcftools query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%SB\t%DP\t%AF\t%PQ]\n"  $LFVCF > ${FN}.stats.tab
    bcftools query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%AF]\n"  $LFVCF > ${FN}.afs.tab
    bcftools query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%DP]\n"  $LFVCF > ${FN}.dp.tab
    bcftools query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%PQ]\n"  $LFVCF > ${FN}.pq.tab
    paste ${FN}.afs.tab <(cut -f5- ${FN}.dp.tab) <(cut -f5- ${FN}.pq.tab) >  ${FN}.afs.dp.pq.tab
    python ${BASEDIR}/Scripts/get_positions_from_sync.py -a ${FN}.afs.dp.pq.tab -b ${FN}_snpeff.tab --both | cat <(head -1  ${FN}.afs.dp.pq.tab )  -  >  ${FN}.afs.anno.tab ;
done
for i in 0.1 0.05 0.001 ;
do
    LFVCF=lofreq2_all_samp_bed_norm_${i}.vcf
    FN=`basename $LFVCF .vcf`
    LOGFILE=${FN}.log
    ERRORLOG=${FN}.err.log

    sed ' s/DQ361065\.7/NP/g; s/DQ361065\.4/GP/g; s/DQ361066\.4/geneZ/g ; s/DQ361066\.7/geneL/g; s/\[[0-9]*\]//g' <  ${FN}.afs.anno.tab >  ${FN}.afs.anno_alt.tab ;
    
done

for i in 0.1 0.05 0.001 ;
do
    LFVCF=lofreq2_all_samp_bed_norm_${i}.vcf
    FN=`basename $LFVCF .vcf`
    LOGFILE=${FN}.log
    ERRORLOG=${FN}.err.log

    sed 's/DQ361065\.7/NP/g; s/DQ361065\.4/GP/g; s/DQ361066\.4/geneZ/g ; s/DQ361066\.7/geneL/g; s/\[[0-9]*\]//g' <  ${FN}.stats.tab >  ${FN}.stats_alt.tab ;

done

