#!/bin/bash
#----------
# author: Lukas Endler
# date: 17.9.2015 at 16:34
# runs vphaser2 a file with a list of bam files and an optional vcf file from samtools for adding depths and creates vcf files
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses_short.fasta
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
VPHASER2=$BASEDIR/Tools/viral-ngs/intrahost_alt2.py
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$BASEDIR/Tools/VPhaser-2-02112013/bamtools-2.3.0/lib

while read FN; do
	FNB=`basename $FN .bam`
	SAMP=${FNB/_S19*/}
	#picard AddOrReplaceReadGroups I=$FN O=`basename $FN .bam`"_vp.bam" RGID=${SAMP} RGLB=$SAMP RGPL=illumina RGPU=unitx RGSM=$SAMP 2> $FN.picard.err.log
	#FN=`basename $FN .bam`"_vp.bam"
	#samtools index $FN
	SAMPLES=${SAMPLES}${SAMP}" "
	#ALGS=${ALGS}$FN" "
	ISNVS=${ISNVS}${FNB}.tab" "
	LOGFILE=${FNB}.vphaser2.log
	ERRORLOG=${FNB}.vphaser2.err.log
	OMP_NUM_THREADS=8
	echo "start VPHASER2 for " $FN " at" `date` >> $LOGFILE
	echo python $VPHASER2 vphaser_one_sample  --maxBias=100000  --vphaserNumThreads 18 --minReadsEach 2 $FN $REFGENOME ${FNB}.tab >> $LOGFILE
	python $VPHASER2 vphaser_one_sample   --maxBias=100000 --vphaserNumThreads 18 --minReadsEach 2 $FN $REFGENOME ${FNB}.tab 2>> $ERRORLOG >> $LOGFILE
	ES=$?
	echo finished vphaser2 at `date` with exit state $ES >> $LOGFILE
done <$1

samtools faidx $REFGENOME L > L.fasta
samtools faidx $REFGENOME S > S.fasta

for i in $SAMPLES; do
  samtools faidx $REFGENOME L | sed 's/\>L/\>'$i'/g' >> L.fasta
  samtools faidx $REFGENOME S | sed 's/\>S/\>'$i'/g' >> S.fasta    
done

LOGFILE=all.log

echo " merge to vcf " $1 " at" `date` >> all.log
echo python $VPHASER2 merge_to_vcf $REFGENOME all_out.vcf --samples $SAMPLES --isnvs $ISNVS --alignments L.fasta S.fasta  >> all.log
python $VPHASER2 merge_to_vcf $REFGENOME all_out.vcf --samples $SAMPLES --isnvs $ISNVS --alignments L.fasta S.fasta  2>> all.err.log >> all.log

bgzip -c  all_out.vcf >  all_out.vcf.gz
tabix -p vcf  all_out.vcf.gz

DPS=$2
LVCF=all_out
if [[ -e $DPS ]];
then
    ~/Tools/bcftools-1.3.1/bcftools annotate -a $DPS -c "+FORMAT/DP" all_out.vcf.gz > all_out_dps.vcf
    LVCF=all_out_dps
fi

for i in 0.1 0.05 0.01 ;
do
    echo bcftools norm  -f $REFGENOME -m -any ${LVCF}.vcf \| bcftools view -e 'max(AF[*])<'$i' ' - \| bcftools norm  -f $REFGENOME -m +any - \> ${LVCF}_${i}.vcf >> $LOGFILE
    bcftools norm  -f $REFGENOME -m - ${LVCF}.vcf | bcftools view -e 'max(AF[*])<'$i' ' - | bcftools norm  -f $REFGENOME -m +any - > ${LVCF}_${i}.vcf
    echo starting running snpeff at `date` >> $LOGFILE
    echo snpeff lcmv -no-intergenic -no "INTRAGENIC" -no-downstream -no-upstream -stats ${LVCF}_${i}_snpeff_log.html ${LVCF}_${i}.vcf \> ${LVCF}_${i}_anno.vcf >> $LOGFILE
    snpeff lcmv -no-intergenic -no "INTRAGENIC" -no-downstream -no-upstream -stats ${LVCF}_${i}_snpeff_log.html ${LVCF}_${i}.vcf > ${LVCF}_${i}_anno.vcf
    java -jar $BASEDIR/Tools/snpEff/SnpSift.jar extractFields ${LVCF}_${i}_anno.vcf CHROM POS REF "ANN[0].GENE" "ANN[0].ALLELE" "ANN[0].EFFECT" "ANN[0].AA" "ANN[1].ALLELE" "ANN[1].EFFECT" "ANN[1].AA" > ${LVCF}_${i}_anno.tab
    bcftools query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%LB\t%DP\t%AF]\n" ${LVCF}_${i}_anno.vcf > ${LVCF}_${i}.stats.tab
    bcftools query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%AF]\n"  ${LVCF}_${i}_anno.vcf > ${LVCF}_${i}.afs.tab
    bcftools query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%DP]\n"  ${LVCF}_${i}_anno.vcf > ${LVCF}_${i}.dp.tab
    paste ${LVCF}_${i}.afs.tab <(cut -f5- ${LVCF}_${i}.dp.tab)  >  ${LVCF}_${i}.afs.dp.tab
    python ${BASEDIR}/Scripts/get_positions_from_sync.py -a ${LVCF}_${i}.afs.dp.tab -b ${LVCF}_${i}_anno.tab --both | cat <(head -1  ${LVCF}_${i}.afs.dp.tab )  -  >  ${LVCF}_${i}.afs.anno.tab ;
done


for i in 0.1 0.05 0.01 ;
do    

    sed 's/DQ361065\.7/NP/g; s/DQ361065\.4/GP/g; s/DQ361066\.4/geneZ/g ; s/DQ361066\.7/geneL/g; s/\[[0-9]*\]//g' <  ${LVCF}_${i}.afs.anno.tab >  ${LVCF}_${i}.afs.anno_alt.tab ;

done    

for i in 0.1 0.05 0.01 ;
do    

    sed 's/DQ361065\.7/NP/g; s/DQ361065\.4/GP/g; s/DQ361066\.4/geneZ/g ; s/DQ361066\.7/geneL/g; s/\[[0-9]*\]//g' <  ${LVCF}_${i}.stats.tab >  ${LVCF}_${i}.stats_alt.tab ;

done

exit 0


if [ ! -d VPHASER2/$FN ]; then
    echo mkdir >> $LOGFILE
    mkdir -p VPHASER2/$FN
fi





/Volumes/Temp/Lukas/LCMV_project/Tools/viral-ngs/intrahost.py merge_to_vcf  $REFGENOMNE --samples "ND_0,ND_ARM_infected_BHK21_S11759" --isnvs "ND_0_noprime_sorted_viral_sh_max_cov_10000_short_RG.tab,ND_ARM_infected_BHK21_S11759_noprime_sorted_viral_sh_max_cov_10000_short_RG.tab"  --alignments "../ND_0_noprime_sorted_viral_sh_max_cov_10000_short_RG.bam,../ND_ARM_infected_BHK21_S11759_noprime_sorted_viral_sh_max_cov_10000_short_RG.bam"


