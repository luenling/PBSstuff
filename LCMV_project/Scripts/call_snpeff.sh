#!/bin/bash
#----------
# author: Lukas Endler
# date: 23.10.2015 at 16:46
# takes a vcf file adds snpeff annotations for lcmv and creates output tables
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses_short.fasta
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
GATK=/Volumes/Temp/Lukas/LCMV_project/Tools/GenomeAnalysisTK-3.4-46.jar
SAMTOOLS=/usr/local/bin/samtools
VARSCAN=$BASEDIR/Tools/varscan/VarScan.v2.4.0.jar
FN=`basename $1 .vcf`
LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log

echo startingrunning snpeff at `date` >> $LOGFILE
echo cut -f 1-5 $1 \| snpeff lcmv -no-intergenic -no "INTRAGENIC" -no-downstream -no-upstream -stats ${FN}_snpeff_log.html - \> ${FN}_snpeff.vcf >> $LOGFILE
cut -f 1-5 $1 | snpeff lcmv -no-intergenic  -no "INTRAGENIC" -no-downstream -no-upstream -stats ${FN}_snpeff_log.html - > ${FN}_snpeff.vcf

ES=$?
echo finished snpeff at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES

java -jar $BASEDIR/Tools/snpEff/SnpSift.jar extractFields ${FN}_snpeff.vcf CHROM POS REF "ANN[0].GENE" "ANN[0].ALLELE" "ANN[0].EFFECT" "ANN[0].AA" "ANN[1].ALLELE" "ANN[1].EFFECT" "ANN[1].AA" > ${FN}_snpeff.tab
bcftools query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%SDP\t%DP\t%RD\t%AD\t%PVAL\t%RBQ\t%ABQ\t%AF]\n"  $1 > ${FN}.stats.tab
bcftools query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%AF]\n"  $1 > ${FN}.afs.tab
bcftools query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%SDP]\n"  $1 > ${FN}.sdps.tab
paste ${FN}.afs.tab <(cut -f5- ${FN}.sdps.tab)  >  ${FN}.afs.sdps.tab
python ${BASEDIR}/Scripts/get_positions_from_sync.py -a ${FN}.afs.sdps.tab -b ${FN}_snpeff.tab --both | cat <(head -1  ${FN}.afs.sdps.tab )  -  >  ${FN}.afs.anno.tab

sed 's/gi\|86440167\|gb\|DQ361066\.1\|/L/g; s/gi\|116563461\|gb\|DQ361065.2\|/S/g; s/DQ361065\.7/NP/g; s/DQ361065\.4/GP/g; s/DQ361066\.4/geneZ/g ; s/DQ361066\.7/geneL/g; s/\[[0-9]*\]//g' <  ${FN}.afs.anno.tab >  ${FN}.afs.anno_alt.tab
sed 's/gi\|86440167\|gb\|DQ361066\.1\|/L/g; s/gi\|116563461\|gb\|DQ361065.2\|/S/g; s/DQ361065\.7/NP/g; s/DQ361065\.4/GP/g; s/DQ361066\.4/geneZ/g ; s/DQ361066\.7/geneL/g; s/\[[0-9]*\]//g' <  ${FN}.stats.tab >  ${FN}.stats_alt.tab
