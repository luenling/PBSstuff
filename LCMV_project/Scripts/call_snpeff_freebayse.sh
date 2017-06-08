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

echo starting running snpeff at `date` >> $LOGFILE
echo cut -f 1-5 $1 \| snpeff lcmv -no-intergenic -no "INTRAGENIC" -no-downstream -no-upstream -stats ${FN}_snpeff_log.html - \> ${FN}_snpeff.vcf >> $LOGFILE
cut -f 1-5 $1 | snpeff lcmv -no-intergenic -no-downstream -no "INTRAGENIC" -no-upstream -stats ${FN}_snpeff_log.html - > ${FN}_snpeff.vcf

ES=$?
echo finished snpeff at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES

java -jar $BASEDIR/Tools/snpEff/SnpSift.jar extractFields ${FN}_snpeff.vcf CHROM POS REF "ANN[0].GENE" "ANN[0].ALLELE" "ANN[0].EFFECT" "ANN[0].AA" "ANN[1].ALLELE" "ANN[1].EFFECT" "ANN[1].AA" "ANN[2].ALLELE" "ANN[2].EFFECT" "ANN[2].AA" "ANN[3].ALLELE" "ANN[3].EFFECT" "ANN[3].AA" "ANN[4].ALLELE" "ANN[4].EFFECT" "ANN[4].AA" > ${FN}_snpeff.tab
bcftools query -Hf "%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%MQMR\t%MQM\t%FSB\t%SB[\t%DP\t%RO\t%AO\t%AF]\n"  $1 > ${FN}.stats.tab
bcftools query -Hf "%CHROM\t%POS\t%REF\t%ALT\t%TYPE[\t%DP][\t%AF]\n"  $1 > ${FN}.afs.tab
python ${BASEDIR}/Scripts/get_positions_from_sync.py -a ${FN}.afs.tab -b ${FN}_snpeff.tab --both | cat <(head -1  ${FN}.afs.tab )  -  >  ${FN}.afs.anno.tab

sed 's/_S194[^:]*//g; s/_trimmed//g; s/gi\|86440167\|gb\|DQ361066\.1\|/L/g; s/gi\|116563461\|gb\|DQ361065.2\|/S/g; s/DQ361065\.7/NP/g; s/DQ361065\.4/GP/g; s/DQ361066\.4/geneZ/g ; s/DQ361066\.7/geneL/g; s/\[[0-9]*\]//g' <  ${FN}.afs.anno.tab >  ${FN}.afs.anno_alt.tab
sed 's/_S194[^:]*//g; s/_trimmed//g; s/gi\|86440167\|gb\|DQ361066\.1\|/L/g; s/gi\|116563461\|gb\|DQ361065.2\|/S/g; s/DQ361065\.7/NP/g; s/DQ361065\.4/GP/g; s/DQ361066\.4/geneZ/g ; s/DQ361066\.7/geneL/g; s/\[[0-9]*\]//g' <  ${FN}.stats.tab >  ${FN}.stats_alt.tab
