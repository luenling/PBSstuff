FN=`basename $1 .vcf`

bcftools query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%AF][\t%SDP]\n"  $1 | sed 's/gi\|86440167\|gb\|DQ361066\.1\|/L/g; s/gi\|116563461\|gb\|DQ361065.2\|/S/g; s/DQ361065\.7/NP/g; s/DQ361065\.4/GP/g; s/DQ361066\.4/geneZ/g ; s/DQ361066\.7/geneL/g; s/\[[0-9]*\]//g' >  ${FN}.afs.tab
bcftools query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%SDP\t%DP\t%RD\t%AD\t%PVAL\t%RBQ\t%ABQ\t%AF]\n"  $1 | sed 's/gi\|86440167\|gb\|DQ361066\.1\|/L/g; s/gi\|116563461\|gb\|DQ361065.2\|/S/g; s/DQ361065\.7/NP/g; s/DQ361065\.4/GP/g; s/DQ361066\.4/geneZ/g ; s/DQ361066\.7/geneL/g; s/\[[0-9]*\]//g' >  ${FN}.stats_alt.tab
