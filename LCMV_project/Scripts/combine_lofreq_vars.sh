#!/bin/bash
#----------
# author: Lukas Endler
# date: 20.9.2015 at 16:46
# combines vcf files from lofreq and creates a bed for calling with
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses_short.fasta


for i in *_lofreq.vcf; do
    SMP=${i%_lofreq*}
    awk -v OFS="\t" '/^\#\#[^I]/ {print} /^\#\#INFO/ {print $0; sub("INFO","FORMAT",$0); print $0; } /\#CH/ {print $0,"FORMAT","'$SMP'"} !/^\#/ {form=$8;  gsub(/=[^A-Z]+/,":",form); gsub(/;/,":",form); sub(/:$/,"",form); sub(/INDEL:/,"",form);  samp=$8; gsub(/[A-Z4]+=/,"",samp); gsub(/;/,":",samp); sub(/INDEL:/,"",samp); print $0,form,samp}' $i | bgzip -c > ${SMP}_samp.lofreq.vcf.gz
    tabix -p vcf ${SMP}_samp.lofreq.vcf.gz  
done

bcftools merge -O v *_samp.lofreq.vcf.gz > all_samp_bcf.vcf
awk -v OFS="\t" '!/^\#/ {print $1,$2-1,$2}' all_samp_bcf.vcf > all_samp_bcf.bed
