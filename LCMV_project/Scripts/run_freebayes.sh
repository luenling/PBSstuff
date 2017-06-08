#!/bin/bash
#----------
# author: Lukas Endler
# date: 20.9.2015 at 16:46
# takes a bam file calls variants with freebayes
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses_short.fasta
PICARD=/usr/local/Cellar/picard-tools/2.5.0/share/java/picard.jar
GATK=/Volumes/Temp/Lukas/LCMV_project/Tools/GenomeAnalysisTK-3.4-46.jar
SAMTOOLS=/usr/local/bin/samtools
#FB=$BASEDIR/Tools/freebayes/bin/freebayes
FB=/usr/local/bin/freebayes
FN=`basename $1 .list`
MAXCOV=1000000
LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log
INF=$1

if [ $2 ] ; then
    if [[ $2 == "clip" ]] ; then
	echo clipping overlaps at  `date` >> $LOGFILE
	while read BF ;
	      do
		  echo bam clipOverlap --in $BF --out `basename $BF .bam`_co.bam --stats >> $LOGFILE
		  bam clipOverlap --in $BF --out  `basename $BF .bam`_co.bam --stats >> $LOGFILE 2>> $ERRORLOG
		  samtools index `basename $BF .bam`_co.bam
		  echo `basename $BF .bam`_co.bam >> ${FN}_co.list
	done < $1
	INF=${FN}_co.list
    else	
	MAXCOV=$2
    fi
fi
echo "start freebayes at" `date` >> $LOGFILE
echo $FB --fasta-reference $REFGENOME --hwe-priors-off --pooled-continuous --allele-balance-priors-off \
     -u --haplotype-length 0 --use-best-n-alleles 4 --bam-list $INF --min-coverage 50\
        --min-alternate-count 3  --min-alternate-fraction 0.001  -Q 10  --read-max-mismatch-fraction 0.10 \
        --read-snp-limit 4 --read-indel-limit 1 --min-mapping-quality 30 --min-base-quality 30 \> ${FN}_freebayes_nocmpl.vcf >> $LOGFILE
$FB --fasta-reference $REFGENOME --hwe-priors-off --pooled-continuous --allele-balance-priors-off \
     -u --haplotype-length 0 --use-best-n-alleles 4 --bam-list $INF --min-coverage 50\
        --min-alternate-count 3  --min-alternate-fraction 0.001  -Q 10  --read-max-mismatch-fraction 0.10 \
        --read-snp-limit 4 --read-indel-limit 1 --min-mapping-quality 30 --min-base-quality 30  > ${FN}_freebayes_nocmpl.vcf 2>> $ERRORLOG
echo "finished with freebayes at " `date` >> $LOGFILE 



exit 0

for i in 0.01 0.05 0.1; do
    python ~/LCMV_project/Scripts/add_afs_freebayes.py -i ${FN}_freebayes_nocmpl.vcf -m $i > ${FN}_freebayes_nocmpl_$i.vcf
    bash ~/LCMV_project/Scripts/call_snpeff_freebayse.sh ${FN}_freebayes_nocmpl_$i.vcf
done

input filters:

   -4 --use-duplicate-reads
                   Include duplicate-marked alignments in the analysis.
                   default: exclude duplicates marked as such in alignments
   -m 30
                   Exclude alignments from analysis if they have a mapping
                   quality less than Q.  default: 1
   -q 30
                   Exclude alleles from analysis if their supporting base
                   quality is less than Q.  default: 0
   -R --min-supporting-allele-qsum Q
                   Consider any allele in which the sum of qualities of supporting
                   observations is at least Q.  default: 0
   -Y --min-supporting-mapping-qsum Q
                   Consider any allele in which and the sum of mapping qualities of
                   supporting reads is at least Q.  default: 0
   -Q 10
                   Count mismatches toward --read-mismatch-limit if the base
                   quality of the mismatch is >= Q.  default: 10
   -U 
                   Exclude reads with more than N mismatches where each mismatch
                   has base quality >= mismatch-base-quality-threshold.
                   default: ~unbounded
   -z --read-max-mismatch-fraction 0.10
                   Exclude reads with more than N [0,1] fraction of mismatches where
                   each mismatch has base quality >= mismatch-base-quality-threshold
                   default: 1.0
--read-snp-limit 4
                   Exclude reads with more than N base mismatches, ignoring gaps
                   with quality >= mismatch-base-quality-threshold.
                   default: ~unbounded
   -e --read-indel-limit 1
                   Exclude reads with more than N separate gaps.
                   default: ~unbounded
   -0 --standard-filters  Use stringent input base and mapping quality filters
                   Equivalent to -m 30 -q 20 -R 0 -S 0
   -F --min-alternate-fraction 0.005
                   Require at least this fraction of observations supporting
                   an alternate allele within a single individual in the
                   in order to evaluate the position.  default: 0.2
   -C --min-alternate-count 3
                   Require at least this count of observations supporting
                   an alternate allele within a single individual in order
                   to evaluate the position.  default: 2
   -3 --min-alternate-qsum N
                   Require at least this sum of quality of observations supporting
                   an alternate allele within a single individual in order
                   to evaluate the position.  default: 0
   -G --min-alternate-total N
                   Require at least this count of observations supporting
                   an alternate allele within the total population in order
                   to use the allele in analysis.  default: 1
   --min-coverage 50
                   Require at least this coverage to process a site. default: 0
   --max-coverage N
                   Do not process sites with greater than this coverage. default: no limit

	
