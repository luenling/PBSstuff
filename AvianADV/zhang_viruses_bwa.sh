#!/bin/bash
#----------
# author: Lukas Endler
# Time-stamp: <2016-07-07 14:30:45 lukasendler>
# date: 03.03.2016 at 20:31
# takes two fastq files runs bwa with 12 threads and outputs a bam file called like the fastq prefix against gallus gallus
# only keeps the non mapped reads and puts them into fastq files
# trims them
# call with bash command blub_1.fq blub_2.fq >> logfile.log 2>> log.error.log
#--------------


BASEDIR=/Volumes/Temp/Anna.Schachner
REFPHI=$BASEDIR/References/phix.fasta
REFGENOME=$BASEDIR/Zhang_new/Hep_Viruses/avian_HPV_phiX.fa
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
BBDUK=/Volumes/Temp/Lukas/LCMV_project/Tools/bbmap/bbduk.sh
BBMERGE=/Volumes/Temp/Lukas/LCMV_project/Tools/bbmap/bbmerge.sh
BWA=/usr/local/Cellar/bwa/0.7.13/bin/bwa
ADAPTS=/Volumes/Temp/Lukas/LCMV_project/Tools/bbmap/resources/adapters.fa

FN=`basename $1 .gz`
FN=`basename $FN _nogal_12.fq`
FN=`basename $FN _1.fq`
FN=`basename $FN _1.fastq`
# just to replace the whole bloody name tag
R2=`echo $1 | sed 's/_1\.f/_2.f/'`
FN=${FN/*_H325VBBXX_5_/}
SAMPLE=$FN
LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log

RG='@RG\tID:'$FN'\tSM:'$FN
echo "start bwa mem to phiX  at" `date` >> $LOGFILE
echo $BWA mem -R $RG -M -t 5 $REFPHI $1  $R2 \| samtools view -Shb -F 256 - \| $SAMTOOLS sort -T ${FN}_temp - \> ${FN}_phi.bam  >> $LOGFILE

$BWA mem -R $RG -M -t 5 $REFPHI $1 $R2 2>> $ERRORLOG | samtools view -Shb -F 256 - | $SAMTOOLS sort -T ${FN}_temp - > ${FN}_phi.bam
ES=$?

echo finished bwa mem mapping to phix at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
echo "remove phiX at " `date` >> $LOGFILE
echo $SAMTOOLS view -bf 13  -F 256 ${FN}_phi".bam" \|  $SAMTOOLS bam2fq - \| gzip -c - \> ${FN}_nophiX_12.fq.gz >> $LOGFILE
$SAMTOOLS view -b -f 13  -F 256 ${FN}_phi".bam" |  $SAMTOOLS bam2fq - | gzip -c - > ${FN}_nophiX_12.fq.gz
echo $SAMTOOLS index  ${FN}_phi.bam >> $LOGFILE
$SAMTOOLS index  ${FN}_phi.bam
echo $SAMTOOLS flagstat  ${FN}_phi.bam >> $LOGFILE
$SAMTOOLS flagstat  ${FN}_phi.bam >> $LOGFILE
echo $SAMTOOLS idxstats  ${FN}_phi.bam >> $LOGFILE
$SAMTOOLS idxstats  ${FN}_phi.bam >> $LOGFILE


echo trimming and removing adaptors at `date` >> ${LOGFILE}
echo $BBDUK -Xmx4g in1=${FN}_nophiX_12.fq.gz out1=${FN}_trim_1.fq.gz out2=${FN}_trim_2.fq.gz qtrim=r minlen=100 ktrim=r ref=\"$ADAPTS\" k=20 mink=15 trimq=10 overwrite=t tbo=t tpe=t lhist=${FN}_rlen.hist stats=${FN}_trim.stats \>\> ${LOGFILE} 2\>\> ${ERRORLOG} >> ${LOGFILE}
 $BBDUK -Xmx4g in1=${FN}_nophiX_12.fq.gz out1=${FN}_trim_1.fq.gz out2=${FN}_trim_2.fq.gz qtrim=r minlen=100 ktrim=r ref=\"$ADAPTS\" k=20 mink=15 trimq=10 overwrite=t tbo=t tpe=t lhist=${FN}_rlen.hist stats=${FN}_trim.stats \>\> ${LOGFILE} 2\>\> ${ERRORLOG}  >> ${LOGFILE} 2>> ${ERRORLOG}
ES=$?
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES

fastq-stats ${FN}_trim_2.fq.gz > ${FN}_trim_2.fq.stats &
fastq-stats ${FN}_trim_1.fq.gz > ${FN}_trim_1.fq.stats &

echo merging reads with flash at `date` >>  ${LOGFILE}
echo flash -t 10 -d Flash --max-overlap 300 -p 33 -z -o ${FN}_trim_flash -x 0.25 ${FN}_trim_1.fq.gz ${FN}_trim_2.fq.gz \>\> ${FN}_sub_flash.log >>  ${LOGFILE}

flash -t 10 -d Flash --max-overlap 300 -p 33 -z -o ${FN}_trim_flash -x 0.25 ${FN}_trim_1.fq.gz ${FN}_trim_2.fq.gz >> ${FN}_sub_flash.log

READNUM=100000
SEED=$RANDOM
echo sample fastq reads at `date` >> ${LOGFILE}
fastq-sample -s $SEED -o  Flash/${FN}_trim_flash_sub -n $READNUM Flash/${FN}_trim_flash.extendedFrags.fastq.gz


exit 0

fastq-stats $1 > ${1}.stats &
fastq-stats $R2 > ${R2}.stats &

RG='@RG\tID:'$FN'\tSM:'$FN
echo "start bwa mem to HepV and phiX  at" `date` >> $LOGFILE
echo $BWA mem -R $RG -M -t 5 $REFGENOME $1  $R2 \| samtools view -Shb -F 256 - \| $SAMTOOLS sort -T ${FN}_temp - \> ${FN}_hpv_phi.bam  >> $LOGFILE

$BWA mem -R $RG -M -t 5 $REFGENOME $1  $R2 2>> $ERRORLOG |samtools view -Shb -F 256 - | $SAMTOOLS sort -T ${FN}_temp - > ${FN}_hpv_phi.bam
ES=$?
echo finished bwa mem mapping to phix at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES

echo $SAMTOOLS index  ${FN}_hpv_phi.bam >> $LOGFILE
$SAMTOOLS index  ${FN}_hpv_phi.bam
echo $SAMTOOLS flagstat  ${FN}_hpv_phi.bam >> $LOGFILE
$SAMTOOLS flagstat  ${FN}_hpv_phi.bam >> $LOGFILE
echo $SAMTOOLS idxstats  ${FN}_hpv_phi.bam >> $LOGFILE
$SAMTOOLS idxstats  ${FN}_hpv_phi.bam >> $LOGFILE

exit 0


BASEDIR=/Volumes/Temp/Anna.Schachner
REFPHI=$BASEDIR/References/phix.fasta
REFGENOME=$BASEDIR/Samples/Zhang/Hep_Viruses/avian_hep_viruses.fa
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
BBDUK=/Volumes/Temp/Lukas/LCMV_project/Tools/bbmap/bbduk.sh
BBMERGE=/Volumes/Temp/Lukas/LCMV_project/Tools/bbmap/bbmerge.sh
BWA=/usr/local/Cellar/bwa/0.7.12/bin/bwa
ADAPTS=/Volumes/Temp/Lukas/LCMV_project/Tools/bbmap/resources/adapters.fa


FN=`basename $1 .gz`
FN=`basename $FN _nogal_12.fq`
FN=`basename $FN _1.fq`
FN=`basename $FN _1.fastq`
# just to replace the whole bloody name tag
R2=`echo $1 | sed 's/_1\.f/_2.f/'`
FN=${FN/*_H325VBBXX_5_/}
SAMPLE=$FN
LOGFILE=${FN}.log
ERRORLOG=${FN}.err.log

RG='@RG\tID:'$FN'\tSM:'$FN
echo "start bwa mem to phiX  at" `date` >> $LOGFILE
echo $BWA mem -R $RG -M -t 5 $REFPHI $1  $R2 \| samtools view -Shb -F 256 - \| $SAMTOOLS sort -T ${FN}_temp - \> ${FN}_phi.bam  >> $LOGFILE

$BWA mem -R $RG -M -t 5 $REFPHI $1 $R2 2>> $ERRORLOG |samtools view -Shb -F 256 - | $SAMTOOLS sort -T ${FN}_temp - > ${FN}_phi.bam
ES=$?
echo finished bwa mem mapping to phix at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
echo "remove phiX at " `date` >> $LOGFILE
echo $SAMTOOLS view -bf 13  -F 256 ${FN}_phi".bam" \|  $SAMTOOLS bam2fq - \| gzip -c - \> ${FN}_nophiX_12.fq.gz >> $LOGFILE
$SAMTOOLS view -b -f 13  -F 256 ${FN}_phi".bam" |  $SAMTOOLS bam2fq - | gzip -c - > ${FN}_nophiX_12.fq.gz
echo $SAMTOOLS index  ${FN}_phi.bam >> $LOGFILE
$SAMTOOLS index  ${FN}_phi.bam
echo $SAMTOOLS flagstat  ${FN}_phi.bam >> $LOGFILE
$SAMTOOLS flagstat  ${FN}_phi.bam >> $LOGFILE
echo $SAMTOOLS idxstats  ${FN}_phi.bam >> $LOGFILE
$SAMTOOLS idxstats  ${FN}_phi.bam >> $LOGFILE

fastq-stats $1 > ${1}.stats &
fastq-stats $R2 > ${R2}.stats &


echo trimming and removing adaptors at `date` >> ${LOGFILE}
echo $BBDUK -Xmx4g in1=${FN}_nophiX_12.fq.gz out1=${FN}_trim_1.fq out2=${FN}_trim_2.fq qtrim=r minlen=50 ktrim=r ref=\"$ADAPTS\" k=20 mink=15 trimq=15 overwrite=t tbo=t tpe=t lhist=${FN}_rlen.hist stats=${FN}_trim.stats \>\> ${LOGFILE} 2\>\> ${ERRORLOG} >> ${LOGFILE}
$BBDUK -Xmx4g in1=${FN}_nophiX_12.fq.gz out1=${FN}_trim_1.fq out2=${FN}_trim_2.fq qtrim=r minlen=50 ktrim=r ref=\"$ADAPTS\" k=20 mink=15 trimq=15 overwrite=t tbo=t tpe=t lhist=${FN}_rlen.hist stats=${FN}_trim.stats  >> ${LOGFILE} 2>> ${ERRORLOG}
ES=$?
echo finished at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES

fastq-stats ${FN}_trim_2.fq > ${FN}_trim_2.fq.stats &
fastq-stats ${FN}_trim_1.fq > ${FN}_trim_1.fq.stats &


RG='@RG\tID:'$FN'\tSM:'$FN
echo "start bwa mem to hep  at" `date` >> $LOGFILE
echo $BWA mem -R $RG -M -t 5 $REFGENOME ${FN}_trim_1.fq  ${FN}_trim_2.fq \| samtools view -Shb -F 256 - \| $SAMTOOLS sort -T ${FN}_temp - \> ${FN}_hep.bam  >> $LOGFILE
$BWA mem -R $RG -M  -t 5 $REFGENOME ${FN}_trim_1.fq  ${FN}_trim_2.fq 2>> $ERRORLOG | samtools view -Shb -F 256 - | $SAMTOOLS sort -T ${FN}_temp - > ${FN}_hep.bam
ES=$?
echo finished bwa mem mapping to hep at `date` with exit state $ES >> $LOGFILE
[ $ES -eq 0 ] || exit $ES
echo "check hep reads at " `date` >> $LOGFILE
echo $SAMTOOLS index  ${FN}_hep.bam >> $LOGFILE
$SAMTOOLS index  ${FN}_hep.bam
echo $SAMTOOLS flagstat  ${FN}_hep.bam >> $LOGFILE
$SAMTOOLS flagstat  ${FN}_hep.bam >> $LOGFILE
echo $SAMTOOLS idxstats  ${FN}_hep.bam >> $LOGFILE
$SAMTOOLS idxstats  ${FN}_hep.bam >> $LOGFILE

exit

blastn  -db ../Hep_Viruses/avian_hep_viruses -query ${FN} -evalue 1e-4 -num_threads 5 -outfmt '7 qseqid sseqid evalue qlen qstart qend length sstart send slen pident stitle qcovs' -out ${FN}_blast.out
or FN in sample_id_3557*.fa; do blastn  -db ../Hep_Viruses/avian_hep_viruses -query ${FN} -evalue 1e-2 -num_threads 1 -outfmt '7 qseqid sseqid evalue qlen qstart qend length sstart send slen pident stitle qcovs' -out ${FN}_blast.out; done
makeblastdb -in $REFPHI -input_type fasta -dbtype nucl -parse_seqids -out phix_db
for FN in sample_id_3557*.gz.fa; do blastn  -db phix_db -query ${FN} -evalue 1e-2 -num_threads 1 -outfmt '7 qseqid sseqid evalue qlen qstart qend length sstart send slen pident stitle qcovs' -out ${FN}_phi_blast.out; done
for FN in sample_id_3557*.gz.fa; do blastn  -db phix_db -query ${FN} -evalue 1e-2 -num_threads 1 -outfmt '7 qseqid sseqid evalue qlen qstart qend length sstart send slen pident stitle qcovs' -out ${FN}_phi_blast.out; done


makeblastdb -in hep_phix.fna -input_type fasta -dbtype nucl -parse_seqids -out hep_phix
for FN in sample_id_3557*.gz.fa; do blastn  -db hep_phix -query ${FN} -evalue 1e-1 -num_threads 10 -outfmt '7 qseqid sseqid evalue qlen qstart qend length sstart send slen pident stitle qcovs' -out ${FN}_hepphix_blast.out; done

wget -O avian_hep_viruses.fa "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AM943647,JN597006,AY535004,EF206691,AM943646,GU954430,JN997392,KF511797&rettype=fasta&retmode=text"
wget -O avian_ibv.fa "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001451.1,M95169.1&rettype=fasta&retmode=text"

wget -O avian_ibv.fa "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001451.1,M95169.1&rettype=fasta&retmode=text"
cat avian_ibv.fa hep_phix.fna > ibv_hep.phix.fna
makeblastdb -in ibv_hep.phix.fna -input_type fasta -dbtype nucl -parse_seqids -out ibv_hep.phix
for FN in sample_id_3557*.gz.fa; do blastn  -db ibv_hep.phix -query ${FN} -evalue 1e-1 -num_threads 5 -outfmt '7 qseqid sseqid evalue qlen qstart qend length sstart send slen pident stitle qcovs' -out ${FN}_ibvhepphix_blast.out; done
#only phiX
grep  "^[^#]" sample_id_3557*.fq.gz.fa_ibvhepphix_blast.out | grep -v "phiX174" |  less

export BLASTDB=/Volumes/Temp/BlastDBs;  { for FN in sample_id_3557*_trim_[12].fa; do blastn  -db nt -query ${FN} -evalue 1e-2 -num_threads 10 -outfmt '7 qseqid sseqid evalue qlen qstart qend length sstart send slen pident stitle qcovs  sscinames scomnames sskingdoms' -num_alignments 1 -out ${FN}_nt_blast.out; done } &

# get common names and count them
for FN in *1.fa_nt_blast.out; do
    FBN=${FN/_trim_[12].fa_nt_blast.out/}
    echo reads : `wc -l ${FBN}_trim_1.fa ${FBN}_trim_2.fa` >  ${FBN}.species
    echo file ${FBN}_trim_1.fa hits found: `grep -c "\# [^0].* hits found" ${FBN}_trim_1.fa_nt_blast.out` >>  ${FBN}.species
    grep -A 1 "\# [^0].* hits found" $FN | grep -v "^[\#-]" | cut -f 15 | sort | uniq -c | sort -rn -k1 >> ${FBN}.species
    echo file ${FBN}_trim_2.fa hits found: `grep -c "\# [^0].* hits found" ${FBN}_trim_2.fa_nt_blast.out` >>  ${FBN}.species
    grep -A 1 "\# [^0].* hits found" ${FBN}_trim_2.fa_nt_blast.out | grep -v "^[\#-]" | cut -f 15 | sort | uniq -c  | sort -nr -k1 >> ${FBN}.species
done
    
wget ftp://ftp.ncbi.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
wget ftp://ftp.ncbi.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz


#kraken
for i in ../Flash/new_sample_id_*_trim_flash.extendedFrags.fastq.gz; do
    BN=`basename $i _trim_flash.extendedFrags.fastq.gz`
    kraken --db minikraken_20141208 --fastq-input --gzip-compressed $i --output ${BN}_kraken.out --preload --threads 10 2>> kraken.log
    #kraken-report --db minikraken_20141208 --show-zeros ${BN}_kraken.out > ${BN}.full_report
    kraken-report --db minikraken_20141208 ${BN}_kraken.out > ${BN}.report
done

for i in ../*_1.fq.gz; do
    BN=`basename $i _1.fq.gz`
    kraken --db minikraken_20141208 --fastq-input --paired --gzip-compressed $i ../${BN}_2.fq.gz --output ${BN}_PE_kraken.out --preload --threads 10 2>> kraken_PE.log
    #kraken-report --db minikraken_20141208 --show-zeros ${BN}_kraken.out > ${BN}.full_report
    kraken-report --db minikraken_20141208 ${BN}_PE_kraken.out > ${BN}_PE.report
done


Trinity --seqType fq --max_memory 20G --CPU 10 --single ../Flash/new_sample_id_35573_trim_flash.extendedFrags.fastq.gz --output Trinity_sample_id_35573_trim_flash --min_contig_length 200
# blast trinity results for krona
blastn -task megablast -db /Volumes/Temp/BlastDBs/nt -evalue 1e-3 -num_threads 4 -query Trinity.fasta -outfmt 7 -max_target_seqs 5 > trinity_megablast.tab &

blastn -task megablast -db /Volumes/Temp/BlastDBs/nt -evalue 1e-3 -num_threads 8 -query Trinity.fasta -outfmt 7 -max_target_seqs 5 > trinity_megablast.tab &


SPADES=/Volumes/Temp/Lukas/Tools/SPAdes-3.7.1-Darwin/bin/spades.py
for i in ../new_sample_id_355*_trim_1.fq.gz; do
	 BN=`basename $i _1.fq.gz`
	 $SPADES --meta -k 21,33,55 --pe1-1 $i  --pe1-2 ../${BN}_2.fq.gz -o Spade_${BN} >> ${BN}_spades_meta.log 2>>  ${BN}_spades_meta.error.log
done


~/Tools/CLARKSCV1.2.3/exe/CLARK -k 19 -D ./DBS_virus -n 10 -P ../new_sample_id_35573_1.fq.gz ../new_sample_id_35573_2.fq.gz -T ~/Tools/CLARKSCV1.2.3/DBS_virus/targets.txt -R new_sample_id_35573.k19.results

 ~/LCMV_project/Tools/bbmap/bbduk.sh in=new_sample_id_35573_1.fq.gz in2=new_sample_id_35573_2.fq.gz literal="AATGTGCTGCGGGGTGTCAA,CATCTGGTACCGTGCGAGTA" out=unmatched.fq outm=matched.fq k=15 hdist=1 stats=stats.txt
 fastq_to_fasta -i matched.fq -o matched_primers_sample_id_35573.fna
 blastn -task blastn -db /Volumes/Temp/BlastDBs/nt -evalue 1e-3 -num_threads 8 -query matched_primers_sample_id_35573.fna -outfmt 7 -max_target_seqs 5 > blastn.tab &


for i in ../new_sample_id_355*_trim_1.fq.gz; do
    BN=`basename $i _1.fq.gz`
    ~/LCMV_project/Tools/bbmap/bbduk.sh -Xmx10g in=$i in2=../${BN}_2.fq.gz k=20 hdist=2 mkh=3 ref=~/Anna.Schachner/Zhang_new/Hep_Viruses/avian_hep_E_ibv.fa out=/dev/null minrskip=5 outm=${BN}_matched.fasta stats=${BN}_stats.txt overwrite=true 2> ${BN}.log
done

blastn -task blastn -db ~/Anna.Schachner/Zhang_new/Hep_Viruses/avian_hep_E_ibv -evalue 1e-2 -num_threads 12 -query ${BN}_matched.fasta -outfmt 7 -num_alignments 2 >  ${BN}_matched.blastn.tab &

blat -t=dna -q=dna -tileSize=10 -minMatch=1  ~/Anna.Schachner/Zhang_new/Hep_Viruses/avian_hep_E_ibv.fa ${BN}_matched.fasta  ${BN}_matched.blat.psl

DATABASE     avianHVE    ~/Anna.Schachner/Zhang_new/FastQ_screen/avian_hep_viruses
DATABASE     avianIBV    ~/Anna.Schachner/Zhang_new/FastQ_screen/avian_ibv


for i in ../new_sample_id_355*_trim_1.fq.gz; do
    BN=`basename $i _1.fq.gz`
    ~/LCMV_project/Tools/bbmap/bbsplit.sh -Xmx10g in=$i in2=../${BN}_2.fq.gz ref=avian_ibv.fa,avian_hep_viruses.fa,gallus_gallus_rna.fa.gz basename=${BN}_o%.fq.gz slow=t >> ${BN}.log 2>&1
done


~/Tools/CLARKSCV1.2.3/exe/CLARK -k 19 -D ./DBS_avian_virus -n 10 -P ../new_sample_id_35573_1.fq.gz ../new_sample_id_35573_2.fq.gz -T ./avian.targets.txt -R new_sample_id_35573.k19.results


GALLUSGENOME=/Volumes/Temp/Anna.Schachner/References/Gallus_gallus.Galgal4.dna.toplevel.fa.gz
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools
BBDUK=/Volumes/Temp/Lukas/LCMV_project/Tools/bbmap/bbduk.sh
BBMERGE=/Volumes/Temp/Lukas/LCMV_project/Tools/bbmap/bbmerge.sh
BWA=/usr/local/Cellar/bwa/0.7.13/bin/bwa
FN=new_sample_id_35571
LOGFILE=$FN.log
echo $BWA mem -t 15 $GALLUSGENOME ${FN}_trim_1.fq.gz  ${FN}_trim_2.fq.gz \| samtools view -Shb - \| $SAMTOOLS sort -T ${FN}_temp - \> ${FN}_gallus.bam  >> $LOGFILE
$BWA mem -t 15 $GALLUSGENOME ${FN}_trim_1.fq.gz  ${FN}_trim_2.fq.gz | samtools view -Shb - | $SAMTOOLS sort -T ${FN}_temp - > ${FN}_gallus.bam
samtools flagstat ${FN}_gallus.bam >> $LOGFILE


for i in ../new_sample_id_355*_trim_1.fq.gz; do
    BN=`basename $i _1.fq.gz`
    SMP=${BN#*id_}
    SMP=${SMP/_trim*}
    REFFN=`ls *${SMP}*.fa`
    ~/LCMV_project/Tools/bbmap/bbduk.sh -Xmx10g in=$i in2=../${BN}_2.fq.gz k=25 hdist=2 hdist2=0 ref=$REFFN out=/dev/null minrskip=5 outm=${BN}_matched_1.fq outm2=${BN}_matched_2.fq stats=${BN}_stats.txt overwrite=true 2> ${BN}.log
done
