#!/bin/bash
#----------
# author: Lukas Endler
# authored: 2016-05-02 12:01:03 
# Time-stamp: <2016-09-06 09:54:02 lukasendler>
# command lines
# usage: command DATADIR BAMSUFFIX
# default: DATADIR: /Volumes/Temp/Hammer/STAR_2nd_pass/ BAMSUFFIX: _Aligned.sorted.bam
#--------------

DATA=/Volumes/Temp/Lukas/Hammer/STAR_2nd_pass/
SCRIPTS=/Volumes/Temp/Lukas/Hammer/scripts/
BRIDGER=${HOME}/Tools/Bridger_r2014-12-01/Bridger.pl
TRINITY=${HOME}/Tools/trinityrnaseq-2.2.0/Trinity
BBDUK=${HOME}/Tools/bbmap/bbduk.sh
BBFIX=${HOME}/Tools/bbmap/repair.sh
ADAPTS=${HOME}/Tools/bbmap/resources/adapters.fa
BBMERGE=${HOME}/Tools/bbmap/bbmerge.sh
GTF=/Volumes/Temp/Lukas/Hammer/annotations/Sus_scrofa.Sscrofa10.2.84.with.chr.no_gene.gtf
GTF_DIR=/Volumes/Temp/Hammer/annotations/
REFGENOME=/Volumes/Temp/Lukas/Hammer/annotations/susScr3.fa
REF_DIR=/Volumes/Temp/Lukas/Hammer/annotations/
SUFFIX=_Aligned.sorted.bam
TCRD=chr7:81500000-83000000
TCRG=chr9:119505000-119687000
#TCRD=chr7:81500000-82300000
#TCRG=chr9:119500000-120000000
BINPACKER=${HOME}/Tools/BinPacker_1.1/BinPacker

# BAM=$1
# BN=`basename $1 .bam`
# LOGFILE=${BN}".log"

if [[ -d $1 ]]
then
    DATA=$1
fi


OUTDIR=$DATA/Assemblies_alt/
mkdir -p $OUTDIR

if [ $2 ]
then
    SUFFIX=$2
fi

LOGFILE=${OUTDIR}/"logfile.log"
ERRORFILE=${OUTDIR}/"logfile.err.log"

function add_chr_correctly () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME on $BN at "`date`
    echo "-----------------------------------------------------------------"
    sed -E '/^[0-9XYM]/ s/^/chr/g; s/chrMT/chrM/g; s/(^[A-Z0-9]*)\.([1234])/\1-\2/g ' $GTF_DIR/Sus_scrofa.Sscrofa10.2.84.gtf > $GTF_DIR/Sus_scrofa.Sscrofa10.2.84.chr_added_correctly.gtf
    awk '$3 != "gene" {print}' $GTF_DIR/Sus_scrofa.Sscrofa10.2.84.chr_added_correctly.gtf > $GTF_DIR/Sus_scrofa.Sscrofa10.2.84.with.chr.no_gene.gtf
}


function get_fastq_files_alt () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME on $BN at "`date`
    echo "-----------------------------------------------------------------"
    echo samtools view -b $DATA/${BN}$SUFFIX $TCRD \> ${OUTDIR}/${BN}_delta.bam >> $LOGFILE
    samtools view -b $DATA/${BN}$SUFFIX $TCRD > ${OUTDIR}/${BN}_delta.bam
    echo samtools view -b $DATA/${BN}$SUFFIX $TCRG \> ${OUTDIR}/${BN}_gamma.bam >> $LOGFILE
    samtools view -b $DATA/${BN}$SUFFIX $TCRG > ${OUTDIR}/${BN}_gamma.bam
    for TCR in delta gamma
    do
	
	echo  python ${SCRIPTS}/fix_bam_set_readfile.py --input $DATA/${BN}$SUFFIX --filter ${OUTDIR}/${BN}_${TCR}.bam --out ${OUTDIR}/${BN}_${TCR} --unmapped >> $LOGFILE
	python ${SCRIPTS}/fix_bam_set_readfile.py --input $DATA/${BN}$SUFFIX --filter ${OUTDIR}/${BN}_${TCR}.bam --out ${OUTDIR}/${BN}_${TCR} --unmapped >> $ERRORFILE	
	echo samtools view -bf 4 -f 8 -F 256 ${OUTDIR}/${BN}_${TCR}_clean.bam  \| samtools fastq -1 ${OUTDIR}/${BN}_${TCR}_1.fq -2 ${OUTDIR}/${BN}_${TCR}_2.fq -0  ${OUTDIR}/singles.fq  - 2\>\> $ERRORFILE >> $LOGFILE
	samtools view -bf 4 -f 8 -F 256 ${OUTDIR}/${BN}_${TCR}_clean.bam  | samtools fastq -1 ${OUTDIR}/${BN}_${TCR}_1.fq -2 ${OUTDIR}/${BN}_${TCR}_2.fq -0  ${OUTDIR}/singles.fq - 2>> $ERRORFILE
	rm -f ${OUTDIR}/${BN}_${TCR}_clean.bam
	echo bwa mem -M -t 10 /Volumes/Temp/Lukas/Hammer/annotations/susScr3_${TCR}.fa ${OUTDIR}/${BN}_${TCR}_1.fq ${OUTDIR}/${BN}_${TCR}_2.fq 2\>\> $ERRORFILE \| samtools view -Sb -F 12  -F 256 - \> ${OUTDIR}/${BN}_${TCR}_clean.bam  >> $LOGFILE
	bwa mem -M -t 10 /Volumes/Temp/Lukas/Hammer/annotations/susScr3_${TCR}.fa ${OUTDIR}/${BN}_${TCR}_1.fq ${OUTDIR}/${BN}_${TCR}_2.fq  2>> $ERRORFILE | samtools view -Sb -F 12 -F 256 - > ${OUTDIR}/${BN}_${TCR}_clean.bam
	echo samtools merge ${OUTDIR}/${BN}_${TCR}_all.bam ${OUTDIR}/${BN}_${TCR}_filt.bam ${OUTDIR}/${BN}_${TCR}_clean.bam >> $LOGFILE
	samtools merge ${OUTDIR}/${BN}_${TCR}_all.bam ${OUTDIR}/${BN}_${TCR}_filt.bam ${OUTDIR}/${BN}_${TCR}_clean.bam
	echo samtools fastq -0 /dev/null -1 ${OUTDIR}/${BN}_${TCR}_1.fq -2 ${OUTDIR}/${BN}_${TCR}_2.fq ${OUTDIR}/${BN}_${TCR}_all.bam >> $LOGFILE
	samtools fastq -0 /dev/null -1 ${OUTDIR}/${BN}_${TCR}_1.fq -2 ${OUTDIR}/${BN}_${TCR}_2.fq ${OUTDIR}/${BN}_${TCR}_all.bam
	$BBFIX in=${OUTDIR}/${BN}_${TCR}_1.fq in2=${OUTDIR}/${BN}_${TCR}_2.fq out=${OUTDIR}/${BN}_${TCR}_1.fq.gz out2=${OUTDIR}/${BN}_${TCR}_2.fq.gz overwrite=t 2>> $ERRORFILE
	rm -f ${OUTDIR}/${BN}_${TCR}_1.fq ${OUTDIR}/${BN}_${TCR}_2.fq
	rm -f ${OUTDIR}/${BN}_${TCR}_clean.bam ${OUTDIR}/${BN}_${TCR}_filt.bam ${OUTDIR}/${BN}_${TCR}_all.bam
	exit
    done
}



function get_fastq_files () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME on $BN at "`date`
    echo "-----------------------------------------------------------------"
    echo samtools view -bf 4 $DATA/${BN}$SUFFIX \> ${OUTDIR}/${BN}_unmapped.bam >> $LOGFILE
    samtools view -bf 4 $DATA/${BN}$SUFFIX > ${OUTDIR}/${BN}_unmapped.bam
    echo samtools view -b $DATA/${BN}$SUFFIX $TCRD \> ${OUTDIR}/${BN}_delta.bam >> $LOGFILE
    samtools view -b $DATA/${BN}$SUFFIX $TCRD > ${OUTDIR}/${BN}_delta.bam
    echo samtools view -b $DATA/${BN}$SUFFIX $TCRG \> ${OUTDIR}/${BN}_gamma.bam >> $LOGFILE
    samtools view -b $DATA/${BN}$SUFFIX $TCRG > ${OUTDIR}/${BN}_gamma.bam
    for TCR in delta gamma
    do
	echo samtools merge ${OUTDIR}/${BN}_${TCR}_all.bam ${OUTDIR}/${BN}_${TCR}.bam ${OUTDIR}/${BN}_unmapped.bam >> $LOGFILE
	samtools merge ${OUTDIR}/${BN}_${TCR}_all.bam ${OUTDIR}/${BN}_${TCR}.bam ${OUTDIR}/${BN}_unmapped.bam
	echo samtools fastq -0 /dev/null -1 ${OUTDIR}/${BN}_${TCR}_1.fq -2 ${OUTDIR}/${BN}_${TCR}_2.fq ${OUTDIR}/${BN}_${TCR}_all.bam >> $LOGFILE
	samtools fastq -0 /dev/null -1 ${OUTDIR}/${BN}_${TCR}_1.fq -2 ${OUTDIR}/${BN}_${TCR}_2.fq ${OUTDIR}/${BN}_${TCR}_all.bam
	$BBFIX in=${OUTDIR}/${BN}_${TCR}_1.fq in2=${OUTDIR}/${BN}_${TCR}_2.fq out=${OUTDIR}/${BN}_${TCR}_1.fq.gz out2=${OUTDIR}/${BN}_${TCR}_2.fq.gz 2>> $ERRORFILE
	rm -f ${OUTDIR}/${BN}_${TCR}_1.fq ${OUTDIR}/${BN}_${TCR}_2.fq ${OUTDIR}/${BN}_${TCR}.bam
    done
    rm -f ${OUTDIR}/${BN}_unmapped.bam
}

function bbduk_trimming () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME on $BN at "`date`
    echo "-----------------------------------------------------------------"
    echo "starting bbduk at:"`date`  >> $LOGFILE
    # run bbduk and trim overlapping parts, remove adaptors
    echo trimming and removing adaptors from $BN1 at `date` >> ${LOGFILE}
    echo $BBDUK -Xmx4g in1=${OUTDIR}/${BN1}_1.fq.gz in2=${OUTDIR}/${BN1}_2.fq.gz out1=${OUTDIR}/${BN1}_na_1.fq.gz out2=${OUTDIR}/${BN1}_na_2.fq.gz minlen=50 ktrim=r ref=\"$ADAPTS\" k=23 mink=11 hdist=1 overwrite=t tbo=t tpe=t \>\> ${LOGFILE} 2\> ${ERRORFILE} >> ${LOGFILE}
     $BBDUK -Xmx4g in1=${OUTDIR}/${BN1}_1.fq.gz in2=${OUTDIR}/${BN1}_2.fq.gz out1=${OUTDIR}/${BN1}_na_1.fq.gz out2=${OUTDIR}/${BN1}_na_2.fq.gz minlen=50 ktrim=r ref=\"$ADAPTS\" k=23 mink=11 hdist=1 overwrite=t tbo=t tpe=t >> ${LOGFILE} 2>> ${ERRORFILE}
    ES=$?
    echo finished at `date` with exit state $ES >> $LOGFILE
    [ $ES -eq 0 ] || exit $ES
    echo finished at `date` with exit state $ES
    rm -f ${OUTDIR}/${BN1}_1.fq.gz ${OUTDIR}/${BN1}_2.fq.gz

}

function check_bam_index () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME on $BAM at "`date`
    echo "-----------------------------------------------------------------"
    if [ ! -e $BAM".bai" ]
    then
	echo "indexing $BAM at `date`" >> $LOGFILE
	echo samtools index $BAM  >> $LOGFILE
	samtools index $BAM
	ES=$?
	echo finished at `date` with exit state $ES >> $LOGFILE
	[ $ES -eq 0 ] || exit $ES
    fi    
}

function trinity () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME on $BAM at "`date`
    echo "-----------------------------------------------------------------"
    echo "starting trinity for $BN1 at:"`date`  >> $LOGFILE
    # not sure but the strandedness should be RF for firststrand (which it seems to be :( )
    echo $TRINITY --bflyHeapSpaceMax 16G --bflyHeapSpaceInit 1G --full_cleanup --seqType fq --max_memory 160G --CPU 40 --SS_lib_type RF --left ${OUTDIR}/${BN1}_na_1.fq.gz --right ${OUTDIR}/${BN1}_na_2.fq.gz --output ${OUTDIR}/${BN1}_trinity 2\>\> ${ERRORFILE}_trinity \>\> ${LOGFILE}_trinity >> $LOGFILE
    $TRINITY --bflyHeapSpaceMax 16G --bflyHeapSpaceInit 1G --full_cleanup --seqType fq --max_memory 160G --CPU 40 --SS_lib_type RF --left ${OUTDIR}/${BN1}_na_1.fq.gz --right ${OUTDIR}/${BN1}_na_2.fq.gz --output ${OUTDIR}/${BN1}_trinity 2>> ${ERRORFILE}_trinity >> ${LOGFILE}_trinity
    ES=$?
    echo finished at `date` with exit state $ES >> $LOGFILE
#	[ $ES -eq 0 ] || exit $ES
    
}

function bridger () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME on $BAM at "`date`
    echo "-----------------------------------------------------------------"
    echo "starting bridger for $BN1 at:"`date`  >> $LOGFILE
    echo $BRIDGER --kmer_length 30  --clean --seqType fq --CPU 40 --SS_lib_type RF --left ${OUTDIR}/${BN1}_na_1.fq.gz --right ${OUTDIR}/${BN1}_na_2.fq.gz --output ${OUTDIR}/${BN1}_bridger 2\>\> ${ERRORFILE}_bridger \>\> ${LOGFILE}_bridger >> $LOGFILE
    $BRIDGER --kmer_length 30  --clean --seqType fq --CPU 40 --SS_lib_type RF --left ${OUTDIR}/${BN1}_na_1.fq.gz --right ${OUTDIR}/${BN1}_na_2.fq.gz --output ${OUTDIR}/${BN1}_bridger 2>> ${ERRORFILE}_bridger >> ${LOGFILE}_bridger
    ES=$?
    echo finished at `date` with exit state $ES >> $LOGFILE
#	[ $ES -eq 0 ] || exit $ES
    
}


function binpacker () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME on $BAM at "`date`
    echo "-----------------------------------------------------------------"
    echo "starting binpacker for $BN1 at:"`date`  >> $LOGFILE
    echo $BINPACKER -p pair -s fq -m RF -l ${OUTDIR}/${BN1}_na_1.fq.gz -r ${OUTDIR}/${BN1}_na_2.fq.gz -o ${OUTDIR}/${BN1}_binpacker \>\> ${LOGFILE}_binpacker 2\>\> ${ERRORFILE}_binpacker >> $LOGFILE
     $BINPACKER -p pair -s fq -m RF -l ${OUTDIR}/${BN1}_na_1.fq.gz -r ${OUTDIR}/${BN1}_na_2.fq.gz -o ${OUTDIR}/${BN1}_binpacker >> ${LOGFILE}_binpacker 2>> ${ERRORFILE}_binpacker
    ES=$?
    echo finished at `date` with exit state $ES >> $LOGFILE
#	[ $ES -eq 0 ] || exit $ES
    
}

function gmap_pig () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME on $BAM at "`date`
    echo "-----------------------------------------------------------------"
    echo "starting gmap for $ASSEMBLY at:"`date`  >> $LOGFILE
    echo gmap -D /Volumes/Temp/Lukas/Hammer/annotations -d susScr3 -n 0 -t 10 -f samse $ASSEMBLY 2\> ${ASSEMBLY}.gmap.error \| samtools view -S - \| samtools sort -o ${ASSEMBLY}.gmap.bam -  >> $LOGFILE
     gmap -D /Volumes/Temp/Lukas/Hammer/annotations -d susScr3 -n 0 -t 10 -f samse $ASSEMBLY 2> ${ASSEMBLY}.gmap.error | samtools view -Sb - | samtools sort -o ${ASSEMBLY}.gmap.bam -
    ES=$?
    echo finished at `date` with exit state $ES >> $LOGFILE
    samtools index ${ASSEMBLY}.gmap.bam
#	[ $ES -eq 0 ] || exit $ES
    
}


function pblat_pig () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME on $BAM at "`date`
    echo "-----------------------------------------------------------------"
    echo "starting pblat for $ASSEMBLY at:"`date`  >> $LOGFILE
    if [[ $TCR == "delta" ]]
    then
	CHR_NAME=susScr3_chr7
    else
	CHR_NAME=susScr3_chr9
    fi
    
    if [[ ! -e $REF_DIR/${CHR_NAME}_11.ooc ]]
    then
	echo pblat -t=dna -q=rna  -threads=25 -out=pslx -makeOoc=$REF_DIR/${CHR_NAME}_11.ooc $REF_DIR/${CHR_NAME}.fa $ASSEMBLY ${ASSEMBLY}.pslx 2\> ${ASSEMBLY}.pblat.error.log >> $LOGFILE
	pblat -t=dna -q=rna  -threads=25 -out=pslx -makeOoc=$REF_DIR/${CHR_NAME}_11.ooc $REF_DIR/${CHR_NAME}.fa $ASSEMBLY ${ASSEMBLY}.pslx 2> ${ASSEMBLY}.pblat.error.log
    fi
    echo pblat  -t=dna -q=rna  -threads=25 -out=pslx -ooc=$REF_DIR/${CHR_NAME}_11.ooc $REF_DIR/${CHR_NAME}.fa $ASSEMBLY ${ASSEMBLY}.pslx  2\>\> ${ASSEMBLY}.pblat.error.log >> $LOGFILE
    pblat  -t=dna -q=rna  -threads=25 -out=pslx -ooc=$REF_DIR/${CHR_NAME}_11.ooc $REF_DIR/${CHR_NAME}.fa $ASSEMBLY ${ASSEMBLY}.pslx  2>> ${ASSEMBLY}.pblat.error.log
    echo pslSort g2g ${ASSEMBLY}.sorted.fa.pslx /tmp/tempDir ${ASSEMBLY}.pslx  2\>\> ${ASSEMBLY}.pblat.error.log >> $LOGFILE
    pslSort g2g ${ASSEMBLY}.sorted.fa.pslx /tmp/tempDir ${ASSEMBLY}.pslx  2>> ${ASSEMBLY}.pblat.error.log
    echo pslReps ${ASSEMBLY}.sorted.fa.pslx ${ASSEMBLY}.best.fa.pslx ${ASSEMBLY}.sorted.fa.psr  2\>\> ${ASSEMBLY}.pblat.error.log >> $LOGFILE
    pslReps ${ASSEMBLY}.sorted.fa.pslx ${ASSEMBLY}.best.fa.pslx ${ASSEMBLY}.sorted.fa.psr  2>> ${ASSEMBLY}.pblat.error.log
    ES=$?
    echo finished at `date` with exit state $ES >> $LOGFILE
#	[ $ES -eq 0 ] || exit $ES
    
}

function cleanseq () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME on $BAM at "`date`
    echo "-----------------------------------------------------------------"
    echo "starting cleanseq for $ASSEMBLY at:"`date`  >> $LOGFILE
    AS_BN=${ASSEMBLY//\.fa*}
    echo ~/Tools/seqclean-x86_64/seqclean $ASSEMBLY -c 15 -r ${AS_BN}.clean.report -o ${AS_BN}.clean.fa  2\> ${ASSEMBLY}.clean.error.log >> $LOGFILE
    ~/Tools/seqclean-x86_64/seqclean $ASSEMBLY -c 15 -r ${AS_BN}.clean.report -o ${AS_BN}.clean.fa 2> ${AS_BN}.clean.error.log
    ES=$?
    echo finished at `date` with exit state $ES >> $LOGFILE
    rm -rf cleaning_*
#	[ $ES -eq 0 ] || exit $ES
    
}


#for BAM in ${DATA}/*$SUFFIX
for BAM in  ${DATA}/[ST]*$SUFFIX
do
    BN=`basename $BAM $SUFFIX`
    LOGFILE=${OUTDIR}/${BN}"_logfile.log"
    ERRORFILE=${OUTDIR}/${BN}"_error.log"
    #check_bam_index
    #get_fastq_files
    #get_fastq_files_alt
    for TCR in delta gamma;
    do
	BN1=${BN}_$TCR
	#bbduk_trimming
	trinity	
	ASSEMBLY=${OUTDIR}/${BN1}_trinity.Trinity.fasta
	cleanseq
	ASSEMBLY=${AS_BN}.clean.fa	
	gmap_pig 
	pblat_pig 
	binpacker
	ASSEMBLY=${OUTDIR}/${BN1}_binpacker/BinPacker.fa
	cleanseq
	ASSEMBLY=${AS_BN}.clean.fa	
	gmap_pig
	pblat_pig 
    done
done


exit

# Assemble all:
# create joined fq files:
find -name "*gamma_na_1.fq.gz" -exec zcat {} + | gzip -c > all_gamma_na_1.fq.gz &
find -name "*delta_na_1.fq.gz" -exec zcat {} + | gzip -c > all_delta_na_1.fq.gz &
find -name "*delta_na_2.fq.gz" -exec zcat {} + | gzip -c > all_delta_na_2.fq.gz &
find -name "*gamma_na_2.fq.gz" -exec zcat {} + | gzip -c > all_gamma_na_2.fq.gz &

#for BAM in ${DATA}/*$SUFFIX
for BAM in "all"
do
    BN=`basename $BAM $SUFFIX`
    LOGFILE=${OUTDIR}/${BN}"_logfile.log"
    ERRORFILE=${OUTDIR}/${BN}"_error.log"
    #check_bam_index
    #get_fastq_files
    #get_fastq_files_alt
    for TCR in delta gamma;
    do
	BN1=${BN}_$TCR
	#bbduk_trimming
	trinity	
	ASSEMBLY=${OUTDIR}/${BN1}_trinity.Trinity.fasta
	cleanseq
	ASSEMBLY=${AS_BN}.clean.fa	
	gmap_pig 
	pblat_pig 
	binpacker
	ASSEMBLY=${OUTDIR}/${BN1}_binpacker/BinPacker.fa
	cleanseq
	ASSEMBLY=${AS_BN}.clean.fa	
	gmap_pig
	pblat_pig 
    done
done

exit


~/Tools/TransDecoder-3.0.0/TransDecoder.LongOrfs -S -t ../CD2min_delta_trinity.Trinity.clean.fa -m 60
# only keep
hmmscan --cpu 16 --domtblout pfam.domtblout /Volumes/Temp/Trinotate/Pfam-A.hmm CD2min_delta_trinity.Trinity.clean.fa.transdecoder_dir/longest_orfs.pep 2> hmmscan.err.out &
~/Tools/TransDecoder-3.0.0/TransDecoder.Predict -t ../CD2min_delta_trinity.Trinity.clean.fa --retain_pfam_hits pfam.domtblout
grep "^>" BinPacker.fa | sed 's/>//g' | awk -F "." '{print $0.$2"\t"$0}' > gene_transcript_map.txt
hmmscan --cpu 8 --domtblout pfam.domtblout /Volumes/Temp/Trinotate/Pfam-A.hmm BinPacker.clean.fa.transdecoder_dir/longest_orfs.pep > hmmscan.err.out &

~/Tools/TransDecoder-3.0.0/TransDecoder.LongOrfs --gene_trans_map gene_transcript_map.txt  -t BinPacker.clean.fa -m 60
cat BinPacker.clean.fa.transdecoder_dir/longest_orfs.pep | parallel -k -j 6 -N 200 --recstart '>' --pipe hmmscan --cpu 8 --domtblout pfam{#}.domtblout -o hmmscan{#}.err.out /Volumes/Temp/Trinotate/Pfam-A.hmm -
cat hmmscan[0-9]*.err.out > hmmscan.err.out
rm -f hmmscan[0-9]*.err.out
cat pfam[0-9]*.domtblout > pfam.domtblout
rm -f pfam[0-9]*.domtblout
~/Tools/TransDecoder-3.0.0/TransDecoder.Predict  --gene_trans_map gene_transcript_map.txt -t BinPacker.clean.fa --retain_pfam_hits pfam.domtblout


# Assemble all:
# create joined fq files:
find -name "*gamma_na_1.fq.gz" -exec zcat {} + | gzip -c > all_gamma_na_1.fq.gz &
find -name "*delta_na_1.fq.gz" -exec zcat {} + | gzip -c > all_delta_na_1.fq.gz &
find -name "*delta_na_2.fq.gz" -exec zcat {} + | gzip -c > all_delta_na_2.fq.gz &
find -name "*gamma_na_2.fq.gz" -exec zcat {} + | gzip -c > all_gamma_na_2.fq.gz &


## create gtf for gamma and delta exons:
grep "ENSSSCG00000023584\|ENSSSCG00000022512" Sus_scrofa.Sscrofa10.2.84.with.chr.no_gene.gtf | awk '$3=="exon"' > ENSSSCG00000023584_ENSSSCG00000022512_exons.gtf

wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=pig[orgn]+NOT+srcdb_refseq_predicted[PROP]&usehistory=y" -O pig_query.xml
less pig_query.xml
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&WebEnv=NCID_1_111950462_130.14.22.215_9001_1471427377_750776235_0MetA0_S_MegaStore_F_1&query_key=1&rettype=fasta&retmode=text" -O pig_proteins.fasta
makeblastdb -help
makeblastdb -in pig_proteins.fasta -dbtype prot -title pig_proteins -hash_index -parse_seqids -out pig_proteins


wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=pig[orgn]+NOT+srcdb_refseq_predicted[PROP]+AND+%22t cell receptor%22[title] &usehistory=y"  -O pig_tcr_query.xml
 wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&WebEnv=NCID_1_256692318_130.14.18.34_9001_1471438861_1446048850_0MetA0_S_MegaStore_F_1&query_key=1&rettype=acc&retmode=text" -O  pig_tcr.acc

blastx -db /Volumes/Temp/BlastDB/pig_proteins -evalue 1e-3 -strand plus -query all_gamma_trinity.Trinity.clean.fa -max_target_seqs 100 -outfmt 7 -num_threads 10 > test.out

grep -f "/Volumes/Temp/BlastDB/pig_tcr.acc" test.out  | cut -f 1 | sort -u > test.readids

samtools view all_gamma_trinity.Trinity.clean.fa.gmap.bam | grep -f test.readids - | awk '{print ">"$1"\t"$3":"$4"\n"$10}' > test.reads.fna
