#!/bin/bash
#----------
# author: Lukas Endler
# authored: 2016-05-02 12:01:03 
# Time-stamp: <2016-07-05 11:55:08 lukasendler>
# command lines for cufflinks
#--------------

DATA=/Volumes/Temp/Hammer/data/
SCRIPTS=/Volumes/Temp/Hammer/scripts/
OUTDIR=/Volumes/Temp/Hammer/data/STAR_2nd_pass
CUFFDIR=/Volumes/Temp/Lukas/Tools/cufflinks-2.2.1/src
REF_GTF=/Volumes/Temp/Hammer/annotations/Sus_scrofa.Sscrofa10.2.84.with.chr.no_gene.gtf
MERGE_GTF=/Volumes/Temp/Hammer/Cufflinks/cuffmerge/merged_asm/merged.gtf
REFGENOME=/Volumes/Temp/Hammer/annotations/susScr3.fa
GEN_DIR=/Volumes/Temp/Hammer/annotations/genomeDir_2pass
LOGFILE=${OUTDIR}/"logfile.log"
ERRLOG=${OUTDIR}/"logfile.err.log"

function PrepareReference_2nd_pass () {
    echo "$FUNCNAME on $BN at "`date` >> $LOGFILE
    CURRDIR=`pwd`
    cd $annotationDir	
    #ftp://ftp.ensembl.org/pub/release-84/gtf/sus_scrofa/Sus_scrofa.Sscrofa10.2.84.gtf.gz 8mar2016
    #gtfFile=Sus_scrofa.Sscrofa10.2.84
    fastaFile=susScr3.fa
    #awk 'BEGIN{FS=OFS="\t"}; NR >6 {print "chr"$1, $2, $3,$4, $5, $6, $7, $8,$9}' $gtfFile.gtf > $gtfFile.with.chr.gtf
    #head $gtfFile.gtf
    #head $gtfFile.with.chr.gtf
    
    if [ -d $GEN_DIR ]; then
	echo rm $GEN_DIR >> $LOGFILE
	rm $GEN_DIR
    fi 
    echo mkdir $GEN_DIR >> $LOGFILE
    mkdir $GEN_DIR
    cat $REF_GTF $MERGE_GTF > $OUTDIR/merged.gtf
    echo 	/Volumes/Temp/Hammer/STAR-2.5.2a/bin/MacOSX_x86_64/STAR --runMode genomeGenerate --genomeDir $GEN_DIR \
		--genomeFastaFiles $REFGENOME \
		--sjdbGTFfile  $OUTDIR/merged.gtf \
		--sjdbOverhang 124 \
		--runThreadN 20 >> $LOGFILE  
    /Volumes/Temp/Hammer/STAR-2.5.2a/bin/MacOSX_x86_64/STAR --runMode genomeGenerate --genomeDir $GEN_DIR \
							    --genomeFastaFiles $REFGENOME \
							    --sjdbGTFfile  $OUTDIR/merged.gtf \
							    --sjdbOverhang 124 \
							    --runThreadN 20 2>> $ERRLOG >> $LOGFILE
    ES=$?
    echo finished at `date` with exit state $ES >> $LOGFILE
    cd $CURRDIR
    [ $ES -eq 0 ] || exit $ES
    
}

function Mapping () {
    echo "FUNCTION: $FUNCNAME at "`date` >> $LOGFILE 
    cd $OUTDIR
    tissues=( Liver_ CD2min_ CD2plus_ IEL_ Lung_ PBMC_ Spleen_ Thymus_ )
    for ID in "${tissues[@]}"
    do
	echo "mapping: $ID at "`date` >> $LOGFILE 	
	echo /Volumes/Temp/Hammer/STAR-2.5.2a/bin/MacOSX_x86_64/./STAR --runMode alignReads \
								  --genomeDir $GEN_DIR \
								  --runThreadN 20 \
								  --readFilesIn $DATA/${ID}R1.fastq $DATA/${ID}R2.fastq \
								  --outFileNamePrefix ${ID} \
								  --outSAMtype BAM Unsorted >> $LOGFILE
								  #--outSAMmapqUnique Integer0to255 \
								  #--quantMode TranscriptomeSAM >> $LOGFILE
	/Volumes/Temp/Hammer/STAR-2.5.2a/bin/MacOSX_x86_64/./STAR --runMode alignReads \
								  --genomeDir $GEN_DIR \
								  --runThreadN 20 \
								  --readFilesIn $DATA/${ID}R1.fastq $DATA/${ID}R2.fastq \
								  --outFileNamePrefix ${ID} \
								  --outSAMtype BAM Unsorted  2>> $ERRLOG >> $LOGFILE
								  #--outSAMmapqUnique Integer0to255 \
								  #--quantMode TranscriptomeSAM >> $LOGFILE
	ES=$?
	echo finished at `date` with exit state $ES >> $LOGFILE
    done
    #on the fly sorting needs far too much RAM 		
    #doesnt work with compressed fastq files also not with the corresponding option
    #do NOT USE the sorting option. needs too much memory! --outSAMtype BAM SortedByCoordinate output sorted by coordinate
}

function Mapping_Chim () {
    echo "FUNCTION: $FUNCNAME at "`date` >> $LOGFILE 
    cd $OUTDIR
    tissues=( Liver_ CD2min_ CD2plus_ IEL_ Lung_ PBMC_ Spleen_ Thymus_ )
    for ID in "${tissues[@]}"
    do
	echo "mapping: $ID at "`date` >> $LOGFILE 	
	echo /Volumes/Temp/Hammer/STAR-2.5.2a/bin/MacOSX_x86_64/./STAR --runMode alignReads \
								  --genomeDir $GEN_DIR \
								  --runThreadN 20 \
								  --readFilesIn $DATA/${ID}R1.fastq $DATA/${ID}R2.fastq \
								  --outFileNamePrefix ${ID} \
								  --chimSegmentMin 30 \
								  --outSAMtype BAM Unsorted >> $LOGFILE
								  #--outSAMmapqUnique Integer0to255 \
								  #--quantMode TranscriptomeSAM >> $LOGFILE
	/Volumes/Temp/Hammer/STAR-2.5.2a/bin/MacOSX_x86_64/./STAR --runMode alignReads \
								  --genomeDir $GEN_DIR \
								  --runThreadN 20 \
								  --readFilesIn $DATA/${ID}R1.fastq $DATA/${ID}R2.fastq \
								  --outFileNamePrefix ${ID} \
								  --chimSegmentMin 30 \
								  --outSAMtype BAM Unsorted  2>> $ERRLOG >> $LOGFILE
								  #--outSAMmapqUnique Integer0to255 \
								  #--quantMode TranscriptomeSAM >> $LOGFILE
	ES=$?
	echo finished mapping at `date` with exit state $ES >> $LOGFILE
	echo at `date`  >> $LOGFILE
	echo samtools view -hSb ${ID}_Chimeric.out.sam \| samtools sort -o ${ID}.Chimeric.sorted.bam >> $LOGFILE
	samtools view -hSb Chimeric.out.sam | samtools sort -o ${ID}.Chimeric.sorted.bam -
	samtools index ${ID}.Chimeric.sorted.bam
    done
    #on the fly sorting needs far too much RAM 		
    #doesnt work with compressed fastq files also not with the corresponding option
    #do NOT USE the sorting option. needs too much memory! --outSAMtype BAM SortedByCoordinate output sorted by coordinate
}

#Mapping
function Sorting () {
    echo "-----------------------------------------------------------------"
    echo "FUNCTION: $FUNCNAME at "`date` >> $LOGFILE 
    echo "-----------------------------------------------------------------"
    cd $OUTDIR
    tissues=( Liver_ CD2min_ CD2plus_ IEL_ Lung_ PBMC_ Spleen_ Thymus_ )
    #tissues=( CD2min_ CD2plus_ IEL_ Lung_ PBMC_ Spleen_ Thymus_ )
    for ID in "${tissues[@]}"
    do
	echo "start samtools sort $ID at"`date`	>> $LOGFILE
	echo samtools sort --threads 4 -m 3G ${ID}Aligned.out.bam -o ${ID}Aligned.sorted.bam >> $LOGFILE
	#attention no .bam for output file -m increases memory to reduce the number of intermediate files produced by samtools that throw an error message
	samtools sort --threads 4 -m 3G ${ID}Aligned.out.bam -o ${ID}Aligned.sorted.bam 2>> $ERRLOG
#	echo samtools sort -m 5000000000 ${ID}Aligned.toTranscriptome.out.bam -o ${ID}Aligned.toTranscriptome.sorted.bam >> $LOGFILE
	#attention no .bam for output file -m increases memory to reduce the number of intermediate files produced by samtools that throw an error message
#	samtools sort -m 5000000000 ${ID}Aligned.toTranscriptome.out.bam -o ${ID}Aligned.toTranscriptome.sorted.bam 2>> $ERRLOG
	ES=$?
	echo finished at `date` with exit state $ES >> $LOGFILE
	samtools index ${ID}Aligned.sorted.bam
    done
}


#PrepareReference_2nd_pass
Mapping_Chim
#Sorting
