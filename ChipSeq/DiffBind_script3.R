################ ChiqSeq with Peakfiles #######################
#setwd("C://Users/Chris/Documents/VetMed_Master/Bioinformatics/Project/")
library("ChIPQC")
library("DiffBind")
library(rtracklayer)
wd<-getwd()
wd<-paste(wd,"data/05_bowtieUcscHg19/prep/bam_peak_calling", sep="/")
wd

# A) generate .txt file > SampleID.txt | Unix
# B) generate .txt file > bam.files.txt | Unix
# C) generate .txt file > peaks.txt | Unix
wd <- "/Volumes/vetlinux01/ChipSeq/PBI_ChipSeq/"
Path.sample<-paste(wd,"data",sep="/") 
Path.bam<-paste(wd,"data",sep="/") 
Path.control<-paste(wd,"data/SRR2000691__SRR2000692__SRR2000693_trimmo.bam",sep="/")
Path.peaks<-paste(wd,"Peaks",sep="/")
Path.general<-wd
#-------------------

bam_files=grep("692",list.files(Path.sample,pattern="*_trimmo.bam$"),invert=TRUE,value=TRUE)
peak_files=grep("692",list.files(Path.peaks,pattern="*.xls$"),invert=TRUE,value=TRUE)


SampleID<-c("ER_e_1","ER_e_2","ER_e_3","ER_ep_1","ER_ep_2","ER_ep_3","PR_e_1","PR_e_2","PR_e_3","PR_ep_1","PR_ep_2","PR_ep_3")

Tissue<-rep("T47D",12)

Factor<-c("ER","ER","ER","ER","ER","ER",
          "PR","PR","PR","PR","PR","PR")

Condition<-rep("Resistant",12)

Treatment<-c("e","e","e","ep","ep","ep","e","e","e","ep","ep","ep")

Replicate<-c("1","2","3","1","2","3","1","2","3","1","2","3")

bamRead<-paste(rep(Path.sample,12),bam_files,sep = "/")

ControlID<-rep("Ctrl.ID",12)

Peaks<-paste(Path.peaks,peak_files,sep = "/")

Peakcaller<-rep("macs",12)

bamControl<-rep(Path.control,12)

datasheet<-data.frame(SampleID=SampleID,Tissue=Tissue,Factor=Factor,Condition=Condition,Treatment=Treatment,
                      Replicate=Replicate,bamReads=bamRead,ControlID=ControlID,bamControl=bamControl,
                      Peaks=Peaks,PeakCaller=Peakcaller)

datasheet$Replicate<-as.integer(datasheet$Replicate)

chip_experiment<-dba(sampleSheet = datasheet) 

# peak overlaps
overlap_rate=dba.overlap(chip_experiment,mode = DBA_OLAP_RATE)
plot(1:length(overlap_rate),overlap_rate,log="y",ty="b")

#dba_object created

#heatmap

pdf("clustering_heatmap.pdf")
plot(chip_experiment)
dev.off()

#pca plot

pdf("pca_plot.pdf")
dba.plotPCA(chip_experiment,attributes=c(DBA_TREATMENT,DBA_FACTOR))
dev.off()

#venn plot

pdf("venn_plot1.pdf")
par(mfrow=c(2,2))
dba.plotVenn(chip_experiment,mask=c(1:3),main="")
dba.plotVenn(chip_experiment,mask=c(4:6),main="")
dba.plotVenn(chip_experiment,mask=c(7:9),main="")
dba.plotVenn(chip_experiment,mask=c(10:12),main="") #attributes=c(DBA_TREATMENT,DBA_FACTOR))
title(main="Binding Sites Overlap among replicates",outer=T)
dev.off()


#diff bind analysis
#occupancy analysis
consensus_peaks<-dba.peakset(chip_experiment,consensus=c(DBA_FACTOR,DBA_TREATMENT),minOverlap = 2)
dba.plotVenn(consensus_peaks,consensus_peaks$masks$Consensus & consensus_peaks$masks$ER)
consensus_dba<-dba(consensus_peaks,mask=consensus_peaks$masks$Consensus)

pdf("venn_qualitative.pdf")
par(mfrow=c(1,2))
dba.plotVenn(consensus_dba,mask=c(1:2),main="")
dba.plotVenn(consensus_dba,mask=c(3:4),main="")
title(main="New binding sites",outer=T)
dev.off()

overlap_ER<-dba.overlap(consensus_dba,mask = c(1:2))
overlap_ER$onlyB
overlap_PR<-dba.overlap(consensus_dba,mask = c(3:4))
overlap_PR$onlyB

consensus_dba<-dba.peakset(DBA=consensus_dba,peaks=overlap_ER$onlyB,sampID="overlap_ER_only_ep",tissue="T47D",factor="ER",condition="Resistant",treatment="ep",replicate=1,bamReads=NA,control=NA,bamControl=NA,peak.caller="macs")

consensus_dba<-dba.peakset(DBA=consensus_dba,peaks=overlap_PR$onlyB,sampID="overlap_PR_only_ep",tissue="T47D",factor="PR",condition="Resistant",treatment="ep",replicate=1,bamReads=NA,control=NA,bamControl=NA,peak.caller="macs")

dba.plotVenn(consensus_dba,mask=c(5:6))

#overlap of native progesterone and EPM medium

pdf("overlap_venn.pdf")
par(mfrow=c(1,1))
dba.plotVenn(new_binding2,mask=c(1:2),main="Most binding sites that overlap are from native progesterone targets")
dev.off()

overlap_prog<-dba.overlap(consensus_dba,mask=c(5:6))
#output the bed file for motif calling
export(overlap_prog$inAll,"progesterone_overlapping_sites.bed",format="BED")


#quantitative analysis
ER_peakset<-dba(chip_experiment,mask=c(1:6))
PR_peakset<-dba(chip_experiment,mask=c(7:12))

ER_peakset1<-dba(chip_experiment,mask=c(1:6))
PR_peakset1<-dba(chip_experiment,mask=c(7:12))

ER_peakset<-dba.peakset(ER_peakset,consensus = DBA_TISSUE,minOverlap = 2)
PR_peakset<-dba.peakset(PR_peakset,consensus = DBA_TISSUE,minOverlap = 2)

ER_peakset_dba<-dba.peakset(ER_peakset,c(7),bRetrieve = T)
PR_peakset_dba<-dba.peakset(PR_peakset,c(7),bRetrieve = T)

#dba.count

ER_peakset_dba<-dba.count(ER_peakset1,peaks=ER_peakset_dba,score = DBA_SCORE_READS)
PR_peakset_dba<-dba.count(PR_peakset1,peaks=PR_peakset_dba,score = DBA_SCORE_READS)

#plot correlation of counted reads in peaks -- before: pvalues

pdf("clustering_after_count.pdf")
plot(ER_peakset_dba)
plot(PR_peakset_dba)
dev.off()

#define contrasts---
ER_peakset_dba<-dba.contrast(ER_peakset_dba,categories = DBA_TREATMENT)
PR_peakset_dba<-dba.contrast(PR_peakset_dba,categories = DBA_TREATMENT)

#analysis of differential binding
ER_peakset_dba<-dba.analyze(ER_peakset_dba,method=DBA_DESEQ2, bSubControl = TRUE, bFullLibrarySize = TRUE, bTagwise = TRUE)
PR_peakset_dba<-dba.analyze(PR_peakset_dba,method=DBA_DESEQ2, bSubControl = TRUE, bFullLibrarySize = TRUE, bTagwise = TRUE)

#get diffbind in report

ER_peakset_report <- dba.report(ER_peakset_dba,contrast=1,method=DBA_DESEQ2,th=0.05,bUsePval=FALSE,
                                fold=log2(2),DataType=DBA_DATA_GRANGES,bNormalized = T)

PR_peakset_report <- dba.report(PR_peakset_dba,contrast=1,method=DBA_DESEQ2,th=0.05,bUsePval=FALSE,
                                fold=log2(2),DataType=DBA_DATA_GRANGES,bNormalized = T)


#extract only upregulated

ER_peakset_report <- ER_peakset_report[mcols(ER_peakset_report)[,4] < -1] 
PR_peakset_report <- PR_peakset_report[mcols(PR_peakset_report)[,4] < -1] 


#overlap cotreated regions
overlap_upreg <- dba.peakset(NULL, ER_peakset_report, sampID = "ER_UP",tissue="T47D",factor="ER",condition="Resistant",
                             treatment="noTreat",replicate = "1-2-3", control="noInp",peak.caller="macs")

overlap_upreg2<-dba.peakset(overlap_upreg, PR_peakset_report, sampID="PR_UP",tissue="T47D",factor="PR",condition="Resistant",
                            treatment="noTreat",replicate="1-2-3",control="noInp",peak.caller = "macs")


pdf("overlap_venn_quantitative.pdf")
par(mfrow=c(1,1))
dba.plotVenn(overlap_upreg2,mask=c(1:2),main="Most of the ER acquired binding sites overlap with the progesterone regulated genes")
dev.off()

# extract overlapping peakset and export as bed.file 
overlap_upreg2_extract <- dba.overlap(overlap_upreg2, c(1:2))

write.table(overlap_upreg2_extract$inAll,"quantitative_overlap",sep="\t",row.names = F,quote = F,col.names = F)
extract<-read.delim("quantitative_overlap",header = F)
extract<-extract[,1:3]
vec2<-seq(1,length(extract$V1),by=1)
vec2<-paste("A",vec2,sep="")
extract<-cbind(extract,vec2)
write.table(extract,"quantitative_overlapping_sites.bed",sep="\t",row.names = F,quote = F,col.names = F)


