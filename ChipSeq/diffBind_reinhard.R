# ---
# title: "diffBind"
# author: "ReinhardGrausenburger"
# date: "March 3, 2017"
# output: html_document
# ---


### load required packages

source("https://bioconductor.org/biocLite.R")
biocLite("ChIPQC")
library(ChIPQC)
browseVignettes("DiffBind")
library(DiffBind)
library(tidyverse)
library(stringr)
library(readxl)
library(tidyr)
library(GenomicRanges)
library(rtracklayer)
biocLite("rGREAT")
library(rGREAT)
### define working directory and directory for output

wdir <- "/home/grausenburgerr/Data/Vorlesungen/PraktischeBioinformatik2017/Downloads/ChIPSeqData/TestAllTillSambamba" 
setwd(wdir)

outDir <- paste(wdir, "12_diffBind", sep = "/")
dir.create(outDir)
readDir <- paste(wdir, "07B_sambOutWithDup", sep = "/")
peakDir <- paste(wdir, "09_macs2", sep = "/")
readFil <- c("SRR2000673__SRR2000674_trimmo.withdups.filtered.bam", 
  "SRR2000675__SRR2000676_trimmo.withdups.filtered.bam", 
  "SRR2000677__SRR2000678_trimmo.withdups.filtered.bam",
  "SRR2000679__SRR2000680_trimmo.withdups.filtered.bam", 
  "SRR2000681__SRR2000682_trimmo.withdups.filtered.bam", 
  "SRR2000683__SRR2000684_trimmo.withdups.filtered.bam", 
  "SRR2000712__SRR2000713_trimmo.withdups.filtered.bam", 
  "SRR2000714__SRR2000715_trimmo.withdups.filtered.bam", 
  "SRR2000716__SRR2000717_trimmo.withdups.filtered.bam", 
  "SRR2000718__SRR2000719_trimmo.withdups.filtered.bam", 
  "SRR2000720__SRR2000721_trimmo.withdups.filtered.bam", 
  "SRR2000722__SRR2000723_trimmo.withdups.filtered.bam")

peakFiles <- c("SRR2000673__SRR2000674_peaks.xls", 
  "SRR2000675__SRR2000676_peaks.xls", 
  "SRR2000677__SRR2000678_peaks.xls",
  "SRR2000679__SRR2000680_peaks.xls", 
  "SRR2000681__SRR2000682_peaks.xls", 
  "SRR2000683__SRR2000684_peaks.xls", 
  "SRR2000712__SRR2000713_peaks.xls", 
  "SRR2000714__SRR2000715_peaks.xls", 
  "SRR2000716__SRR2000717_peaks.xls", 
  "SRR2000718__SRR2000719_peaks.xls", 
  "SRR2000720__SRR2000721_peaks.xls", 
  "SRR2000722__SRR2000723_peaks.xls")
### Make sample Table:
# sampleID: unique identifiers
SampleID <- c("er_e_1", "er_e_2", "er_e_3", "er_ep_1", "er_ep_2", "er_ep_3", "pr_e_1", "pr_e_2", "pr_e_3", "pr_ep_1", "pr_ep_2", "pr_ep_3")
# column tissue, factor, condition, treatment, and replicate is used to define masks
# vector with the number of different groups top level
Tissue <- rep("T47D", 12)
# chIP-ed Factor, one chip - one level
Factor <- c(rep("er", 6), rep("pr", 6))
# factor that subgroups within top level group, e.g. genotypes in tissue could be further subdevided into male and female, 
Condition <- c(rep("e", 3), rep("ep", 3), rep("e", 3), rep("ep", 3))
# third sublevel, mice of different genotype with different sex could be treated or untreated
Treatment <- rep("noTreat", length(SampleID))
# indicates replicate
Replicate <- rep(c(1,2,3), 4)
bamReads <- str_c(readDir, readFil, sep = "/")
ControlID <- c(rep("inp_e", length(SampleID)))
# ControlID <- c(rep("inp_er_e", 3), rep("inp_er_ep", 3), rep("inp_pr_e", 3), rep("inp_pr_ep", 3)) 
bamControl <- rep(str_c(readDir, "SRR2000691__SRR2000692__SRR2000693_trimmo.withdups.filtered.bam", sep = "/"), length(SampleID))
# # bamControl <- c(rep(str_c(readDir, "SRR2000691__SRR2000692__SRR2000693_trimmo.withdups.filtered.bam", sep = "/"), 3), 
#                 rep(str_c(readDir, "SRR2000691__SRR2000692__SRR2000693_trimmo.withdups.filteredC1.bam", sep = "/"), 3),
#                 rep(str_c(readDir, "SRR2000691__SRR2000692__SRR2000693_trimmo.withdups.filteredC2.bam", sep = "/"), 3),
#                 rep(str_c(readDir, "SRR2000691__SRR2000692__SRR2000693_trimmo.withdups.filteredC3.bam", sep = "/"), 3))
                        
Peaks <- str_c(peakDir, peakFiles, sep = "/")

PeakCaller <- rep("macs", length(SampleID))
ScoreCol <- rep(7, 12) # als score nehmen wir -log10(pvalue) for the peak summit

samples <- data.frame(SampleID=SampleID, Tissue=Tissue, Factor=Factor, Condition=Condition, Treatment=Treatment, Replicate=Replicate, bamReads=bamReads, ControlID=ControlID, bamControl=bamControl, Peaks=Peaks, PeakCaller=PeakCaller, ScoreCol=ScoreCol)
View(samples)
# Overlap analysis of replicates from indivudual conditions

#### Define dba object
dbaO <- dba(sampleSheet = samples)
dbaO
#### Plot correlation of -log10 p Values of Peaks to show how samples cluster
plot(dbaO)
#*replicates cluster nicely; again, conditions are closer than the ChIPs with identical antibodies* \

#### Plot Venn diagramms of indivdual conditions to show overlap of replicates, and overlapping rates
#*only overlap of max. 4 samples is possible* 
#*er_e:*
dba.plotVenn(dbaO, c(1:3))
olap.rateWtRep <- dba.overlap(dbaO, c(1:3), mode = DBA_OLAP_RATE)
olap.rateWtRep
plot(olap.rateWtRep,type='b', ylab='# peaks', xlab='Overlap at least this many peaksets')
#*many peaks are only present in one of the three replicates. about 30000 peaks present in two of three. 12000 are in all three*

#*er_ep:*
dba.plotVenn(dbaO, c(4:6))
olap.rateWtRep <- dba.overlap(dbaO, c(4:6), mode = DBA_OLAP_RATE)
olap.rateWtRep
plot(olap.rateWtRep,type='b', ylab='# peaks', xlab='Overlap at least this many peaksets')
#there are more peaks than in er_e, about 50000 peaks present in two of three. 30000 are in all three*

#*pr_e:*
dba.plotVenn(dbaO, c(7:9))
olap.rateWtRep <- dba.overlap(dbaO, c(7:9), mode = DBA_OLAP_RATE)
olap.rateWtRep
plot(olap.rateWtRep,type='b', ylab='# peaks', xlab='Overlap at least this many peaksets')
#there are more peaks than in er_e, about 50000 peaks present in two of three. 30000 are in all three*

#*pr_e:*
dba.plotVenn(dbaO, c(7:9))
olap.rateWtRep <- dba.overlap(dbaO, c(7:9), mode = DBA_OLAP_RATE)
olap.rateWtRep
plot(olap.rateWtRep,type='b', ylab='# peaks', xlab='Overlap at least this many peaksets')
#here are far more peaks than in the samples before, about 120000 are in two, but only 3416 are overlapping in all three. pr_e_2 seems to be an outlier*

#*pr_ep:*
dba.plotVenn(dbaO, c(10:12))
olap.rateWtRep <- dba.overlap(dbaO, c(10:12), mode = DBA_OLAP_RATE)
olap.rateWtRep
plot(olap.rateWtRep,type='b', ylab='# peaks', xlab='Overlap at least this many peaksets')
#there are 200000 peaks in only 1 replicate, 80000 in at least two, and 60000 in all three*

#Overlap analysis from consensus Peak sets ER-ChIP - estrogen only, and ER-ChIP - estrogen/progesteron co-treatment#

## we look for ER binding sites that are aquired by the co-treatment ##

###consensus peakset is intersect of 3 replicates###

#### make consensus peak sets from individual replicates of er_e and er_ep
#### show with Venn Diagrams the Overlap between estrogen only and estrogen/progesterone co-treatment

#er_e und er_ep
dbaEr <- dba(dbaO, c(1:6)) # only select ERalpha ChIPs
dbaErOl3 <- dba.peakset(dbaEr, consensus = -DBA_REPLICATE, minOverlap = 3) # extract consensus peaksets
dba.show(dbaErOl3, mask = dbaErOl3$masks$Consensus) # show consensus peaksets

# this is another, more intutive way of makeing consensus peaksets
dbaErOl3 <- dba.peakset(dbaEr, consensus = DBA_CONDITION, minOverlap = 3)
dba.show(dbaErOl3, mask = dbaErOl3$masks$Consensus)

# plot consensus peaksets
dba.plotVenn(dbaErOl3, mask=dbaErOl3$masks$Consensus)
#it looks as there 3465 sites lost, and 19461 aquired upon co-treatment*

#### do analysis with a more lenient definition for consensus Peaksets
#er_e und er_ep
# this is another, more intutive way of makeing consensus peaksets
dbaErOl2 <- dba.peakset(dbaEr, consensus = DBA_CONDITION, minOverlap = 2)
dba.show(dbaErOl2, mask = dbaErOl2$masks$Consensus)

# plot consensus peaksets
dba.plotVenn(dbaErOl2, mask=dbaErOl3$masks$Consensus)
#it looks as there 11507 sites lost, and 34302 aquired upon co-treatment, so here, the difference between lost and gained binding sites is not so dramatic*

## now find PR acquired binding sites during estrogen/progesteron co-treatment
#pr_e und pr_ep
dbaPr <- dba(dbaO, c(7:12)) # only select PR ChIPs

# this is another, more intutive way of makeing consensus peaksets
dbaPrOl3 <- dba.peakset(dbaPr, consensus = DBA_CONDITION, minOverlap = 3)
dba.show(dbaPrOl3, mask = dbaPrOl3$masks$Consensus)

# plot consensus peaksets
dba.plotVenn(dbaPrOl3, mask=dbaPrOl3$masks$Consensus)
#it looks as there 1625 sites lost, and 59725 sites aquired upon co-treatment*

#### do analysis with a more lenient definition for consensus Peaksets
dbaPrOl2 <- dba.peakset(dbaPr, consensus = DBA_CONDITION, minOverlap = 2)
dba.show(dbaPrOl2, mask = dbaPrOl2$masks$Consensus)

# plot consensus peaksets
dba.plotVenn(dbaPrOl2, mask=dbaPrOl2$masks$Consensus)
#it looks as there 25685 sites lost, and 83578 sites aquired upon co-treatment, so here, the difference between lost and gained binding sites is not so dramatic as with a more strict consensus peak set*

## now check how big the overlap between ERalpha and PR binding sites upon estrogen/progesterone co-treatment is 

### first check the overlap in the more strict consensus peak set (intersection)
# I do this according to this site:
# https://support.bioconductor.org/p/67306/

# retrieve peaksets of more strict consensus peak sets (present in all three)
dbaErOl3Grl <- dba.overlap(dbaErOl3, c(7:8)) # computes overlapping and unioun sites of er_e and er_ep, and output it as a list of GRanges Objects
dbaPrOl3Grl <- dba.overlap(dbaPrOl3, c(7:8))# computes overlapping and unioun sites of pr_e and pr_ep, and output it as a list of GRanges Objects

# construct a new dba object containing consensus peaks of er_ep_unique and pr_ep_unique 
dbaErEpPrEpU <- dba.peakset(NULL, dbaErOl3Grl$onlyB, sampID = "er_ep_u", tissue = "T47D", factor = "er", condition = "ep", treatment = "noTreat", replicate = "1-2-3", control = "noInp", peak.caller = "macs") # retrieve er_ep_unique peaks, and construct with this peak set a new DBA object

dbaErEpPrEpU <- dba.peakset(dbaErEpPrEpU, dbaPrOl3Grl$onlyB, sampID = "pr_ep_u", tissue = "T47D", factor = "pr", condition = "ep", treatment = "noTreat", replicate = "1-2-3", control = "noInp", peak.caller = "macs")

dba.plotVenn(dbaErEpPrEpU, c(1,2)) # plot overlap of upon cotreatment gained pr and er binding sites

# extract overlapping peakset and export as bed.file to perform motif analysis on that peaks
dbaErEpPrEpUnOl <- dba.overlap(dbaErEpPrEpU, c(1:2))

# findMotifsGenome.pl from Homer does not work if the name of the ranges in the bed file are numeric
# therefore add a string suffix
names(dbaErEpPrEpUnOl$inAll) <- paste(names(dbaErEpPrEpUnOl$inAll), "occ_strict", sep = "_")

export.bed(dbaErEpPrEpUnOl$inAll, con = paste(outDir, "occupancy_strict", "occ_strict_ErEpPrEp.bed", sep = "/")) 
#of 59597 PR co-treatment gained peaks, and 19156 ER co-treatment gained peaks, 16985 peaks are overlapping, and thus of interest for further analysis*

### then check the overlap in the more lenient consensus peak set (peaks that are present in at least 2 out of 3 replicates)
# I do this according to this site:
# https://support.bioconductor.org/p/67306/

# retrieve peaksets of more lenient consensus peak sets (present in all three)
dbaErOl2Grl <- dba.overlap(dbaErOl2, c(7:8)) # computes overlapping and unioun sites of er_e and er_ep, and output it as a list of GRanges Objects
dbaPrOl2Grl <- dba.overlap(dbaPrOl2, c(7:8))# computes overlapping and unioun sites of pr_e and pr_ep, and output it as a list of GRanges Objects

# construct a new dba object containing consensus peaks of er_ep_unique and pr_ep_unique 
dbaErEpPrEpUnLe <- dba.peakset(NULL, dbaErOl2Grl$onlyB, sampID = "er_ep_u", tissue = "T47D", factor = "er", condition = "ep", treatment = "noTreat", replicate = "1-2-3", control = "noInp", peak.caller = "macs") # retrieve er_ep_unique peaks, and construct with this peak set a new DBA object

dbaErEpPrEpUnLe <- dba.peakset(dbaErEpPrEpUnLe, dbaPrOl2Grl$onlyB, sampID = "pr_ep_u", tissue = "T47D", factor = "pr", condition = "ep", treatment = "noTreat", replicate = "1-2-3", control = "noInp", peak.caller = "macs")

dba.plotVenn(dbaErEpPrEpUnLe, c(1,2)) # plot overlap of upon cotreatment gained pr and er binding sites

# extract overlapping peakset and export as bed.file to perform motif analysis on that peaks
dbaErEpPrEpUnLeOl <- dba.overlap(dbaErEpPrEpUnLe, c(1:2))
# findMotifsGenome.pl from Homer does not work if the name of the ranges in the bed file are numeric
# therefore add a string suffix
names(dbaErEpPrEpUnLeOl$inAll) <- paste(names(dbaErEpPrEpUnLeOl$inAll), "occ_len", sep = "_")

dir.create(paste(outDir, "occupancy_lenient", sep = "/"))
export.bed(dbaErEpPrEpUnLeOl$inAll, con = paste(outDir, "occupancy_lenient", "occ_lenient_ErEpPrEp.bed", sep = "/")) 
#of 83365 PR co-treatment gained peaks, and 33793 ER co-treatment gained peaks, 23754 peaks are overlapping, and thus of interest for further analysis*

## perform quantitative differential binding analysis:

### first identify which binding sites are significantly enriched by co-treatment

#### as a consensus peakset, use all peaks that are in at least 2 out of all 6 cell lines that where ChIPed with ER, or with PR

#### first do it for ER ChIPs:
dbaErOl2of6 <- dba.peakset(dbaEr, consensus = DBA_FACTOR, minOverlap = 2)
dba.show(dbaErOl2of6, mask = dbaErOl2of6$masks$Consensus) # show consensus peaksets

# convert consensus peakset in a GRanges Object file
dbaErOl2of6Gr <- dba.peakset(dbaErOl2of6, c(7), bRetrieve = TRUE)

# use this GRanges Ojbect, togehter with DBA object that contains single Bam and Peak files of all ER-ChIPed samples, 
# for counting of reads in this regions ----------

dbaErOl2Of6Qa <- dba.count(dbaEr, peaks = dbaErOl2of6Gr, score=DBA_SCORE_READS)
dbaErOl2Of6Qa

# plot correlation heatmap of ER-ChIP samples based on read counts in consensus peakset
plot(dbaErOl2Of6Qa)

# make contrasts, and analyse data --------------
dbaErOl2Of6Qa <- dba.contrast(dbaErOl2Of6Qa, categories = DBA_CONDITION)
dbaErOl2Of6Qa

dbaErOl2Of6Qa <- dba.analyze(dbaErOl2Of6Qa, method = c(DBA_EDGER, DBA_DESEQ2), bSubControl = TRUE, bFullLibrarySize = TRUE, bTagwise = TRUE)
dbaErOl2Of6Qa 
# 24443 sites found by edgeR, and 20177 by DESeq2

# extract GRanges Object with differntially enriched Peaks, Of interest are the peaks that are enriched in er_ep conditions

# # bDB=T generates a report-based DBA object; for subsequent steps, a GRranges Ojb is better
# dbaErOl2Of6QaDeUpRep <- dba.report(dbaErOl2Of6Qa, contrast = 1, method = DBA_DESEQ2, th = 0.05,
#                                 bUsePval = FALSE, fold = log2(2), DataType = DBA_DATA_GRANGES,
#                                 bNormalized = T, bDB=T, bGain = T, bCalled=T)
# dbaErOl2Of6QaDeUpRep$called # mit bCalled kann man nachschauen, ob der range auch als Peak gecalled wurde in der Gruppe1, also er_e, oder der Gruppe 2, also er_ep) 

dbaErOl2Of6QaDeUp <- dba.report(dbaErOl2Of6Qa, contrast = 1, method = DBA_DESEQ2, th = 0.05,
                                bUsePval = FALSE, fold = log2(2), DataType = DBA_DATA_GRANGES,
                                bNormalized = T)

# 20027 peaks are regulated with a fdr < 0.05 and a log2 fold change >1 or < 1

# extract only reads that are upregulated in co-treatment:
dbaErOl2Of6QaDeUp <- dbaErOl2Of6QaDeUp[mcols(dbaErOl2Of6QaDeUp)[,4] < -1]
# 18316 ranges are upregulated in co-treatment

dbaErOl2Of6QaEdgUp <- dba.report(dbaErOl2Of6Qa, contrast = 1, method = DBA_EDGER, th = 0.05,
                                bUsePval = FALSE, fold = log2(2), DataType = DBA_DATA_GRANGES,
                                bNormalized = T, bCalled = T)
# 23859 peaks are regulated under cotreatment with a FDR < 0.05 by EdgeR

# extract only reads that are upregulated in co-treatment:
dbaErOl2Of6QaEdgUp <- dbaErOl2Of6QaEdgUp[mcols(dbaErOl2Of6QaEdgUp)[,4] < -1]
# 14187 ranges are upregulated in co-treatment (edgeR. in contrast to DEseq, also calls a lot of downregulated binding sites)

# keep only peaks that are present in both
dbaErOl2Of6QaUp <- subsetByOverlaps(dbaErOl2Of6QaDeUp, dbaErOl2Of6QaEdgUp)

# 14183 ranges are detected by both Methods: So we have 14183 more enriched sites in er_ep than in er_e
#### then do it for PR ChIPs:
dbaPrOl2of6 <- dba.peakset(dbaPr, consensus = DBA_FACTOR, minOverlap = 2)
dba.show(dbaPrOl2of6, mask = dbaPrOl2of6$masks$Consensus) # show consensus peaksets

# convert consensus peakset in a GRanges Object file
dbaPrOl2of6Gr <- dba.peakset(dbaPrOl2of6, c(7), bRetrieve = TRUE)

# use this GRanges Ojbect, togehter with DBA object that contains single Bam and Peak files of all ER-ChIPed samples, 
# for counting of reads in this regions ----------

dbaPrOl2Of6Qa <- dba.count(dbaPr, peaks = dbaPrOl2of6Gr, score=DBA_SCORE_READS)
dbaPrOl2Of6Qa

# plot correlation heatmap of ER-ChIP samples based on read counts in consensus peakset
plot(dbaPrOl2Of6Qa)

# make contrasts, and analyse data --------------
dbaPrOl2Of6Qa <- dba.contrast(dbaPrOl2Of6Qa, categories = DBA_CONDITION)
dbaPrOl2Of6Qa

dbaPrOl2Of6Qa <- dba.analyze(dbaPrOl2Of6Qa, method = c(DBA_EDGER, DBA_DESEQ2), bSubControl = TRUE, bFullLibrarySize = TRUE, bTagwise = TRUE)
dbaPrOl2Of6Qa 
# 85280 sites found by edgeR, and 63969 by DESeq2

# extract GRanges Object with differntially enriched Peaks, Of interest are the peaks that are enriched in er_ep conditions
dbaPrOl2Of6QaDeUp <- dba.report(dbaPrOl2Of6Qa, contrast = 1, method = DBA_DESEQ2, th = 0.05,
                                bUsePval = FALSE, fold = log2(2), DataType = DBA_DATA_GRANGES,
                                bNormalized = T)

# 63966 peaks are regulated with a fdr < 0.05 and a log2 fold change >1 or < 1

# extract only reads that are upregulated in co-treatment:
dbaPrOl2Of6QaDeUp <- dbaPrOl2Of6QaDeUp[mcols(dbaPrOl2Of6QaDeUp)[,4] < -1]
# 63057 ranges are upregulated in co-treatment

dbaPrOl2Of6QaEdgUp <- dba.report(dbaPrOl2Of6Qa, contrast = 1, method = DBA_EDGER, th = 0.05,
                                bUsePval = FALSE, fold = log2(2), DataType = DBA_DATA_GRANGES,
                                bNormalized = T, bCalled = T)
# 84120 peaks are regulated under cotreatment with a FDR < 0.05 by EdgeR

# extract only reads that are upregulated in co-treatment:
dbaPrOl2Of6QaEdgUp <- dbaPrOl2Of6QaEdgUp[mcols(dbaPrOl2Of6QaEdgUp)[,4] < -1]
# 41362 ranges are upregulated in co-treatment (edgeR. in contrast to DEseq, also calls a lot of downregulated binding sites)

# keep only peaks that are present in both
dbaPrOl2Of6QaUp <- subsetByOverlaps(dbaPrOl2Of6QaDeUp, dbaPrOl2Of6QaEdgUp)

# 41363 ranges are detected by both Methods: So we have 41363 more enriched sites in pr_ep than in pr_e
### overlap ER and PR co-treatment significantly higher enriched regions, and write into file for motif analysis 
dbaErPrOl2Of6QaUp <- dba.peakset(NULL, dbaErOl2Of6QaUp, sampID = "er_ep_up", tissue = "T47D", factor = "er", condition = "ep", treatment = "noTreat", replicate = "1-2-3", control = "noInp", peak.caller = "macs")

dbaErPrOl2Of6QaUp <- dba.peakset(dbaErPrOl2Of6QaUp, dbaPrOl2Of6QaUp, sampID = "pr_ep_up", tissue = "T47D", factor = "pr", condition = "ep", treatment = "noTreat", replicate = "1-2-3", control = "noInp", peak.caller = "macs")

dba.plotVenn(dbaErPrOl2Of6QaUp, c(1,2)) # plot overlap of upon cotreatment gained pr and er binding sites
# all but 681 of the er_ep aquired peaks, so 13308, overlap with the pr_ep_up regulated genes

# extract overlapping peakset and export as bed.file to perform motif analysis on that peaks
dbaErPrOl2Of6QaUpOl <- dba.overlap(dbaErPrOl2Of6QaUp, c(1:2))
dir.create(paste(outDir, "quantitative", sep = "/"))
# findMotifsGenome.pl from Homer does not work if the name of the ranges in the bed file are numeric
# therefore add a string suffix
names(dbaErPrOl2Of6QaUpOl$inAll) <- paste(names(dbaErPrOl2Of6QaUpOl$inAll), "quant", sep = "_")
export.bed(dbaErPrOl2Of6QaUpOl$inAll, con = paste(outDir, "quantitative", "quantitativeErEpPrEpUp.bed", sep = "/")) 
###finally see how the overlap is in the peaks obtained from occupancy, and von quantitative, analysis
olOccQuanErEpPrEpUp <- dba.peakset(NULL, dbaErEpPrEpUnOl$inAll, sampID = "occ_strict_ErPrUpOl", tissue = "T47D", factor = "ErPr", condition = "ep", treatment = "noTreat", replicate = "1-2-3", control = "noInp", peak.caller = "macs")

olOccQuanErEpPrEpUp <- dba.peakset(olOccQuanErEpPrEpUp, dbaErEpPrEpUnLeOl$inAll, sampID = "occ_lenient_ErPrUpOl", tissue = "T47D", factor = "ErPr", condition = "ep", treatment = "noTreat", replicate = "1-2-3", control = "noInp", peak.caller = "macs")

olOccQuanErEpPrEpUp <- dba.peakset(olOccQuanErEpPrEpUp, dbaErPrOl2Of6QaUpOl$inAll, sampID = "quant_ErPrUpOl", tissue = "T47D", factor = "ErPr", condition = "ep", treatment = "noTreat", replicate = "1-2-3", control = "noInp", peak.caller = "macs")

dba.plotVenn(olOccQuanErEpPrEpUp, c(1:3))
dba.plotVenn(olOccQuanErEpPrEpUp, c(1,3))
dba.plotVenn(olOccQuanErEpPrEpUp, c(2:3))
#9702 sites are reported by all three methods; 1115 sites are only reported by the quantitative affinity method*

#### scatterplot of pvalue vs. fold change to see whether there are many peaks with a small p value and a very small log2 fold change 
# for that, we have to extract the ER-PR-ep upregulated peaks from all analyzed ER_ep and PR_ep reports, since after the overlap the metadata
# columns get deleted
#### do great analysis
job = submitGreatJob(dbaErPrOl2Of6QaUpOl$inAll, species = "hg19")
tb = getEnrichmentTables(job)
names(tb)
#### save.image
save.image(paste(outDir, "diffBind.RData", sep = "/"))
