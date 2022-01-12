# ---
# title: "diffBind"
# author: "ReinhardGrausenburger"
# date: "March 3, 2017"
# output: html_document
# ---


### load required packages

#source("https://bioconductor.org/biocLite.R")
#biocLite("ChIPQC")
library(ChIPQC)
browseVignettes("DiffBind")
library(DiffBind)
library(GenomicRanges)
library(rGREAT)
### define working directory and directory for output

getwd()
wdir <- "/Volumes/Temp/PBI_ChipSeq/data" 
setwd(wdir)

outDir <- paste(wdir, "12_diffBind", sep = "/")
dir.create(outDir)



# Overlap analysis of replicates from indivudual conditions

#### Define dba object
dbaO <- dba(sampleSheet = "~/samples.black.csv")
dbaO
#### Plot correlation of -log10 p Values of Peaks to show how samples cluster
plot(dbaO)
#*replicates cluster nicely; again, conditions are closer than the ChIPs with identical antibodies* \

#### do analysis with a more lenient definition for consensus Peaksets
#er_e und er_ep
# this is another, more intutive way of makeing consensus peaksets
dbaEr <- dba(dbaO, c(1:6)) # only select ERalpha ChIPs
dbaErOl2 <- dba.peakset(dbaEr, consensus = DBA_CONDITION, minOverlap = 2)
dba.show(dbaErOl2, mask = dbaErOl2$masks$Consensus)

# plot consensus peaksets
dba.plotVenn(dbaErOl2, mask=dbaErOl2$masks$Consensus)
#it looks as there 11507 sites lost, and 34302 aquired upon co-treatment, so here, the difference between lost and gained binding sites is not so dramatic*

## now find PR acquired binding sites during estrogen/progesteron co-treatment
#### do analysis with a more lenient definition for consensus Peaksets
dbaPr <- dba(dbaO, c(7:12)) # only select PR ChIPs
dbaPrOl2 <- dba.peakset(dbaPr, consensus = DBA_CONDITION, minOverlap = 2)
dba.show(dbaPrOl2, mask = dbaPrOl2$masks$Consensus)
# plot consensus peaksets
dba.plotVenn(dbaPrOl2, mask=dbaPrOl2$masks$Consensus)
#it looks as there 25685 sites lost, and 83578 sites aquired upon co-treatment, so here, the difference between lost and gained binding sites is not so dramatic as with a more strict consensus peak set*

## now check how big the overlap between ERalpha and PR binding sites upon estrogen/progesterone co-treatment is 

### first check the overlap in the more strict consensus peak set (intersection)
# I do this according to this site:
# https://support.bioconductor.org/p/67306/

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

# use this GRanges Ojbect, together with DBA object that contains single Bam and Peak files of all ER-ChIPed samples, 
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
dbaErOl2Of6QaDeUp
# 20027 peaks are regulated with a fdr < 0.05 and a log2 fold change >1 or < 1

# extract only reads that are upregulated in co-treatment:
dbaErOl2Of6QaDeUp <- dbaErOl2Of6QaDeUp[mcols(dbaErOl2Of6QaDeUp)[,4] < -1]
dbaErOl2Of6QaDeUp
# 18316 ranges are upregulated in co-treatment

dbaErOl2Of6QaEdgUp <- dba.report(dbaErOl2Of6Qa, contrast = 1, method = DBA_EDGER, th = 0.05,
                                bUsePval = FALSE, fold = log2(2), DataType = DBA_DATA_GRANGES,
                                bNormalized = T, bCalled = T)
# 23859 peaks are regulated under cotreatment with a FDR < 0.05 by EdgeR

# extract only reads that are upregulated in co-treatment:
dbaErOl2Of6QaEdgUp <- dbaErOl2Of6QaEdgUp[mcols(dbaErOl2Of6QaEdgUp)[,4] < -1]
dbaErOl2Of6QaEdgUp
# 14187 ranges are upregulated in co-treatment (edgeR. in contrast to DEseq, also calls a lot of downregulated binding sites)

# keep only peaks that are present in both
dbaErOl2Of6QaUp <- subsetByOverlaps(dbaErOl2Of6QaDeUp, dbaErOl2Of6QaEdgUp)
dbaErOl2Of6QaUp

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
dbaPrOl2Of6QaDeUp

# 63966 peaks are regulated with a fdr < 0.05 and a log2 fold change >1 or < 1

# extract only reads that are upregulated in co-treatment:
dbaPrOl2Of6QaDeUp <- dbaPrOl2Of6QaDeUp[mcols(dbaPrOl2Of6QaDeUp)[,4] < -1]
dbaPrOl2Of6QaDeUp

# 63057 ranges are upregulated in co-treatment

dbaPrOl2Of6QaEdgUp <- dba.report(dbaPrOl2Of6Qa, contrast = 1, method = DBA_EDGER, th = 0.05,
                                bUsePval = FALSE, fold = log2(2), DataType = DBA_DATA_GRANGES,
                                bNormalized = T, bCalled = T)
dbaPrOl2Of6QaEdgUp

# 84120 peaks are regulated under cotreatment with a FDR < 0.05 by EdgeR

# extract only reads that are upregulated in co-treatment:
dbaPrOl2Of6QaEdgUp <- dbaPrOl2Of6QaEdgUp[mcols(dbaPrOl2Of6QaEdgUp)[,4] < -1]
dbaPrOl2Of6QaEdgUp
# 41362 ranges are upregulated in co-treatment (edgeR. in contrast to DEseq, also calls a lot of downregulated binding sites)

# keep only peaks that are present in both
dbaPrOl2Of6QaUp <- subsetByOverlaps(dbaPrOl2Of6QaDeUp, dbaPrOl2Of6QaEdgUp)
dbaPrOl2Of6QaUp

# 41363 ranges are detected by both Methods: So we have 41363 more enriched sites in pr_ep than in pr_e
### overlap ER and PR co-treatment significantly higher enriched regions, and write into file for motif analysis 
dbaErPrOl2Of6QaUp <- dba.peakset(NULL, dbaErOl2Of6QaUp, sampID = "er_ep_up", tissue = "T47D", factor = "er", condition = "ep", treatment = "noTreat", replicate = "1-2-3", control = "noInp", peak.caller = "macs")

dbaErPrOl2Of6QaUp <- dba.peakset(dbaErPrOl2Of6QaUp, dbaPrOl2Of6QaUp, sampID = "pr_ep_up", tissue = "T47D", factor = "pr", condition = "ep", treatment = "noTreat", replicate = "1-2-3", control = "noInp", peak.caller = "macs")
dbaErPrOl2Of6QaUp

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

#dbaErEpPrEpU <- dba.peakset(NULL, dbaErOl3Grl$onlyB, sampID = "er_ep_u", tissue = "T47D", factor = "er", condition = "ep", treatment = "noTreat", replicate = "1-2-3", control = "noInp", peak.caller = "macs") # retrieve er_ep_unique peaks, and construct with this peak set a new DBA object

#dbaErEpPrEpU <- dba.peakset(dbaErEpPrEpU, dbaPrOl3Grl$onlyB, sampID = "pr_ep_u", tissue = "T47D", factor = "pr", condition = "ep", treatment = "noTreat", replicate = "1-2-3", control = "noInp", peak.caller = "macs")

#dba.plotVenn(dbaErEpPrEpU, c(1,2)) # plot overlap of upon cotreatment gained pr and er binding sites

dbaErEpPrEpUnOl <- dba.overlap(dbaErEpPrEpUnLe, c(1:2))

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
