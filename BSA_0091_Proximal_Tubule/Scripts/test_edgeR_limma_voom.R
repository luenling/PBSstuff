setwd("/Volumes/Temp/Lukas/Data/Test/")
#setwd("/Volumes/vetgrid01/Data/Test/")
library(limma)
library(edgeR)

#filename="/Volumes/Temp/Lukas/BSA_0091_Proximal_Tubule/mm10//featCounts/mm10.1000nt.featurecount.counts_edgeR.table"
filename="/Volumes/Temp/Lukas/BSA_0091_Proximal_Tubule/mm10//featCounts/new.test/ucsc_mm10_ensembl_exons_1000nt_most_highly.counts_edgeR.table"
dropcols=c("WT_KO4_8h_DT","Klotho_VDR_KO2_2h_PT","Klotho_VDR_KO2_2h_DT","Klotho_VDR_KO4_8h_PT","Klotho_VDR_KO4_8h_DT","Klotho_VDR_KO3_8h_DT","VDR_KO4_8h_PT")
drop_reg="^VDR_KO"
a = create.dge(filename,dropcols=dropcols,drop_reg=drop_reg,both=T)
a[["dge"]] <- calcNormFactors(a[["dge"]])
a[["dge"]]$samples
plotMDS(a[["dge"]])

dges.1000nt = create.dge(filename,dropcols=dropcols,drop_reg=drop_reg,both=F)
# get sample matrices
sample.mat.PT=dges.1000nt[["sample.mat.PT"]]
sample.mat.DT=dges.1000nt[["sample.mat.DT"]]
# create design matrices
design.PT=model.matrix( ~ sample.mat.PT$GT + sample.mat.PT$time + sample.mat.PT$GT:sample.mat.PT$time )
rownames(design.PT)=rownames(sample.mat.PT)
design.DT=model.matrix( ~ sample.mat.DT$GT + sample.mat.DT$time + sample.mat.DT$GT:sample.mat.DT$time )
rownames(design.DT)=rownames(sample.mat.DT)
# set nicer variables for dge objects
dge.PT=dges.1000nt[["dge.PT"]] 
dge.DT=dges.1000nt[["dge.DT"]]
# if you do not want a safety copy of the list rm it:
#rm(dges.1000nt)


# different dispersion estimates
dge.PT<- estimateGLMCommonDisp(dge.PT, design.PT, verbose=TRUE)
dge.PT<- estimateGLMTrendedDisp(dge.PT, design.PT, verbose=F)

dge.PT<- estimateGLMTagwiseDisp(dge.PT, design.PT)
plotBCV(dge.PT)
plotMeanVar(dge.PT,show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)
# from edgeR tutorial http://www.bioconductor.org/help/course-materials/2014/BioC2014/BioC2014_edgeR_voom.html
#par(mfrow=c(1,2))
plotSmear(dge.PT, pair=c("FGF23_VDR_0h_PT","FGF23_VDR_2h_PT"), ylim=c(-5,5))
abline(h=0,col="red")
f <- glmFit(dge.PT,design.PT)
con <- makeContrasts("DE-ES"=grpDE-grpES,levels=colnames(mm))
con
# fit big design matrix
fitdges.1000nt <- glmFit(dge.PT, design.PT)
lrt <- glmLRT(fitdges.1000nt,coef=3)
topTags(lrt)
# plot the cpms as barplots
cps <- cpm(dge.PT)
cols<-as.numeric(dge.PT$samples$group)
o <- order(dge.PT$samples$group)
barplot( cps["ENSMUSG00000006724",o], col=cols[o], las=2,main=dge.PT$genes[c("ENSMUSG00000006724"),c("gene_symbol")])
# for limma/voom
# from limma user guide
dge=dges.1000nt[["dge.PT"]] 
dge <- calcNormFactors(dge)
v <- voom(dge,design.PT,plot=TRUE)
# without TMM normalisation
v <- voom(dge$counts,design.PT,plot=TRUE)
# voom with quality weights
cols=as.numeric(dge$samples$group)
vQW = voomWithQualityWeights(dge,design=design.PT,plot=T,col=cols)
# from http://bioinf.wehi.edu.au/RNAseqCaseStudy/
#just if you want rpkms:
dge_rpkm <- rpkm(dge,dge$genes$Length)
# no, no - this is wrong, voom wants raw counts , while rpkm looks actually better it is not waht is meant to be
#v <- voom(dge_rpkm,design.PT,plot=TRUE)
# get rid of lowly expressed stuff, not sure whether this is necessary, we already removed all < min_count
# cpm = 0.5 cpm is still a high threshold with 30 mio fragments per sample and only counting around 40-50% of reads (around 7.5 counts) 
isexpr <- rowSums(cpm(dge) > 0.15) >= 15
y <- voom(dge[isexpr,],design.PT,plot=TRUE)
plotMDS(dge[isexpr,],xlim=c(-2.5,2.5))
fit <- eBayes(lmFit(v,design.PT))
topTable(fit,coef=3)
# remove batcheffect: no effect?
v2 <- v
v2$E <- removeBatchEffect(v,design=design.PT)
fitQW <- eBayes(lmFit(vQW,design.PT))
topTable(fitQW,coef=3)
decTQW<-decideTests(fitQW,method = "separate")
colnames(decTQW)
cols <- as.numeric(v$targets$group)
plotMDS(v,col=cols)
hist(log(cpm(dge)))

plotMDS(v,col=cols, top=2000, main="2000 / lLFC", gene.selection = "common")




create.dge <- function(filename,mincount=15,dropcols=c("VDR_KO4_8h_PT"),both=FALSE,drop_reg=FALSE){
  # gets a filename, the minimal count per gene (over all samples in the returned count matrices), lists of columns to drop, a regular expression for columns to drop and the boolean "both", that it true creates a edger matrix for DT and one for PT
  # returns a list with a dge object and a sampl.mat dataframe for either all samples or for PT and DT each
  mm10=read.table(filename,header=T,row.names="Geneid")
  # drop with regular expression
  if (drop_reg != FALSE) {    
    dropcols=c(dropcols,names(mm10)[grep(drop_reg,names(mm10),perl=TRUE)])
  }
  # drop samples from list
  mm10=mm10[,!(names(mm10) %in% dropcols)]
  # drop rows with counts less than mincount
  mm10=mm10[!(rowSums(mm10[,5:ncol(mm10)]) < mincount), ]
  # create sample matrix
  count.cols=5:ncol(mm10)
  samples=colnames(mm10)[count.cols]
  sample.mat=data.frame(GT=sub('^(\\w*VDR|WT).*','\\1',colnames(mm10)[count.cols],perl=TRUE),
                       samples=samples,
                       loc=sub('.*_([PDT]+)$','\\1',samples,perl=TRUE),
                       group=sub('^(\\w*VDR|WT)(?:_KO)?\\d_','\\1_',samples,perl=TRUE),
                       time=sub('.*_(\\d+)h_.*','\\1',samples,perl=TRUE))
  rownames(sample.mat)=samples
  sample.mat$GT=relevel(sample.mat$GT,"WT")
  sample.mat$time=relevel(sample.mat$time,"0")
  sample.mat$loc=relevel(sample.mat$loc,"DT")
  # obtain gene identifiers and lengths
  genes=mm10[,1:4]
  genes$gene_symbol <- convertIDs(row.names(genes), "ENSEMBL", "SYMBOL", org.Mm.eg.db)
  genes$entrezgene <- convertIDs(row.names(genes), "ENSEMBL", "ENTREZID", org.Mm.eg.db)
  # take the counts part of the df
  mm10_counts=mm10[,count.cols]  
  if (both){
    # create dge object to return
    dge <- DGEList(counts=mm10_counts, group=sample.mat$group, genes=genes)
    return(list( "dge" = dge, "sample.mat"  = sample.mat ))
  }
  else{
    # samples that are PT
    keep_PT <- sample.mat$loc == "PT"
    # create the reduced sample matrices
    sample.mat.PT = droplevels(sample.mat[keep_PT,])
    sample.mat.DT = droplevels(sample.mat[ ! keep_PT,])
    # drop rows with counts less than mincount for PT and DT
    rkeep_PT = ! (rowSums(mm10_counts[,keep_PT]) < mincount)
    rkeep_DT = ! (rowSums(mm10_counts[, ! keep_PT]) < mincount)
    mm10_counts_PT=mm10_counts[ rkeep_PT , keep_PT]
    mm10_counts_DT=mm10_counts[ rkeep_DT , ! keep_PT]
    # create dge objects
    dge_DT <- DGEList(counts=mm10_counts_DT, group=sample.mat.DT$group, genes=genes[rkeep_DT,])
    dge_PT <- DGEList(counts=mm10_counts_PT, group=sample.mat.PT$group, genes=genes[rkeep_PT,])
    return(list( "dge.PT" = dge_PT, "sample.mat.PT" = sample.mat.PT ,"dge.DT" = dge_DT,"sample.mat.DT" = sample.mat.DT ))  
  }  
}

library("AnnotationDbi")
library("org.Mm.eg.db")
# use the convertID function from http://www.bioconductor.org/help/workflows/rnaseqGene/
convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}


##### Analysis Limma

### Read files and create object for DESEQ and LIMMA, normalise voom
### PT and DT together for voom normalisation

filename="/Volumes/Temp/Lukas/BSA_0091_Proximal_Tubule/mm10//featCounts/new.test/ucsc_mm10_ensembl_exons_1000nt_most_highly.counts_edgeR.table"
dropcols=c("WT_KO4_8h_DT","Klotho_VDR_KO2_2h_PT","Klotho_VDR_KO2_2h_DT","Klotho_VDR_KO4_8h_PT","Klotho_VDR_KO4_8h_DT","Klotho_VDR_KO3_8h_DT","VDR_KO4_8h_PT")
drop_reg="^VDR_KO"
all.samps = create.dge(filename,dropcols=dropcols,drop_reg=drop_reg,both=T)
design.mt =model.matrix( ~ GT + time + loc, data=all.samps$sample.mat ) 
dupcor1=duplicateCorrelation(all.samps$dge[[1]],design=model.matrix( ~ all.samps$sample.mat$group))
dupcor2=duplicateCorrelation(all.samps$dge[[1]],design=design.mt)
# calculate the TMM normalisation factors for library sizes
all.samps$dge=calcNormFactors(all.samps$dge)
a1.test = voom(all.samps$dge, plot=TRUE)
a2.test = voom(all.samps$dge, design= design.mt, plot=TRUE)

## create list of subobjects with only relevant samples for GT and time
GT.times.dges = get.GTandTime.dges(all.samps$dge,all.samps$sample.mat,mincount=10)
 

get.GTandTime.dges <- function(dge,sample.mat,mincount=10){
  ## takes a dge and a sample matrix and returns a list of four lists containing the dges,sample and design matrices for all genotypes or time-points separated for PT and DT 
  ## the dges are already normalised and have voom w. qw run
  ## the voom plots are stored as pdfs
  dge.GTs=list()
  dge.times=list()
  sm.GTs=list()
  sm.times=list()
  design.GTs=list()
  design.times=list()
  
  # get locations
  for (j in levels(sample.mat$loc)){
    dge.GTs[[j]]=list()
    sm.GTs[[j]]=list()
    design.GTs[[j]]=list()
    # get the dges for genotypes
    for (i in levels(sample.mat$GT)){
      keep.GT <- sample.mat$loc == j & sample.mat$GT == i
      sample.mat.GT = droplevels(sample.mat[keep.GT,])
      rkeep.GT = ! (rowSums(dge[[1]][,keep.GT]) < mincount)
      dge.GTs[[j]][[i]] <- DGEList(counts=dge[[1]][rkeep.GT,keep.GT], group=sample.mat.GT$group, genes=dge[[3]][rkeep.GT,])
      sm.GTs[[j]][[i]] <- sample.mat.GT
      design.GTs[[j]][[i]] <- model.matrix( ~ time,data=sample.mat.GT)
      # normalisation factors
      dge.GTs[[j]][[i]] <- calcNormFactors(dge.GTs[[j]][[i]])
       filename=paste("voom_qw",j,"GT",i,"pdf",sep = ".")
      pdf(file=filename)
      # calculate with quality weights
      dge.GTs[[j]][[i]] <- voomWithQualityWeights(dge.GTs[[j]][[i]], design=design.GTs[[j]][[i]], normalization="none", plot=TRUE)
      dev.off()
    }
    dge.times[[j]]=list()
    sm.times[[j]]=list()
    design.times[[j]]=list()
    # get the dges for times
    for (i in levels(sample.mat$time)){
      keep.time <- sample.mat$loc == j & sample.mat$time == i
      sample.mat.time = droplevels(sample.mat[keep.time,])
      rkeep.time = ! (rowSums(dge[[1]][,keep.time]) < mincount)
      dge.times[[j]][[i]] <- DGEList(counts=dge[[1]][rkeep.time,keep.time], group=sample.mat.time$group, genes=dge[[3]][rkeep.time,])
      sm.times[[j]][[i]] <- sample.mat.time
      design.times[[j]][[i]] <- model.matrix( ~ GT, data = sample.mat.time)
      # normalisation factors
      dge.times[[j]][[i]] <- calcNormFactors(dge.times[[j]][[i]])
      filename=paste("voom_qw",j,"time",i,"pdf",sep = ".")
      pdf(file=filename)
      # calculate with quality weights
      dge.times[[j]][[i]] <- voomWithQualityWeights(dge.times[[j]][[i]], design=design.times[[j]][[i]], normalization="none", plot=TRUE)
      dev.off()
    }  
  }
  return(list( "dge.GTs" = dge.GTs, "sm.GTs" = sm.GTs , "design.GTs" = design.GTs, "dge.times" = dge.times, "sm.times" = sm.times, "design.times" = design.times ))
}

for(loc in names(GT.times.dges[["dge.GTs"]]) ){
  for(gt in names(GT.times.dges[["dge.GTs"]][["DT"]]) ) {
    sample.mat.gt <- GT.times.dges[["sm.GTs"]][[loc]][[gt]]
    GT.times.dges[["design.GTs"]][[loc]][[gt]] <- model.matrix( ~ time, data = sample.mat.gt)
  }
  for(time in names(GT.times.dges[["dge.times"]][["DT"]]) ) {
    sample.mat.times <- GT.times.dges[["sm.times"]][[loc]][[time]]
    GT.times.dges[["design.times"]][[loc]][[time]] <- model.matrix( ~ GT, data = sample.mat.times)
  }
}

# define contrasts
cont.wt <- makeContrasts(t0t2="-time2",t0t8="-time8",
                         t2t8="time2-time8",levels=GT.times.dges[["design.GTs"]][["PT"]][["WT"]])
colnames(GT.times.dges[["design.GTs"]][["PT"]][["WT"]])=c("Intercept","time2","time8")
fit <- lmFit(GT.times.dges[["dge.GTs"]][["PT"]][["WT"]], design=GT.times.dges[["design.GTs"]][["PT"]][["WT"]])
fit2 <- contrasts.fit(fit, cont.wt)
fit2 <- eBayes(fit2)
head(fit2)
toptable(fit2,coef="t0t8",genelist = fit2$genes[,c(5)])

str(GT.times.dges[["dge.GTs"]][["PT"]][["WT"]])
str(GT.times.dges[["design.GTs"]][["PT"]][["WT"]])

# calculate model


# get results





