#setwd("/Volumes/Temp/Lukas/Data/Test/")
setwd("/Volumes/vetgrid01/Data/Test/")
library(limma)
library(edgeR)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
require(grDevices)



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
  rkeep = ! ( apply( mm10[,5:ncol(mm10)], 1, max, na.rm=TRUE) < mincount )
  mm10=mm10[ rkeep , ]
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
    rkeep_PT = ! ( apply( mm10_counts[,keep_PT], 1, max, na.rm=TRUE) < mincount ) 
    rkeep_DT = ! ( apply( mm10_counts[, ! keep_PT], 1, max, na.rm=TRUE)  < mincount)
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
      rkeep.GT = ! ( apply(dge[[1]][,keep.GT],1,max,na.rm = TRUE) < mincount)
      #rkeep.GT = ! (rowSums(dge[[1]][,keep.GT]) < mincount)
      dge.GTs[[j]][[i]] <- DGEList(counts=dge[[1]][rkeep.GT,keep.GT], group=sample.mat.GT$group, genes=dge[[3]][rkeep.GT,])
      sm.GTs[[j]][[i]] <- sample.mat.GT
      design.GTs[[j]][[i]] <- model.matrix( ~ 0 + time,data=sample.mat.GT)
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
      design.times[[j]][[i]] <- model.matrix( ~ 0 + GT, data = sample.mat.time)
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

plot.Evals <- function(gene,dge.sub=GT.times.dges[["dge.GTs"]][["PT"]][["WT"]], ... ) {
  # plot the different replicates separately for a given gene
  # get sample groups from dge design matrix
  if (class(dge.sub) == "DGEList") {
    samp.groups=as.data.frame(rownames(dge.sub$samples))
    samp.groups=cbind(samp.groups,dge.sub$samples$group)
    colnames(samp.groups)=c("value","L1")
    evals=melt(dge.sub$counts[gene,])
  } else {
    samp.groups=melt(lapply(as.data.frame(dge.sub$design),function(x) rownames(dge.sub$design[x == 1,])))
    evals=melt(dge.sub$E[gene,])    
  }
  evals=cbind(evals,as.factor(samp.groups$L1[sapply(rownames(evals),function (x) which(x == samp.groups$value))]))
  colnames(evals) = c("E","samples")
  qplot(samples, E, data=evals, main=gene, ...) + theme_bw() 
  # or for boxplots
  # plot(E ~ samples, data=evals)  
}

vulcano.plot <- function(fit2, contrast, method="BH", pch=16,cex=0.45, cols=c("red","green","orange"), ...) {
  plot(fit2[["coefficients"]][,contrast],-log10(fit2[["p.value"]][,contrast]), pch=pch,cex=cex, ...)
  sig= p.adjust(fit2[["p.value"]][,contrast],method=method) <= 0.05
  lfgt2= abs(fit2[["coefficients"]][,contrast]) > 1.0
  points(fit2[["coefficients"]][sig & ! lfgt2,contrast],-log10(fit2[["p.value"]][sig & ! lfgt2,contrast]), pch=pch,cex=cex,col=cols[1])
  points(fit2[["coefficients"]][lfgt2 & ! sig,contrast],-log10(fit2[["p.value"]][lfgt2 & ! sig,contrast]), pch=pch,cex=cex,col=cols[2])
  points(fit2[["coefficients"]][sig & lfgt2,contrast],-log10(fit2[["p.value"]][sig & lfgt2,contrast]), pch=pch,cex=cex,col=cols[3])
}

##### Analysis Limma

### Read files and create object for EdgeR and LIMMA, normalise voom
### PT and DT together for voom normalisation

filename="/Volumes/vetgrid01/BSA_0091_Proximal_Tubule/mm10//featCounts/new.test/ucsc_mm10_ensembl_exons_1000nt_most_highly.counts_edgeR.table"
dropcols=c("WT_KO4_8h_DT","Klotho_VDR_KO2_2h_PT","Klotho_VDR_KO2_2h_DT","Klotho_VDR_KO4_8h_PT","Klotho_VDR_KO4_8h_DT","Klotho_VDR_KO3_8h_DT","VDR_KO4_8h_PT")
drop_reg="^VDR_KO"
all.samps = create.dge(filename,dropcols=dropcols,drop_reg=drop_reg,both=T,mincount=15)
design.mt =model.matrix( ~ 0 + GT + time + loc, data=all.samps$sample.mat ) 
#dupcor1=duplicateCorrelation(all.samps$dge[[1]],design=model.matrix( ~ 0 + all.samps$sample.mat$group))
#dupcor2=duplicateCorrelation(all.samps$dge[[1]],design=design.mt)
# calculate the TMM normalisation factors for library sizes
all.samps$dge=calcNormFactors(all.samps$dge)
#a1.test = voom(all.samps$dge, plot=TRUE)
#a2.test = voom(all.samps$dge, design= design.mt, plot=TRUE)

library(reshape2)
df = melt(as.data.frame(log2(all.samps$dge[[1]]+1)), variable.name = "Samples")
boxplot(df$value ~ df$Samples)
# all look very similar - median at around 2^5=32
## create list of subobjects with only relevant samples for GT and time
GT.times.dges = get.GTandTime.dges(all.samps$dge,all.samps$sample.mat,mincount=25)
 


# define the model matrices for the different GTs and times (already correctly done in the function get.GTandTime.dges)
for(loc in names(GT.times.dges[["dge.GTs"]]) ){
  for(gt in names(GT.times.dges[["dge.GTs"]][["DT"]]) ) {
    sample.mat.gt <- GT.times.dges[["sm.GTs"]][[loc]][[gt]]
    GT.times.dges[["design.GTs"]][[loc]][[gt]] <- model.matrix( ~  0 + time, data = sample.mat.gt)
  }
  for(time in names(GT.times.dges[["dge.times"]][["DT"]]) ) {
    sample.mat.times <- GT.times.dges[["sm.times"]][[loc]][[time]]
    GT.times.dges[["design.times"]][[loc]][[time]] <- model.matrix( ~ 0 + GT, data = sample.mat.times)
  }
}

#save.image("wed30092015.RData")
#load("wed30092015.RData")
# loop over all GTs

fits=list(list(list()))
for (loc in names(GT.times.dges[["dge.GTs"]]) ){
  for (gt in names(GT.times.dges[["dge.GTs"]][["DT"]]) ) {
    conts <- makeContrasts(t0t2="time0-time2",t0t8="time0-time8",
                             t2t8="time2-time8",levels=GT.times.dges[["design.GTs"]][[loc]][[gt]])
    fit <- lmFit(GT.times.dges[["dge.GTs"]][[loc]][[gt]], design=GT.times.dges[["design.GTs"]][[loc]][[gt]])
    fit <- contrasts.fit(fit, conts)
    fit <- eBayes(fit)
    fits[["GTs"]][[loc]][[gt]] <- fit
  }
}
for (loc in names(GT.times.dges[["dge.times"]]) ){
  for (time in names(GT.times.dges[["dge.times"]][["DT"]]) ) {
    conts <- makeContrasts(WT_FGF23="GTWT-GTFGF23_VDR",WT_Klotho="GTWT-GTKlotho_VDR",
                           FGF23_Klotho="GTFGF23_VDR-GTKlotho_VDR",levels=GT.times.dges[["design.times"]][[loc]][[time]])
    fit <- lmFit(GT.times.dges[["dge.times"]][[loc]][[time]], design=GT.times.dges[["design.times"]][[loc]][[time]])
    fit <- contrasts.fit(fit, conts)
    fit <- eBayes(fit)
    fits[["times"]][[loc]][[time]] <- fit
  }
}

for (loc in names(GT.times.dges[["dge.times"]]) ){
  for (gt in names(GT.times.dges[["dge.GTs"]][["DT"]])){
    outspec=paste0("location_",loc,"_genotype_",gt)
    fit=fits[["GTs"]][[loc]][[gt]]
    conts=colnames(fit$contrasts)
    # create directory for saving files
    dir.create(gt,showWarnings ='F')    
    # get toptable for each contrast, merge and write to file
    tts=list()
    for (cont in conts){
      tts[[cont]]= topTable(fit,coef=cont,adjust.method="BH",number=nrow(fit$genes),sort.by="none")
      names(tts[[cont]])[7:12]=paste(names(tts[[cont]])[7:12],cont,sep="_")
    }
    ttmerge=merge(tts[[1]],tts[[2]],by=c(1:6))    
    ttmerge=merge(ttmerge,tts[[3]],by=c(1:6))
    # get decide tests over contrasts
    decT=decideTests(fit,adjust.method="BH",method="separate")
    ttmerge=cbind(ttmerge,decT)
    rownames(ttmerge)=rownames(tts[[1]])
    ttmerge=cbind(rownames(ttmerge),ttmerge)
    colnames(ttmerge)[1]="ENSEMBLEID"
    # write merged table
    write.table(ttmerge,file=paste0(gt,"/results_",outspec,"_all_times.tab"),row.names=F,sep=",")
    # plot VennDiagram for decide tests
    outfile=paste0(gt,"/sig_de_",outspec)
    pdf(file=paste0(outfile,"_updown.pdf"))
    vennDiagram(decT,include="both", main=paste(loc,gt,"diff. expressed"))
    dev.off()
    pdf(file=paste0(outfile,"_up.pdf"))
    vennDiagram(decT,include="up", main=paste(loc,gt,"up regulated"))
    dev.off()
    pdf(file=paste0(outfile,"_down.pdf"))
    vennDiagram(decT,include="down", main=paste(loc,gt,"down regulated"))
    dev.off()
    pdf(file=paste0(outfile,"_heatdiagram.pdf"))
    heatDiagram(decT,fit$coef,names=rownames(fit$genes))
    dev.off()    
    for (cont in conts){
      # plot MA plot for each contrast
      pdf(file=paste0(gt,"/MA_",outspec,cont,".pdf"))
      adj.pvals=p.adjust(fit$p.value[,cont],method="BH")
      limma::plotMA(fit,coef=cont,status=as.factor(c("not sig.","sig.")[(adj.pvals <= 0.05) + 1]) ,col=c("red"),cex=0.3,main= paste(loc,gt,cont))
      dev.off()
      # volcano plot fo each contrast
      pdf(file=paste0(gt,"/vulcano_",outspec,cont,".pdf"))
      vulcano.plot(fit,cont,main= paste(loc,gt,cont),ylab=expression(-plain(log)[10](italic(P))), xlab=expression(plain(log)[2](italic(FC))) )
      dev.off()
    }
    # plot geneas for each contrast combination
    genas.cont=matrix(c(c(1,2),c(1,3),c(2,3)),nrow=2,ncol=3)
    genas.results=c()
    outfile=paste0(gt,"/genas_",outspec)
    for (i in 1:ncol(genas.cont)){
      genas.res=c(conts[genas.cont[,i]],unlist(genas(fit,coef=geneas.cont[,i],plot=F,subset="Fpval"))[c(1,6:9)])
      names(genas.res)[1:2]=c("cont1","cont2")
      if (as.numeric(genas.res["p.value"]) < 0.05) {
        pdf(file=paste0(outfile,"_contrasts_",paste0(conts[genas.cont[,i]],collapse="_"),".pdf"))
        bla=genas(fit,coef=geneas.cont[,i],plot=T,subset="Fpval")
        dev.off()
      } 
      if (i==1){
        genas.results=genas.res
      }
      else{
        genas.results=rbind(genas.results,genas.res)
      }
      
    }
    genas.results=data.frame(genas.results)
    write.table(genas.results,file=paste0(outfile,"_res.tab"),row.names=F,sep="\t")    
  }
  
  for (time in names(GT.times.dges[["dge.times"]][["DT"]]) ) {
    outspec=paste0("location_",loc,"_timepoint_",time)
    fit=fits[["times"]][[loc]][[time]]
    conts=colnames(fit$contrasts)
    timepoint=paste0("timepoint_",time)
    # create directory for saving files
    dir.create(timepoint,showWarnings ='F')    
    # get toptable for each contrast, merge and write to file
    tts=list()
    for (cont in conts){
      tts[[cont]]= topTable(fit,coef=cont,adjust.method="BH",number=nrow(fit$genes),sort.by="none")
      names(tts[[cont]])[7:12]=paste(names(tts[[cont]])[7:12],cont,sep="_")
    }
    ttmerge=merge(tts[[1]],tts[[2]],by=c(1:6))    
    ttmerge=merge(ttmerge,tts[[3]],by=c(1:6))
    # get decide tests over contrasts
    decT=decideTests(fit,adjust.method="BH",method="separate")
    ttmerge=cbind(ttmerge,decT)
    rownames(ttmerge)=rownames(tts[[1]])
    ttmerge=cbind(rownames(ttmerge),ttmerge)
    colnames(ttmerge)[1]="ENSEMBLEID"
    # write merged table
    write.table(ttmerge,file=paste0(timepoint,"/results_",outspec,"_all_GTs.tab"),row.names=F,sep=",")
    # plot VennDiagram for decide tests
    outfile=paste0(timepoint,"/sig_de_",outspec)
    pdf(file=paste0(outfile,"_updown.pdf"))
    vennDiagram(decT,include="both", main=paste(loc,"timepoint:",time,"diff. expressed"))
    dev.off()
    pdf(file=paste0(outfile,"_up.pdf"))
    vennDiagram(decT,include="up",main=paste(loc,"timepoint:",time,"up regulated"))
    dev.off()
    pdf(file=paste0(outfile,"_down.pdf"))
    vennDiagram(decT,include="down",,main=paste(loc,"timepoint:",time,"down regulated"))
    dev.off()
    pdf(file=paste0(outfile,"_heatdiagram.pdf"))
    heatDiagram(decT,fit$coef,names=rownames(fit$genes))
    dev.off()    
    for (cont in conts){
      # plot MA plot for each contrast
      pdf(file=paste0(timepoint,"/MA_",outspec,cont,".pdf"))
      adj.pvals=p.adjust(fit$p.value[,cont],method="BH")
      limma::plotMA(fit,coef=cont,status=as.factor(c("not sig.","sig.")[(adj.pvals <= 0.05) + 1]) ,col=c("red"),cex=0.3,main= paste(loc,timepoint,cont))
      dev.off()
      # volcano plot fo each contrast
      pdf(file=paste0(timepoint,"/vulcano_",outspec,cont,".pdf"))
      vulcano.plot(fit,cont,main= paste(loc,timepoint,cont),ylab=expression(-plain(log)[10](italic(P))), xlab=expression(plain(log)[2](italic(FC))) )
      dev.off()
    }
    # plot geneas for each contrast combination
    genas.cont=matrix(c(c(1,2),c(1,3),c(2,3)),nrow=2,ncol=3)
    genas.results=c()
    outfile=paste0(timepoint,"/genas_",outspec)
    for (i in 1:ncol(genas.cont)){
      genas.res=c(conts[genas.cont[,i]],unlist(genas(fit,coef=geneas.cont[,i],plot=F,subset="Fpval"))[c(1,6:9)])
      names(genas.res)[1:2]=c("cont1","cont2")
      if (as.numeric(genas.res["p.value"]) < 0.05) {
        pdf(file=paste0(outfile,"_contrasts_",paste0(conts[genas.cont[,i]],collapse="_"),".pdf"))
        bla=genas(fit,coef=geneas.cont[,i],plot=T,subset="Fpval")
        dev.off()
      } 
      if (i==1){
        genas.results=genas.res
      }
      else{
        genas.results=rbind(genas.results,genas.res)
      }
      
    }
    genas.results=data.frame(genas.results)
    write.table(genas.results,file=paste0(outfile,"_res.tab"),row.names=F,sep="\t")    
  }
}



quit(save="no")

save.image("thu01102015.RData")
#load("thu01102015.RData")




# gene set association
# gene signature association
contrast="t0t2"
loc="PT"
gt="WT"
propTrueNull(fits[["GTs"]][[loc]][[gt]]$p.value)
genas(fits[["GTs"]][[loc]][[gt]],coef=c(1,2),plot=T,subset="p.union")
genas(fits[["GTs"]][[loc]][[gt]],coef=c(1,3),plot=T,subset="p.union")
genas(fits[["GTs"]][[loc]][[gt]],coef=c(2,3),plot=T,subset="p.union")
genas(fits[["GTs"]][[loc]][[gt]],coef=c(1,2),plot=T,subset="all")
genas(fits[["GTs"]][[loc]][[gt]],coef=c(1,2),plot=T,subset="Fpval")
c1=topTable(fits[["GTs"]][[loc]][[gt]],coef=1,adjust.method="BH")
head(topTable(fits[["GTs"]][[loc]][[gt]],coef=1,adjust.method="BH"))
# for getting the contrast names
colnames(fits[["GTs"]][[loc]][[gt]]$p.value)
# finding up and down regulated genes in all contrasts
apply(decideTests(fit,adjust.method="BH",method="nestedF"),2,table)
resWT=decideTests(fits[["GTs"]][[loc]][[gt]],adjust.method="BH",method="separate")
# getting breakdownand plotting it
vennCounts(resWT,include="both")
vennDiagram(resWT,include="both")
# heat diagram of breakdown
heatDiagram(resWT,fits[["GTs"]][[loc]][[gt]]$coef,names=rownames(fits[["GTs"]][[loc]][[gt]]$genes))
summary(resWT)
# variance stabilisation over expression
plotSA(fits[["GTs"]][[loc]][[gt]])
# plot as lines .. not too helpful, I can do as well
plotlines(fits[["GTs"]][[loc]][[gt]])



plot.Evals("ENSMUSG00000025492",GT.times.dges[["dge.GTs"]][[loc]][[gt]])
plot.Evals("ENSMUSG00000025492",all.samps$dge,ylim=c(0,3000))
dge.sub=all.samps$dge
gene="ENSMUSG00000025492" 

contrast="t0t2"
ct=names(fits[["GTs"]][[loc]][[gt]]$contrasts[fits[["GTs"]][[loc]][[gt]]$contrasts[,contrast] != 0,1])
samps1=names(fits[["GTs"]][[loc]][[gt]]$design[fits[["GTs"]][[loc]][[gt]]$design[,ct[1]] == 1,1])
samps2=names(fits[["GTs"]][[loc]][[gt]]$design[fits[["GTs"]][[loc]][[gt]]$design[,ct[2]] == 1,1])
adj.pvals=p.adjust(fits[["GTs"]][[loc]][[gt]]$p.value[,contrast],method="BH")
library(mixOmics)
evals=cbind(GT.times.dges[["dge.GTs"]][[loc]][[gt]]$E[adj.pvals < 0.05,samps1],GT.times.dges[["dge.GTs"]][[loc]][[gt]]$E[adj.pvals < 0.05,samps2])
cim(t(evals), dendrogram = "column", xlab = "Genes", ylab = "Samples",
    col = colorRampPalette(c("red","lightgreen","blue"))(255), symkey = FALSE, lhei = c(1, 3))

melt(evals[1,])

# doing different flavors of volcano
volcanoplot(fits[["GTs"]][[loc]][[gt]])
vulcano.plot(fits[["GTs"]][[loc]][[gt]],contrast, main= paste(gt,loc,contrast),ylab=expression(-plain(log)[10](italic(P))), xlab=expression(plain(log)[2](italic(FC))) )

# MAplots
adj.pvals=p.adjust(fits[["GTs"]][[loc]][[gt]]$p.value[,contrast],method="BH")
limma::plotMA(fits[["GTs"]][[loc]][[gt]],coef=contrast,status=as.factor(c("not sig.","sig.")[(adj.pvals <= 0.05) + 1]) ,col=c("red"),cex=0.3)

# do a mighty mds plot
plotMDS(all.samps$dge, top = 500,pch=as.numeric(all.samps$dge$samples$group),gene.selection="pairwise", col=as.numeric(all.samps$dge$samples$group))
## Just for GT WT!
# define contrasts 
cont.wt <- makeContrasts(t0t2="time0-time2",t0t8="time0-time8",
                         t2t8="time2-time8",levels=GT.times.dges[["design.GTs"]][["PT"]][["WT"]])

fit <- lmFit(GT.times.dges[["dge.GTs"]][["PT"]][["WT"]], design=GT.times.dges[["design.GTs"]][["PT"]][["WT"]])
fit <- eBayes(fit)
toptable(fit)
dectest=decideTests(fit,method="global",adjust.method = "fdr")
table(apply(dectest,1,function(x) paste0(x,collapse="")))

fit2 <- contrasts.fit(fit, cont.wt)
fit2 <- eBayes(fit2)
head(fit2)
contrast="t0t8"
limma::plotMA(fit2,coef=contrast,status=as.factor(c("not sig.","sig.")[(fit2$p.value[,contrast] <= 0.05) + 1]) ,col=c("red"),cex=0.3)


# get DE genes for every contrast
contrasts=colnames(fit2$coefficients)
ttWT=list()
for (cont in contrasts){
  #ttWT[[length(ttWT)+1]]=""
  #names(ttWT[[length(ttWT)]])=cont
  ttWT[[cont]]=toptable(fit2,coef=cont,genelist = fit2$genes[,c(5)],number=nrow(fit2$genes),sort.by="none")
}

contrast="t0t8"
WTt0t8=toptable(fit2,coef="t0t8",genelist = fit2$genes[,c(5)],number=nrow(fit2$genes),sort.by="none")
limma::plotMA(fit2,coef=contrast,status=as.factor(c("not sig.","sig.")[(WTt0t8$adj.P.Val <= 0.05) + 1]) ,col=c("red"),cex=0.3)

# vulcano plot
contrast="t0t2"
vulcano.plot(fit2,contrast,ylab="log10(P)",xlab="log2(FC)",main=paste("WT contrast:",contrast))


contrast="t0t8"
WTt0t8=toptable(fit2,coef="t0t8",genelist = fit2$genes[,c(5)],number=nrow(fit2$genes),sort.by="none",method="global")
limma::plotMA(fit2,coef=contrast,status=as.factor(c("not sig.","sig.")[(WTt0t8$adj.P.Val <= 0.05) + 1]) ,col=c("red"),cex=0.3)


str(GT.times.dges[["dge.GTs"]][["PT"]][["WT"]])
str(GT.times.dges[["design.GTs"]][["PT"]][["WT"]])
tt=toptable(fit2,coef="t0t8",genelist = fit2$genes[,c(5)],number=nrow(fit2$genes),sort.by="none")
head(tt)
dectest=decideTests(fit2,method="global",adjust.method = "fdr")
summary(dectest)
table(apply(dectest,1,function(x) paste0(x,collapse="")))
dectest=decideTests(fit2,method="separate",adjust.method = "holm")

plot(fit2$Amean,-log10(fit2$F.p.value))

qq(fit2$F.p.value)

qqplot2(runif(length(fit$F.p.value)),fit$F.p.value)
qqline(fit$F.p.value,distribution = qunif)
qq(fit$F.p.value)
qq(WTt0t8$P.Value)
library(qqman)
qqplot2 <- function(op,dp,...){
  nl_oP <- -log10(op)
  nl_dP <- -log10(dp)
  s_nl_oP <- sort(nl_oP,decreasing=F)
  s_nl_dP <- sort(nl_dP,decreasing=F)
  plot(s_nl_dP,s_nl_oP,...)
}

# calculate model


# get results


#### playing around
#filename="/Volumes/Temp/Lukas/BSA_0091_Proximal_Tubule/mm10//featCounts/mm10.1000nt.featurecount.counts_edgeR.table"
filename="/Volumes/vetgrid01/BSA_0091_Proximal_Tubule/mm10//featCounts/new.test/ucsc_mm10_ensembl_exons_1000nt_most_highly.counts_edgeR.table"
dropcols=c("WT_KO4_8h_DT","Klotho_VDR_KO2_2h_PT","Klotho_VDR_KO2_2h_DT","Klotho_VDR_KO4_8h_PT","Klotho_VDR_KO4_8h_DT","Klotho_VDR_KO3_8h_DT","VDR_KO4_8h_PT")
drop_reg="^VDR_KO"
a = create.dge(filename,dropcols=dropcols,drop_reg=drop_reg,both=T)
a[["dge"]] <- calcNormFactors(a[["dge"]])
a[["dge"]]$samples
plotMDS(a[["dge"]])all.var[,5:27]

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

# load files to get enriched categories

# get list of files with enrichments:
enr.fn = list.files(path = "Highest1000/Enrichments", pattern = "*.tab")
# get categories for each file (split filename on _ but not before VDR)
enr.cats = gsub(c(".tab|enr_"),"",enr.fn)
enr.cats = strsplit(enr.cats, "_(?!VDR|BP)" ,perl= TRUE)
# in each list entry: 1: enr. category, 2: location, GT or time, 3: "t", 4: contr1, 5: contr2

enr.cats.list = list( list( list() ) )
names(enr.cats.list) = names(GT.times.dges[["dge.GTs"]])
#enr.cats.list = list()
enr.cats.list = setNames( as.list(rep.int("", length(names(GT.times.dges[["dge.GTs"]])))), names(GT.times.dges[["dge.GTs"]])  )

# set up list holding enrichment categories
gt.names = names(GT.times.dges[["dge.GTs"]][["DT"]])
time.names = names(GT.times.dges[["dge.times"]][["DT"]])
cat.names=unique(sapply(enr.cats, function(x) x[1]))
enr.cats.list = list()
enr.cats.list = setNames( as.list(rep.int("", length(names(GT.times.dges[["dge.GTs"]])))), names(GT.times.dges[["dge.GTs"]])  )
for( i in names(enr.cats.list)) {
  enr.cats.list[[i]] = setNames( as.list(rep.int("", length(c(gt.names,time.names)) ) ), c(gt.names,time.names)   )
  for (j in c(gt.names,time.names)) {
    enr.cats.list[[i]][[j]] = setNames( as.list(rep.int(list( list() ), length(cat.names) ) ), cat.names   )
  }
}

path.name="Highest1000/Enrichments"
for (i in 1:length(enr.cats)) {
  # read file into table
  set.names = read.table(paste(path.name,enr.fn[i],sep="/"),header=TRUE)
  # get FDR < 0.05 categories
  set.names =  as.character(set.names$SETNAME[set.names$FDR < 0.05])
  enr.cats.list[[ enr.cats[[i]][2] ]][[ enr.cats[[i]][3] ]][[ enr.cats[[i]][1] ]] = unique(append( 
    enr.cats.list[[ enr.cats[[i]][2] ]][[ enr.cats[[i]][3] ]][[ enr.cats[[i]][1] ]], set.names ) ) 
}

library("xlsx")
library("WriteXLS")

# use WriteXLS: create a list of data.frames, each with the category name, set name, gene names, contrast results
# WriteXLS(x, ExcelFileName = "R.xls", SheetNames = names(x), AdjWidth = TRUE, BoldHeaderRow = FALSE,FreezeRow = 1, FreezeCol = 5)
write_gene_lists <- function (loc,gt,fit,outspec) {
  conts=colnames(fit$contrasts)
  # create directory for saving files
  #dir.create(gt,showWarnings ='F') 
  # get toptable for each contrast, merge and write to file
  g.names=cbind(rownames(fit$genes),fit$genes[,c("gene_symbol","entrezgene")])
  colnames(g.names)[1]="ENSEMBLEID"
  decT=decideTests(fit,adjust.method="BH",method="separate")
  decT=cbind(rownames(decT),decT)
  colnames(decT)[1]="ENSEMBLEID"
  decTm=merge(g.names,decT,by="ENSEMBLEID")
  decTm=decTm[ apply(decTm[,c(4,5,6)], 1, function(x) any(x != "0") ),]
  enr.lists = list()
  xlsbname=paste("enr",outspec,"genes.xls",sep="_")
  for ( cat in c("GO_BP","H","KEGG") ) {
    # if empty, continue
    if (length(enr.cats.list[[loc]][[gt]][[cat]]) == 0) {
      next
    }
    # get merged vector of categories
    all.res = data.frame()
    for (set.name in unlist(enr.cats.list[[loc]][[gt]][[cat]])) {
      if (cat == "GO_BP"){
        gene.indx=which(decTm$ENSEMBLEID %in% Mm_GObp[[set.name]])
      }
      else{
        gene.indx=which(decTm$entrezgene %in% Mm.Migs[[cat]][[set.name]])
      }
      if (length(gene.indx) == 0)  { next } 
      res.df = cbind(set.name,decTm[gene.indx ,])
      all.res=rbind(all.res, res.df)
    }
    if (nrow(all.res) > 0){
      enr.lists[[ cat ]] = all.res
    }
  }
  if (length(enr.lists) == 0) return() 
  WriteXLS(enr.lists, ExcelFileName = paste0(path.name,"/",xlsbname), SheetNames = names(enr.lists), AdjWidth = TRUE, BoldHeaderRow = TRUE,FreezeRow = 1, FreezeCol = 4)
}

get_de_lists <- function (loc,gt,fit,outspec) {
  conts=colnames(fit$contrasts)
  de.list=list()
  # create directory for saving files
  #dir.create(gt,showWarnings ='F') 
  # get toptable for each contrast, merge and write to file
  tts=list()
  conts=colnames(fit$contrasts)
  for (cont in conts){
    tts[[cont]]= topTable(fit,coef=cont,adjust.method="BH",number=nrow(fit$genes),sort.by="none")
    tts[[cont]]= tts[[cont]][,c("Chr","gene_symbol","entrezgene","AveExpr","logFC","P.Value","adj.P.Val") ]
    #tts[[cont]] = tts[[cont]][,! (names(tts[[cont]]) %in% c("Start", "End", "Length", "t", "B")) ]
    names(tts[[cont]])[c(5:7)]=paste(names(tts[[cont]])[c(5:7)],cont,sep="_")
    tts[[cont]]=cbind(rownames(tts[[cont]]),tts[[cont]])
    colnames(tts[[cont]])[1]="ENSEMBLEID"      
  }
  ttmerge=merge(tts[[1]],tts[[2]][,-c(2:5)],by="ENSEMBLEID")    
  ttmerge=merge(ttmerge,tts[[3]][,-c(2:5)],by="ENSEMBLEID")
  rownames(ttmerge)=rownames(tts[[1]])
  decT=decideTests(fit,adjust.method="BH",method="separate")
  decTt=cbind(rownames(decT),decT)
  colnames(decTt)[1]="ENSEMBLEID"
  decT.df=data.frame()
  decT.df=rbind(decT.df,summary(decT)[,])
  # create merged 
  decT.df=cbind(rownames(decT.df),decT.df)
  colnames(decT.df)[1]="Regulation"
  # create merged DE list
  ttmerge=merge(ttmerge,decTt,by="ENSEMBLEID")
  de.index = which( ( abs(ttmerge[ , 6]) >= 1.0 &  ttmerge[ , 8] <= 0.05 ) | ( abs(ttmerge[ , 9]) >= 1.0 &  ttmerge[ , 11] <= 0.05 ) 
                    | ( abs(ttmerge[ , 12]) >= 1.0 &  ttmerge[ , 14] <= 0.05 ) )
  de.list[[ paste0("summary_",outspec) ]] = decT.df
  de.list[[ paste0("Sign_2x_",outspec) ]] = ttmerge[de.index,]
  de.list[[ paste0("all_genes_",outspec) ]] = ttmerge
  return(de.list)
}

path.name="Highest1000/Enrichments"
for (loc in names(GT.times.dges[["dge.times"]]) ){
  for (gt in names(GT.times.dges[["dge.GTs"]][["DT"]])){
    print(paste("going through",loc,gt))
    outspec=paste0("location_",loc,"_genotype_",gt)
    fit=fits[["GTs"]][[loc]][[gt]]
    write_gene_lists(loc,gt,fit,outspec)
  }
  for (time in names(GT.times.dges[["dge.times"]][[ loc]]) ) {
    outspec=paste0("location_",loc,"_timepoint_",time)
    fit=fits[["times"]][[loc]][[time]]
    conts=colnames(fit$contrasts)
    write_gene_lists(loc,time,fit,outspec)
  }
}

path.name="Highest1000"
for (gt in names(GT.times.dges[["dge.GTs"]][["DT"]])){
  loc.list=list()
  for (loc in names(GT.times.dges[["dge.times"]]) ){
    print(paste("going through",loc,gt))
    outspec=paste0("loc_",loc,"_gt_",gt)
    fit=fits[["GTs"]][[loc]][[gt]]
    loc.list=append(loc.list,get_de_lists(loc,gt,fit,outspec))
  }
  xlsbname=paste("de_genes_genotype",gt,".xls",sep="_")
  WriteXLS(loc.list, ExcelFileName = paste0(path.name,"/",gt,"/",xlsbname), SheetNames = names(loc.list), AdjWidth = TRUE, BoldHeaderRow = TRUE,FreezeRow = 1, FreezeCol = 1)
}
for (time in names(GT.times.dges[["dge.times"]][["DT"]]) ) {
  loc.list=list()
  timepoint=paste0("timepoint_",time)
  for (loc in names(GT.times.dges[["dge.times"]]) ){
    outspec=paste0("loc_",loc,"_time_",time)
    fit=fits[["times"]][[loc]][[time]]
    conts=colnames(fit$contrasts)
    loc.list=append(loc.list,get_de_lists(loc,time,fit,outspec))
  }
  xlsbname=paste("de_genes_timepoint_",time,".xls",sep="_")
  WriteXLS(loc.list, ExcelFileName = paste0(path.name,"/",timepoint,"/",xlsbname), SheetNames = names(loc.list), AdjWidth = TRUE, BoldHeaderRow = TRUE,FreezeRow = 1, FreezeCol = 1)
  
}
