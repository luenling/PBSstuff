#setwd("/Volumes/Temp/Lukas/Data/Test/")
setwd("/Volumes/Temp/Lukas/Data/Test/")
library(limma)
library(edgeR)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library("biomaRt")
library(gplots)
library(EnrichmentBrowser)
library(KEGGREST)
require(grDevices)
require(Biobase)
#library(mixOmics)
library("WriteXLS")

do_gs_ebrowser<- function(fit,dge,gt,loc,perm=10000){
  # create eset for each contrast
  db.list = list()
  gt.name=gt
  if ( nchar(gt.name) < 2 ) {
    gt.name=paste("time",gt,sep="_")
  } 
  for (cont in colnames(fit$contrasts)) {
    featureData <- topTable(fit,coef=cont,adjust.method="BH",number=nrow(fit$genes),sort.by="none")[,c("logFC","adj.P.Val","t")]
    colnames(featureData) = c("FC","ADJ.PVAL","limma.STAT")
    featureData <- new("AnnotatedDataFrame", data=featureData)
    exprs=as.matrix(dge)
    rownames(exprs) = rownames(dge)
    exprs=exprs[rownames(featureData),]
    pdata=as.data.frame(dge$targets)
    colnames(pdata)=c("GROUP","lib.size","norm.factors")
    pdata = new("AnnotatedDataFrame", data=pdata)
    eset<-new("ExpressionSet", exprs=exprs, featureData = featureData, phenoData = pdata)
    annotation(eset) <- "org.Mm.eg"
    eset <- map.ids(eset, org="mmu", from="ENSEMBL", to="ENTREZID")
    # do enrichments
    bname=paste("Enrichments/enr",loc,gt.name,sep="_")
    for( db in names(Mm.Migs)){
      if ( is.null(db.list[[db]]) ) {db.list[[db]] = list()}
      sbea.gsenr <- sbea( method = use.gsenr , eset=eset, gs=Mm.Migs[[db]] , alpha = 0.05, perm = perm, padj.method = "BH")
      sbea.gsenr$method <- paste(db,"enr",loc,gt,cont,sep="_")
      sbea.gsenr$res.tbl$NR.SIG.GENES <- get.nsig(sbea.gsenr)
      sbea.gsenr$res.tbl$NR.GENES <- get.ngenes(sbea.gsenr)
      sbea.gsenr$res.tbl <- sbea.gsenr$res.tbl[,c( "GENE.SET","NR.GENES","NR.SIG.GENES","P.VALUE" )]
      if (! is.null(gs.ranking(sbea.gsenr)) ){
        dir.create(paste(getwd(),bname,sep="/"),showWarnings ='F',recursive=T) 
        config.ebrowser( "OUTDIR.DEFAULT", value = paste(getwd(),bname,sep="/") )
        annotation(sbea.gsenr$eset) = "mmu"
        #        ea.browse(sbea.gsenr, html.only = TRUE)
        db.list[[db]][[paste(db,loc,gt.name,cont,sep="_")]] = as.data.frame(gs.ranking(sbea.gsenr))
      }
    }
  }
  all.db = list()
  for(db in names(db.list)){
    if (length(db.list[[db]]) < 1) {
      next
    }
    ebname=paste("Enrichments/enr",db,loc,gt.name,"FDR0.05.xls",sep="_")
    WriteXLS(db.list[[db]], ExcelFileName = ebname, SheetNames = names(db.list[[db]]), AdjWidth = TRUE, BoldHeaderRow = TRUE,FreezeRow = 1, FreezeCol = 0)
    all.db[[db]] = c()
    for (cont in names(db.list[[db]])) {
      all.db[[db]] = unique( append(all.db[[db]],db.list[[db]][[cont]]$GENE.SET) )
    }
  }
  ebname=paste("Enrichments/enr",loc,gt.name,"genes.xls",sep="_")
  g.names=cbind(rownames(fit$genes),fit$genes[,c("gene_symbol","entrezgene")])
  colnames(g.names)[1]="ENSEMBLEID"
  decT=decideTests(fit,adjust.method="BH",method="separate")
  decT=cbind(rownames(decT),decT)
  colnames(decT)[1]="ENSEMBLEID"
  decTm=merge(g.names,decT,by="ENSEMBLEID")
  decTm=decTm[ apply(decTm[,c(4,5,6)], 1, function(x) any(x != "0") ),]
  decTm$Gene.Sets = ""
  de.sets = list()
  for(db in names(all.db)){
    decTm$Gene.Sets = ""
    for(set.name in all.db[[db]]) {
      gene.indx=which(decTm$entrezgene %in% Mm.Migs[[db]][[set.name]])
      for (idx in gene.indx){
        decTm$Gene.Sets[idx] = paste(decTm$Gene.Sets[idx],set.name,sep=" ")
      }
    }
    de.set = decTm[nchar(decTm$Gene.Sets) > 1, ]
    if ( nrow(de.set) > 0) {
      de.sets[[db]] = de.set
    }
  }
  WriteXLS(de.sets, ExcelFileName = ebname, SheetNames = names(de.sets), AdjWidth = TRUE, BoldHeaderRow = TRUE,FreezeRow = 1, FreezeCol = 4)
}  


get_ttmerge <- function(fit){
  tts=list()
  conts=colnames(fit$contrasts)
  for (cont in conts){
    tts[[cont]]= topTable(fit,coef=cont,adjust.method="BH",number=nrow(fit$genes),sort.by="none")
    tts[[cont]]= tts[[cont]][,c("Chr","gene_symbol","entrezgene","AveExpr","logFC","adj.P.Val") ]
    names(tts[[cont]])[c(5:6)]=paste(names(tts[[cont]])[c(5:6)],cont,sep="_")
    tts[[cont]]=cbind(rownames(tts[[cont]]),tts[[cont]])
    colnames(tts[[cont]])[1]="ENSEMBLEID"      
  }
  ttmerge=merge(tts[[1]],tts[[2]][,-c(2:5)],by="ENSEMBLEID")    
  ttmerge=merge(ttmerge,tts[[3]][,-c(2:5)],by="ENSEMBLEID")
  rownames(ttmerge)=ttmerge$ENSEMBLEID
  return(ttmerge)
}


get_de_lists <- function (fit,outspec,dge.sub=NULL) {
  conts=colnames(fit$contrasts)
  de.list=list()
  # create directory for saving files
  #dir.create(gt,showWarnings ='F') 
  # get toptable for each contrast, merge and write to file
  conts=colnames(fit$contrasts)
  ttmerge = get_ttmerge(fit)
  # round fc, signif pvalues and  round ave expr
  ttmerge[,grep("logFC",colnames(ttmerge))]=round(ttmerge[,grep("logFC",colnames(ttmerge))],3)
  ttmerge[,grep("AveExpr",colnames(ttmerge))]=round(ttmerge[,grep("AveExpr",colnames(ttmerge))],2)
  ttmerge[,grep("P\\.Val",colnames(ttmerge))]= signif(ttmerge[,grep("P\\.Val",colnames(ttmerge))],3)
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
  # add RPKMs 
  if ( ! is.null(dge.sub)) {
    ttmerge = cbind(ttmerge, round(get.median.rpkm(dge.sub, genes = as.character(ttmerge$ENSEMBLEID)),2) )
  }
  de.index = which( ( abs(ttmerge[ , 6]) >= 1.0 &  ttmerge[ , 7] <= 0.05 ) | ( abs(ttmerge[ , 8]) >= 1.0 &  ttmerge[ , 9] <= 0.05 ) 
                    | ( abs(ttmerge[ , 10]) >= 1.0 &  ttmerge[ , 11] <= 0.05 ) )
  de.list[[ paste0("summary_",outspec) ]] = decT.df
  de.list[[ paste0("Sign_2x_",outspec) ]] = ttmerge[de.index,]
  de.list[[ paste0("all_genes_",outspec) ]] = ttmerge
  return(de.list)
}

use.gsenr <- function(eset , gs, alpha = 0.05, perm = 9999) {
  genes = rownames(fData(eset))
  t.values = fData(eset)$limma.STAT
  names(t.values) = rownames(fData(eset))
  gs = sapply( gs, function(x) x[ x %in% genes] )
  gStests=sapply(gs, function (x) geneSetTest(x,t.values,alternative="mixed",type="auto",ranks.only=FALSE,nsim=perm) )
  #gStestsBH=p.adjust(gStests,method = "BH")
  return(gStests)
}

get.nsig <- function(res) {
  p.values = fData(res$eset)$ADJ.PVAL
  sig.gn = rownames(fData(res$eset))[p.values < 0.05]
  gns = res$res.tbl$GENE.SET
  NR.SIG.GENES = sapply(gns, function (x) sum( sig.gn  %in% res$gs[[x]] ))
  return(NR.SIG.GENES)
}

get.ngenes <- function(res) {
  gns = res$res.tbl$GENE.SET
  NR.GENES = sapply(gns, function (x) sum( rownames(fData(res$eset))  %in% res$gs[[x]] ))
  return(NR.GENES)
}


get.median.rpkm <- function (dge.sub, genes = c(), all.samples = all.samps$dge){
  # gets a dge for a certain model and only returns the median rpkms for the groups of samples used in the model from all.samps$dge
  # create list of groups containing sample names
  #dge.sub=GT.times.dges$dge.GTs$DT$WT
  if (! "rpkm" %in% names(all.samps$dge) ) {
    all.samples$rpkm=rpkm(all.samps$dge)
  }
  if (length(genes) == 0) {
    genes=rownames(dge.sub$genes)
  }
  smp_names = split(rownames(dge.sub$targets),dge.sub$targets$group)
  return( 
    as.data.frame( sapply(smp_names, 
                          function(x) apply(all.samples$rpkm[genes,x],1,median, na.rm = TRUE ) ) ) )
}

print.sig.heatmaps <- function(fit,loc,gt,colbr=75,dir=".",p.max = 0.05, min.lfc = 0.0, min.Amean = -10){
  cols1=colorRampPalette(c("green",  "red"))(n = colbr)
  cols2=colorRampPalette(c("green",  "red"))(n = colbr)
  dgetype="dge.GTs"
  if (gt %in% c("0","2","8")) {
    dgetype="dge.times"
  }
  bn=paste0(dir,"/sig_genes_heatmap")
  for(t.name in dimnames(fit$p.value)$Contrasts ){
    sig=dimnames(fit$p.value)[[1]][p.adjust(fit$p.value[,t.name]) < p.max & abs(fit$coefficients[,t.name]) >= min.lfc & fit$Amean >= min.Amean ]
    if (length(sig) == 0) {next()}
    exprs=GT.times.dges[[dgetype]][[loc]][[gt]][["E"]][sig,]
    exprs.z=t(scale(t(exprs)))
    pheatmap(exprs.z, col=cols1, key=TRUE, symkey=FALSE, density.info="none", trace="none", scale="row", keysize=1.0, margin=c(6,6),
             clustering_distance_rows="correlation",clustering_distance_cols="correlation", fontsize_row = 10/(1+2*floor(nrow(exprs)/100)),
             filename=paste(bn,loc,gt,t.name,"scaled.pdf",sep="_"), width = 8 , height = 10)
    pheatmap(exprs, col=cols2,key=TRUE, symkey=FALSE, density.info="none", trace="none", scale="none", keysize=1.0, margin=c(6,6),
             clustering_distance_rows="correlation",clustering_distance_cols="correlation", fontsize_row = 10/(1+2*floor(nrow(exprs)/100)),
             filename=paste(bn,loc,gt,t.name,"raw.pdf",sep="_"), width = 8 , height = 10)
#     pdf(file = paste(bn,loc,gt,t.name,"scaled.pdf",sep="_"))
#     heatmap.2(exprs.z, col=cols1, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5/(1+floor(nrow(exprs)/100)),cexCol=0.75, 
#               scale="row", keysize=1.0, margin=c(6,6))
#     dev.off()
#     pdf(file = paste(bn,loc,gt,t.name,"raw.pdf",sep="_"))
#     heatmap.2(exprs, col=cols2,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5/(1+floor(nrow(exprs)/100)), 
#               cexCol=0.75, scale="none", symbreaks=F, keysize=1.0, margin=c(6,6))
#     dev.off()
  }
}



do_gs_enrichments <- function(ttmerge,gt,loc) {
  t.names=names(ttmerge[grep("^t_",names(ttmerge))])  
  for( db in names(Mm.Migs)){
    ebname=paste("Enrichments/enr",db,loc,gt,sep="_")
    indices=ids2indices(Mm.Migs[[db]],ttmerge$entrezgene)
    for( t.name in t.names){
        gser_test(ttmerge,indices,t.name,ebname,loc,gt)
      }
  }
  indices=ids2indices(Mm_GObp,ttmerge$ENSEMBLEID)
  db="GO_BP"
  ebname=paste("Enrichments/enr",db,loc,gt,sep="_") 
  for( t.name in t.names){
    gser_test(ttmerge,indices,t.name,ebname,loc,gt)
  }
}

gser_test <- function(ttmerge,indices,t.name,ebname,loc,gt) {
  gStests=c()
  for( tset in names(indices)){
    pv=geneSetTest(indices[[tset]],ttmerge[,t.name],alternative="mixed",type="auto",ranks.only=FALSE,nsim=9999)
    gStests=append(gStests,pv)
    names(gStests)[length(gStests)]=tset
  }
  gStestsBH=p.adjust(gStests,method = "BH")
  gSdf=data.frame(SETNAME=names(gStests),Pvalue=gStests, FDR=gStestsBH)
  gSdf=gSdf[order(gSdf$Pvalue),]
  if (min(gSdf$FDR) > 0.05) { return()}
  else{
    # write table
    write.table(gSdf,file=paste0(ebname,"_",t.name,".tab"),row.names=F,sep="\t")
    bc.plots = gSdf$SETNAME[gSdf$FDR > 0.05]
    if ( length(bc.plots) > 3) {
      bc.plots=bc.plots[1:3]
    }
    for (sn in bc.plots){
      fn=gsub('[^a-zA-Z0-9_]', '_', sn)
      pdf(file=paste0(ebname,"_",t.name,"_",fn,"_bcplot.pdf"))
      barcodeplot(ttmerge[,t.name],indices[[sn]],main=paste(sn,loc,gt,t.name)  )
      
      dev.off()
    }
  }
}
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
  # use first, so as to not loose any potential annotations
  genes$gene_symbol <- convertIDs(row.names(genes), "ENSEMBL", "SYMBOL", org.Mm.eg.db,"useFirst")
  genes$entrezgene <- convertIDs(row.names(genes), "ENSEMBL", "ENTREZID", org.Mm.eg.db,"useFirst")
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

get.GTandTime.dges <- function(dge,sample.mat,mincount=10, combine=FALSE){
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
  if (class(dge.sub) == "DGEList"){
    counts = TRUE
    samp.groups=as.data.frame(rownames(dge.sub$samples))
    samp.groups=cbind(samp.groups,dge.sub$samples$group)
    colnames(samp.groups)=c("value","L1")
    evals=melt(rpkm(dge.sub)[gene,])
  }
  else{
    counts = FALSE
    samp.groups=melt(lapply(as.data.frame(dge.sub$design),function(x) rownames(dge.sub$design[x == 1,])))
    evals=melt(dge.sub$E[gene,])    
  }
  evals=cbind(evals,as.factor(samp.groups$L1[sapply(rownames(evals),function (x) which(x == samp.groups$value))]))
  colnames(evals) = c("E","samples")
  if (counts) {
    ylab="RPKM"
  }
  else {
    ylab="E"
  }
  if ("genes" %in% names(dge.sub) && ! is.na(dge.sub$genes[gene,"gene_symbol"]) ){
    gene = paste0( gene, " (",dge.sub$genes[gene,"gene_symbol"],")")
  }
  qplot(samples, E, ylab= ylab, data=evals, main=gene, ...) + theme_bw()  +  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  
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

args<-commandArgs(trailingOnly=TRUE)

filename=args[1]
#dir=arg[2]
dir.create(args[2],showWarnings ='F')
setwd(paste0("/Volumes/Temp/Lukas/Data/Test/",args[2]))

#filename="/Volumes/Temp/Lukas/BSA_0091_Proximal_Tubule/mm10//featCounts/new.test/ucsc_mm10_ensembl_exons_1000nt_most_highly.counts_edgeR.table.gz"
#filename="/Volumes/vetgrid01/BSA_0091_Proximal_Tubule/mm10//featCounts/new.test/ucsc_mm10_ensembl_exons_1000nt_most_highly.counts_edgeR.table.gz"
#
dropcols=c("WT_KO4_8h_DT","Klotho_VDR_KO2_2h_PT","Klotho_VDR_KO2_2h_DT","Klotho_VDR_KO4_8h_PT","Klotho_VDR_KO4_8h_DT","Klotho_VDR_KO3_8h_DT","VDR_KO4_8h_PT")
drop_reg="^VDR_KO"
all.samps = create.dge(filename,dropcols=dropcols,drop_reg=drop_reg,both=T,mincount=15)
#all.samps.test=create.dge(filename,dropcols=dropcols,drop_reg=drop_reg,both=T,mincount=15)
design.mt =model.matrix( ~ 0 + GT + time + loc, data=all.samps$sample.mat ) 
#dupcor1=duplicateCorrelation(all.samps$dge[[1]],design=model.matrix( ~ 0 + all.samps$sample.mat$group))
#dupcor2=duplicateCorrelation(all.samps$dge[[1]],design=design.mt)
# calculate the TMM normalisation factors for library sizes
all.samps$dge=calcNormFactors(all.samps$dge)
#a1.test = voom(all.samps$dge, plot=TRUE)
#a2.test = voom(all.samps$dge, design= design.mt, plot=TRUE)
enrich=FALSE
if (length(args > 2) && args[3] == "enrichment") {
  dir.create("Enrichments",showWarnings ='F')
  print("Loading DB for enrichment")
  # setting up databases for enrichment analysis
  load("/Volumes/Temp/Lukas/Data/Test/MSigDB//mouse_H_v5.rdata")
  load("/Volumes/Temp/Lukas/Data/Test/MSigDB//mouse_c2_v5.rdata")
  load("/Volumes/Temp/Lukas/Data/Test/MSigDB//mouse_c3_v5.rdata")
  #load("MSigDB//mouse_c4_v5.rdata")
  #load("MSigDB//mouse_c6_v5.rdata")
  #load("MSigDB//mouse_c7_v5.rdata")
  Mm.Migs=list()
  Mm.Migs[["H"]]=Mm.H
  #KEGG
  Mm.Migs[["KEGG"]]=Mm.c2[grep("KEGG",names(Mm.c2))]
  Mm.Migs[["KEGG"]]=Mm.Migs[["KEGG"]][ sapply(Mm.Migs[["KEGG"]], function(x) (length(x) >= 4 )) ]
  #REACTOME
  Mm.Migs[["REACT"]]=Mm.c2[grep("REACTOME",names(Mm.c2),ignore.case = TRUE)]
  Mm.Migs[["REACT"]]=Mm.Migs[["REACT"]][ sapply(Mm.Migs[["REACT"]], function(x) (length(x) >= 4 )) ]
  #MIR
  Mm.Migs[["MIR"]]=Mm.c3[grep("mir",names(Mm.c3),ignore.case = TRUE)]
  Mm.Migs[["MIR"]]=Mm.Migs[["MIR"]][ sapply(Mm.Migs[["MIR"]], function(x) (length(x) >= 3 )) ]
  #TF
  Mm.Migs[["TF"]]=Mm.c3[ - grep("mir",names(Mm.c3),ignore.case = TRUE)]
  Mm.Migs[["TF"]]=Mm.Migs[["TF"]][ - grep("UNKNOWN",names(Mm.Migs[["TF"]])) ]
  Mm.Migs[["TF"]]=Mm.Migs[["TF"]][ sapply(Mm.Migs[["TF"]], function(x) (length(x) >= 3  )) ]
  print("setting up GO files")
  # main mart is dead use this one
  ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host = "jul2015.archive.ensembl.org")
  datasets <- listDatasets(ensembl)
  ensembl <- useDataset("mmusculus_gene_ensembl",ensembl)
  Mm_GO <- getBM(attributes = c("go_id","ensembl_gene_id","namespace_1003","name_1006"),filters="ensembl_gene_id", values=rownames(all.samps$dge$genes), mart = ensembl)
  Mm_GO$namespace_1003=as.factor(Mm_GO$namespace_1003)
  Mm_GO=Mm_GO[Mm_GO$namespace_1003 == "biological_process",]
  Mm_GO$name_1006=toupper(gsub(' ',"_",Mm_GO$name_1006))
  Mm_GO$name_1006=gsub('-',"_",Mm_GO$name_1006)
  Mm_GObp = split(Mm_GO$ensembl_gene_id,as.factor(Mm_GO$name_1006))
  Mm_GObp = Mm_GObp[ sapply(Mm_GObp, function(x) (length(x) >= 4 )) ]
  enrich=TRUE  
}

#library(reshape2)
df = melt(as.data.frame(log2(all.samps$dge[[1]]+1)), variable.name = "Samples")
#boxplot(df$value ~ df$Samples)
# all look very similar - median at around 2^5=32
## create list of subobjects with only relevant samples for GT and time
print(paste("creating DGE tables"))
GT.times.dges = get.GTandTime.dges(all.samps$dge,all.samps$sample.mat,mincount=25)



# define the model matrices for the different GTs and times (already correctly done in the function get.GTandTime.dges)
# for(loc in names(GT.times.dges[["dge.GTs"]]) ){
#   for(gt in names(GT.times.dges[["dge.GTs"]][["DT"]]) ) {
#     sample.mat.gt <- GT.times.dges[["sm.GTs"]][[loc]][[gt]]
#     GT.times.dges[["design.GTs"]][[loc]][[gt]] <- model.matrix( ~  0 + time, data = sample.mat.gt)
#   }
#   for(time in names(GT.times.dges[["dge.times"]][["DT"]]) ) {
#     sample.mat.times <- GT.times.dges[["sm.times"]][[loc]][[time]]
#     GT.times.dges[["design.times"]][[loc]][[time]] <- model.matrix( ~ 0 + GT, data = sample.mat.times)
#   }
# }

#save.image("wed30092015.RData")
#load("wed30092015.RData")
# correct gene names
# print("correcting genotype gene names")
# locs = names(GT.times.dges[["dge.GTs"]])
# for (loc in locs ) {
#   for (gt in names(GT.times.dges[["dge.GTs"]][["DT"]]) ) {
#     genes=fits[["GTs"]][[loc]][[gt]][["genes"]]
#     genes$gene_symbol <- convertIDs(row.names(genes), "ENSEMBL", "SYMBOL", org.Mm.eg.db,"useFirst")
#     genes$entrezgene <- convertIDs(row.names(genes), "ENSEMBL", "ENTREZID", org.Mm.eg.db,"useFirst")
#     if (sum(rownames(genes) != rownames(GT.times.dges[["dge.GTs"]][[loc]][[gt]][["genes"]])) > 0) {
#       print("problem, shit")
#     }
#     fits[["GTs"]][[loc]][[gt]][["genes"]] = genes
#     GT.times.dges[["dge.GTs"]][[loc]][[gt]][["genes"]] = genes
#   }
# }
# print("creating fitted objects for timepoints")
# for (loc in locs) {
#   for (time in names(GT.times.dges[["dge.times"]][["DT"]]) ) {
#     genes=fits[["times"]][[loc]][[time]][["genes"]]
#     genes$gene_symbol <- convertIDs(row.names(genes), "ENSEMBL", "SYMBOL", org.Mm.eg.db,"useFirst")
#     genes$entrezgene <- convertIDs(row.names(genes), "ENSEMBL", "ENTREZID", org.Mm.eg.db,"useFirst")
#     if (sum(rownames(genes) != rownames(GT.times.dges[["dge.times"]][[loc]][[time]][["genes"]])) > 0) {
#       print("problem, shit")
#     }
#     fits[["times"]][[loc]][[time]][["genes"]] = genes
#     GT.times.dges[["dge.times"]][[loc]][[time]][["genes"]] = genes
#   }
# }



# loop over all GTs
print("creating fitted objects for genotypes")
locs = names(GT.times.dges[["dge.GTs"]])
if (length(args > 3) &&  args[4] %in% names(GT.times.dges[["dge.GTs"]]) )  {
  locs = c( args[4] )  
}
fits=list(list(list()))
for (loc in locs ) {
  for (gt in names(GT.times.dges[["dge.GTs"]][["DT"]]) ) {
    conts <- makeContrasts(t0t2="time0-time2",t0t8="time0-time8",
                           t2t8="time2-time8",levels=GT.times.dges[["design.GTs"]][[loc]][[gt]])
    fit <- lmFit(GT.times.dges[["dge.GTs"]][[loc]][[gt]], design=GT.times.dges[["design.GTs"]][[loc]][[gt]])
    fit <- contrasts.fit(fit, conts)
    fit <- eBayes(fit)
    fits[["GTs"]][[loc]][[gt]] <- fit
  }
}
print("creating fitted objects for timepoints")
for (loc in locs) {
  for (time in names(GT.times.dges[["dge.times"]][["DT"]]) ) {
    conts <- makeContrasts(WT_FGF23="GTWT-GTFGF23_VDR",WT_Klotho="GTWT-GTKlotho_VDR",
                           FGF23_Klotho="GTFGF23_VDR-GTKlotho_VDR",levels=GT.times.dges[["design.times"]][[loc]][[time]])
    fit <- lmFit(GT.times.dges[["dge.times"]][[loc]][[time]], design=GT.times.dges[["design.times"]][[loc]][[time]])
    fit <- contrasts.fit(fit, conts)
    fit <- eBayes(fit)
    fits[["times"]][[loc]][[time]] <- fit
  }
}

#save.image("wed20012016.RData")
#load("wed20012016.RData")
print("going through genotypes")

for (loc in locs) {
  for (gt in names(GT.times.dges[["dge.GTs"]][[loc]])){
    print(paste("going through",loc,gt))
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
      tts[[cont]]=cbind(rownames(tts[[cont]]),tts[[cont]])
      colnames(tts[[cont]])[1]="ENSEMBLEID"      
    }
    ttmerge=merge(tts[[1]],tts[[2]][,-c(2:7)],by="ENSEMBLEID")    
    ttmerge=merge(ttmerge,tts[[3]][,-c(2:7)],by="ENSEMBLEID")
    rownames(ttmerge)=ttmerge$ENSEMBLEID
    if (enrich) {
      do_gs_enrichments(ttmerge,gt,loc)
    }
    else {
      # get decide tests over contrasts
      decT=decideTests(fit,adjust.method="BH",method="separate")
      decTt=cbind(rownames(decT),decT)
      write.table(file = paste0(gt,"/decT_",outspec,"_output.tab"),summary(decT))
      next()
      colnames(decTt)[1]="ENSEMBLEID"
      ttmerge=merge(ttmerge,decTt,by="ENSEMBLEID")
      rownames(ttmerge)=ttmerge$ENSEMBLEID
      #     ttmerge=cbind(ttmerge,decT)
      #     ttmerge=cbind(rownames(ttmerge),ttmerge)
      #     colnames(ttmerge)[1]="ENSEMBLEID"
      # write merged table
      print(paste("heat map for sig genes for",loc,gt))
      print.sig.heatmaps(fit,loc,gt,dir=gt)
      print(paste("writing all times table for",loc,gt))
      write.table(ttmerge,file=paste0(gt,"/results_",outspec,"_all_times.tab"),row.names=F,sep=",")
      # plot VennDiagram for decide tests
      outfile=paste0(gt,"/sig_de_",outspec)
      pdf(file=paste0(outfile,"_updown.pdf"))
      print(paste("venn diagram for diff",loc,gt))
      vennDiagram(decT,include="both", main=paste(loc,gt,"diff. expressed"))
      dev.off()
      pdf(file=paste0(outfile,"_up.pdf"))
      print(paste("venn diagram for diff up",loc,gt))
      vennDiagram(decT,include="up", main=paste(loc,gt,"up regulated"))
      dev.off()
      pdf(file=paste0(outfile,"_down.pdf"))
      print(paste("venn diagram for diff up",loc,gt))
      vennDiagram(decT,include="down", main=paste(loc,gt,"down regulated"))
      dev.off()
      pdf(file=paste0(outfile,"_heatdiagram.pdf"))
      print(paste("heat diagram for",loc,gt))
      heatDiagram(decT,fit$coef,names=rownames(fit$genes))
      dev.off()
      print(paste("ma plots and volcanos for contrasts for",loc,gt))
      for (cont in conts){
        # plot MA plot for each contrast
        pdf(file=paste0(gt,"/MA_",outspec,cont,".pdf"))
        adj.pvals=p.adjust(fit$p.value[,cont],method="BH")
        limma::plotMA(fit,coef=cont,status=as.factor(c("not sig.","sig.")[(adj.pvals <= 0.05) + 1]) 
                      ,col=c("red"),cex=0.3,main= paste(loc,gt,cont))
        dev.off()
        # volcano plot fo each contrast
        pdf(file=paste0(gt,"/vulcano_",outspec,cont,".pdf"))
        vulcano.plot(fit,cont,main= paste(loc,gt,cont),ylab=expression(-plain(log)[10](italic(P))), xlab=expression(plain(log)[2](italic(FC))) )
        dev.off()
      }
      # plot genas for each contrast combination
      print(paste("create geneas plots for each contrast", loc,gt))
      genas.cont=matrix(c(c(1,2),c(1,3),c(2,3)),nrow=2,ncol=3)
      genas.results=c()
      outfile=paste0(gt,"/genas_",outspec)
      for (i in 1:ncol(genas.cont)){
        genas.res=c(conts[genas.cont[,i]],unlist(genas(fit,coef=genas.cont[,i],plot=F,subset="Fpval"))[c(1,6:9)])
        names(genas.res)[1:2]=c("cont1","cont2")
        if (as.numeric(genas.res["p.value"]) < 0.05) {
          pdf(file=paste0(outfile,"_contrasts_",paste0(conts[genas.cont[,i]],collapse="_"),".pdf"))
          bla=genas(fit,coef=genas.cont[,i],plot=T,subset="Fpval")
          dev.off()
        } 
        if (i==1){
          genas.results=genas.res
        }
        else{
          genas.results=rbind(genas.results,genas.res)
        }
        
      }
      row.names(genas.results)=NULL
      genas.results=data.frame(genas.results)
      write.table(genas.results,file=paste0(outfile,"_res.tab"),row.names=F,sep="\t")
    }
  }
  
  for (time in names(GT.times.dges[["dge.times"]][[ loc]]) ) {
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
      tts[[cont]]=cbind(rownames(tts[[cont]]),tts[[cont]])
      colnames(tts[[cont]])[1]="ENSEMBLEID"
      
    }
    ttmerge=merge(tts[[1]],tts[[2]][,-c(2:7)],by="ENSEMBLEID")    
    ttmerge=merge(ttmerge,tts[[3]][,-c(2:7)],by="ENSEMBLEID")
    rownames(ttmerge)=rownames(tts[[1]])
    if (enrich) {
      do_gs_enrichments(ttmerge,time,loc)
    }
    else {
      # get decide tests over contrasts
      decT=decideTests(fit,adjust.method="BH",method="separate")
      decTt=cbind(rownames(decT),decT)
      write.table(file = paste0(timepoint,"/decT_",outspec,"_output.tab"),summary(decT))
      next()
      # create heat maps
      print(paste("heat map for sig genes for",loc,time))
      print.sig.heatmaps(fit,loc,time,dir=timepoint)
      
      colnames(decTt)[1]="ENSEMBLEID"
      ttmerge=merge(ttmerge,decTt,by="ENSEMBLEID")
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
      vennDiagram(decT,include="down",main=paste(loc,"timepoint:",time,"down regulated"))
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
      # plot genas for each contrast combination
      genas.cont=matrix(c(c(1,2),c(1,3),c(2,3)),nrow=2,ncol=3)
      genas.results=c()
      outfile=paste0(timepoint,"/genas_",outspec)
      for (i in 1:ncol(genas.cont)){
        genas.res=c(conts[genas.cont[,i]],unlist(genas(fit,coef=genas.cont[,i],plot=F,subset="Fpval"))[c(1,6:9)])
        names(genas.res)[1:2]=c("cont1","cont2")
        if (as.numeric(genas.res["p.value"]) < 0.05) {
          pdf(file=paste0(outfile,"_contrasts_",paste0(conts[genas.cont[,i]],collapse="_"),".pdf"))
          bla=genas(fit,coef=genas.cont[,i],plot=T,subset="Fpval")
          dev.off()
        } 
        if (i==1){
          genas.results=genas.res
        }
        else{
          genas.results=rbind(genas.results,genas.res)
        }
        
      }
      row.names(genas.results)=NULL
      genas.results=data.frame(genas.results)
      write.table(genas.results,file=paste0(outfile,"_res.tab"),row.names=F,sep="\t")    
    }
  }
}



quit(save="no")

# heatmaps
library(pheatmap)
library(gplots)
library(RColorBrewer)
gt="WT"
gt="Klotho_VDR"
loc="PT"
fit=fits[["GTs"]][[loc]][[gt]]
#cols1=redgreen(colbr)
#cols2=topo.colors(colbr)

print.sig.heatmaps <- function(fit,loc,gt,colbr=75, p.max = 0.05, min.lfc = 0.0, min.Amean = -10){
  # create two heatmaps for each contrast with only the genes significant, with minimum log2 fold change and min Amean values given
  cols1=colorRampPalette(c("green",  "red"))(n = colbr)
  cols2=colorRampPalette(c("green",  "red"))(n = colbr)
  for(t.name in dimnames(fit$p.value)$Contrasts ){
    sig=dimnames(fit$p.value)[[1]][p.adjust(fit$p.value[,t.name]) < p.max & abs(fit$coefficients[,t.name]) >= min.lfc & fit$Amean >= min.Amean ]
    if (length(sig) == 0) {next()}
    exprs=GT.times.dges[["dge.GTs"]][[loc]][[gt]][["E"]][sig,]
    exprs.z=t(scale(t(exprs)))
    pheatmap(exprs.z, col=cols1, key=TRUE, symkey=FALSE, density.info="none", trace="none", scale="row", keysize=1.0, margin=c(6,6),
             clustering_distance_rows="correlation",clustering_distance_cols="correlation", fontsize_row = 10/(1+2*floor(nrow(exprs)/100)),
             filename=paste("sig_genes_heatmap",loc,gt,t.name,"scaled.pdf",sep="_"), width = 8 , height = 10)
#    pdf(file = paste("sig_genes_heatmap",loc,gt,t.name,"scaled.pdf",sep="_"), width = 8 , height = 8)
#     heatmap.2(exprs.z, col=cols1, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5/(1+floor(nrow(exprs)/100)),cexCol=0.75, 
#               scale="row", keysize=1.0, margin=c(6,6))
#    dev.off()
    pheatmap(exprs, col=cols2,key=TRUE, symkey=FALSE, density.info="none", trace="none", scale="none", keysize=1.0, margin=c(6,6),
             clustering_distance_rows="correlation",clustering_distance_cols="correlation", fontsize_row = 10/(1+2*floor(nrow(exprs)/100)),
             filename=paste("sig_genes_heatmap",loc,gt,t.name,"scaled.pdf",sep="_"), width = 8 , height = 10)
#    pdf(file = paste("sig_genes_heatmap",loc,gt,t.name,"raw.pdf",sep="_"), width = 8 , height = 8)
#     heatmap.2(exprs, col=cols2,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5/(1+floor(nrow(exprs)/100)), 
#               cexCol=0.75, scale="none", symbreaks=F, keysize=1.0, margin=c(6,6))
#    dev.off()
  }
}
t.name="t0t2"
colbr=75
colrs=colorRampPalette(c("green",  "red"))(n = colbr)
pdf(file = paste0("sig_genes",t.name))
quartz(width = 8 , height = 8)
heatmap.2(exprs.z, col=colrs, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5,cexCol=0.75, 
          scale="row", keysize=1.0, margin=c(6,6))
pheatmap(exprs.z, col=colrs, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5,cexCol=0.75, 
          scale="row", keysize=1.0, margin=c(6,6),fontsize_col=16)
heatmap.2(exprs, col=colrs,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5, cexCol=0.75, scale="none",
          symbreaks=F, keysize=1.0, margin=c(6,6))
# gsea testing
setwd("/Volumes//vetgrid01/Data/Test/")
load("wed30092015.RData")
gt="WT"
fit=fits[["GTs"]][[loc]][[gt]]
conts=colnames(fit$contrasts)
# create directory for saving files
dir.create(gt,showWarnings ='F')    
# get toptable for each contrast, merge and write to file
tts=list()
for (cont in conts){
  tts[[cont]]= topTable(fit,coef=cont,adjust.method="BH",number=nrow(fit$genes),sort.by="none")
  names(tts[[cont]])[7:12]=paste(names(tts[[cont]])[7:12],cont,sep="_")
  tts[[cont]]=cbind(rownames(tts[[cont]]),tts[[cont]])
  colnames(tts[[cont]])[1]="ENSEMBLEID"
  
}
ttmerge=merge(tts[[1]],tts[[2]][,-c(2:7)],by="ENSEMBLEID")    
ttmerge=merge(ttmerge,tts[[3]][,-c(2:7)],by="ENSEMBLEID")
rownames(ttmerge)=rownames(tts[[1]])

# setting up databases for enrichment analysis
load("MSigDB//mouse_H_v5.rdata")
load("MSigDB//mouse_c2_v5.rdata")
load("MSigDB//mouse_c3_v5.rdata")
load("MSigDB//mouse_c4_v5.rdata")
load("MSigDB//mouse_c6_v5.rdata")
load("MSigDB//mouse_c7_v5.rdata")
#KEGG
Mm.KEGG=Mm.c2[grep("KEGG",names(Mm.c2))]
Mm.KEGG=Mm.KEGG[ sapply(Mm.KEGG, function(x) (length(x) >= 4 )) ]
#REACTOME
Mm.REACT=Mm.c2[grep("REACTOME",names(Mm.c2),ignore.case = TRUE)]
Mm.REACT=Mm.REACT[ sapply(Mm.REACT, function(x) (length(x) >= 4 )) ]
#MIR
Mm.MIR=Mm.c3[grep("mir",names(Mm.c3),ignore.case = TRUE)]
Mm.MIR=Mm.MIR[ sapply(Mm.MIR, function(x) (length(x) >= 3 )) ]
#TF
Mm.TF=Mm.c3[ - grep("mir",names(Mm.c3),ignore.case = TRUE)]
Mm.TF=Mm.TF[ - grep("UNKNOWN",names(Mm.TF)) ]
Mm.TF=Mm.TF[ sapply(Mm.TF, function(x) (length(x) >= 3  )) ]

# main mart is dead use this one
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host = "jul2015.archive.ensembl.org")
datasets <- listDatasets(ensembl)
ensembl <- useDataset("mmusculus_gene_ensembl",ensembl)
Mm_GO <- getBM(attributes = c("go_id","ensembl_gene_id","namespace_1003","name_1006"),filters="ensembl_gene_id", values=rownames(all.samps$dge$genes), mart = ensembl)
Mm_GO$namespace_1003=as.factor(Mm_GO$namespace_1003)
Mm_GO=Mm_GO[Mm_GO$namespace_1003 == "biological_process",]
Mm_GO$name_1006=toupper(gsub(' ',"_",Mm_GO$name_1006))
Mm_GO$name_1006=gsub('-',"_",Mm_GO$name_1006)
Mm_GObp = split(Mm_GO$ensembl_gene_id,as.factor(Mm_GO$name_1006))
Mm_GObp = Mm_GObp[ sapply(Mm_GObp, function(x) (length(x) >= 4 )) ]



# library(GO.db)
# vals = select(GO.db, keys(GO.db, "GOID"), c("TERM", "ONTOLOGY"))
# terms=vals$TERM[ vals$ONTOLOGY == "BP"]
# terms=sub(' ',"_",terms)
# terms=toupper(terms)
#Mm.GObf=Mm.c5[which(names(Mm.c) %in% terms) ] 
#get all indices
indices=ids2indices(Mm.H,ttmerge$entrezgene)
indices=ids2indices(Mm_GObp,ttmerge$ENSEMBLEID)
gStests=c()
for( tset in names(indices)){
  pv=geneSetTest(indices[[tset]],ttmerge$t_t0t2,alternative="mixed",type="auto",ranks.only=FALSE,nsim=1000)
  gStests=append(gStests,pv)
  names(gStests)[length(gStests)]=tset
}
gStestsBH=p.adjust(gStests,method = "BH")
sort(gStestsBH[gStestsBH < 0.05])




conts <- makeContrasts(t0t2="time0-time2",t0t8="time0-time8",
                       t2t8="time2-time8",levels=GT.times.dges[["design.GTs"]][[loc]][[gt]])

indices=ids2indices(Mm.H,GT.times.dges[["dge.GTs"]][[loc]][[gt]]$genes$entrezgene)
indices=ids2indices(Mm_GObp,rownames(GT.times.dges[["dge.GTs"]][[loc]][[gt]]$genes))
a=camera(GT.times.dges[["dge.GTs"]][[loc]][[gt]],indices,GT.times.dges[["design.GTs"]][[loc]][[gt]],contrast=c(1,-1,0))
#interGeneCorrelation(GT.times.dges[["dge.GTs"]][[loc]][[gt]]$E,GT.times.dges[["design.GTs"]][[loc]][[gt]])
a
b=romer(GT.times.dges[["dge.GTs"]][[loc]][[gt]],indices,design=GT.times.dges[["design.GTs"]][[loc]][[gt]],contrast=c(1,-1,0))
b

br=mroast(GT.times.dges[["dge.GTs"]][[loc]][[gt]],indices,design=GT.times.dges[["design.GTs"]][[loc]][[gt]],contrast=c(1,-1,0))
br=mroast(GT.times.dges[["dge.GTs"]][[loc]][[gt]],indices,design=GT.times.dges[["design.GTs"]][[loc]][[gt]],contrast=c(0,1,-1))
br


# trying out GOexpress ... not very convincing for complex experiments, no weights, no nothing :(  

library("GOexpress")
require(Biobase)
# convert to Expression Set
my_set = new("ExpressionSet", exprs=GT.times.dges[["dge.GTs"]][[loc]][[gt]]$E, 
             phenoData=AnnotatedDataFrame(GT.times.dges[["sm.GTs"]][[loc]][[gt]]))
names(my_set)
data(AlvMac)
str(AlvMac)
AlvMac$Treatment
my_set$group
# main mart is dead use this one
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host = "jul2015.archive.ensembl.org")
datasets <- listDatasets(ensembl)
ensembl <- useDataset("mmusculus_gene_ensembl",ensembl)
Mm_GO <- getBM(attributes = c("go_id","ensembl_gene_id","namespace_1003","name_1006"),filters="ensembl_gene_id", values=rownames(all.samps$dge$genes), mart = ensembl)
Mm_GO$namespace_1003=as.factor(Mm_GO$namespace_1003)
Mm_GO=Mm_GO[Mm_GO$namespace_1003 == "biological_process",]
Mm_GO$name_1006=toupper(gsub(' ',"_",Mm_GO$name_1006))
Mm_GO$name_1006=gsub('-',"_",Mm_GO$name_1006)
Mm_GObf = split(Mm_GO$ensembl_gene_id,as.factor(Mm_GO$name_1006))
Mm_GObf = Mm_GObf[ sapply(Mm_GObf, function(x) (length(x) >= 4 )) ]

my_set_anal = GO_analyse(eSet = my_set, f = "group", GO_genes=GO_genes, all_GO=all_GO, all_genes=all_genes)
my_set_anal.pVal = pValue_GO(result=my_set_anal, N=100)
BP.5 <- subset_scores(result = my_set_anal.pVal,namespace = "biological_process", total = 5, p.val=0.05)

expression_plot(gene_id = "ENSMUSG00000025902", result = my_set_anal, eSet=my_set, f="group", x_var = "group", 
                title.size=1.5, legend.title.size=10, legend.text.size=10, legend.key.size=15)

tail(BP.5$GO)
heatmap_GO(go_id = "GO:0010803", result = BP.5, eSet=my_set, cexRow=0.4, cexCol=1, cex.main=1, main.Lsplit=30)
library("biomaRt", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
goids = getBM(attributes=c("ensembl_gene_id","go_id"), filters="ensembl_gene_id", values=rownames(GT.times.dges[["dge.GTs"]][[loc]][[gt]]$E),mart=ensembl)
goids = getBM(attributes=c("ensembl_gene_id","go_id"), filters="ensembl_gene_id", values=rownames(GT.times.dges[["dge.GTs"]][[loc]][[gt]]$E),mart=ensembl)
GO_genes <- getBM(attributes = c("ensembl_gene_id", "go_id"), values=rownames(GT.times.dges[["dge.GTs"]][[loc]][[gt]]$E), mart = ensembl)
colnames(GO_genes)[1]="gene_id"
GO_genes <- GO_genes[GO_genes$go_id != "", ]
GO_genes <- GO_genes[GO_genes$gene_id != "", ]
all_GO <- getBM(attributes = c("go_id", "name_1006", "namespace_1003"), values=rownames(GT.times.dges[["dge.GTs"]][[loc]][[gt]]$E),mart = ensembl)
all_GO <- all_GO[all_GO$go_id != "", ]
all_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"), values=rownames(GT.times.dges[["dge.GTs"]][[loc]][[gt]]$E),mart = ensembl)
colnames(all_genes)[1]="gene_id"
save
df <- data.frame(Expression = exprs(eSet)[gene_id, ], Factor = pData(eSet)[, f], X = pData(eSet)[, x_var])

title="blub"

ggplot(df) + geom_smooth(aes(x = X, y = Expression, group = Factor, color = Factor, fill = Factor)) 
df

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


enrDBs = unique(sapply(enr.cats, function (x) x[1] ))
enrlocs = unique(sapply(enr.cats, function (x) x[2] ))
enrgtps = unique(sapply(enr.cats, function (x) x[3] ))

setwd("/Users/lukasendler/Documents/Projects/Mouse_Kidney_Andrukhova_Erben/Andrukhova_Erben_Results/")
path.name="Highest1000/Enrichments"
for(cat in enrDBs){
  for(loc in enrlocs) {
    for(gt in enrgtps) {
      res=list()
      j=1
      res.names=c()
      for(i in which(sapply(enr.cats, function (x)  x[1] == cat && x[2] == loc && x[3] ==gt  )) ){
        tab=read.table(paste(path.name,enr.fn[i],sep="/"),header=TRUE)
        tab=tab[tab$FDR <= 0.05,]
        if (nrow(tab) < 1) next
        res[[j]] = tab
        res.names=append(res.names,paste(c(cat,loc,gt,"cont",enr.cats[[i]][5:length(enr.cats[[i]])]),sep = "_",collapse ="_"))
        j=j+1
      }
      if (length(res) < 1) next
      names(res)=res.names
      xlsbname=paste("enr",cat,loc,gt,"FDR0.05.xls",sep="_")
      print(xlsbname)
      print(res.names)
      
      WriteXLS(res, ExcelFileName = paste0(path.name,"/",xlsbname), SheetNames = names(res), AdjWidth = FALSE, BoldHeaderRow = TRUE,FreezeRow = 1, FreezeCol = 0)
    }
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

get_de_lists <- function (fit,outspec) {
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
  rownames(ttmerge)=ttmerge$ENSEMBLEID
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
    loc.list=append(loc.list,get_de_lists(fit,outspec))
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
    loc.list=append(loc.list,get_de_lists(fit,outspec))
  }
  xlsbname=paste("de_genes_timepoint_",time,".xls",sep="_")
  WriteXLS(loc.list, ExcelFileName = paste0(path.name,"/",timepoint,"/",xlsbname), SheetNames = names(loc.list), AdjWidth = TRUE, BoldHeaderRow = TRUE,FreezeRow = 1, FreezeCol = 1)
  
}

# All genes DE independent of GT but over time
# create DT PT
dge.locs=list( list())
mincount=25
for (j in levels(all.samps$sample.mat$loc)){
  keep.loc <- all.samps$sample.mat$loc == j 
  sample.mat.loc = droplevels(all.samps$sample.mat[keep.loc,])
  rkeep.loc = ! ( apply(all.samps$dge[[1]][,keep.loc],1,max,na.rm = TRUE) < mincount)
  #rkeep.loc = ! (rowSums(dge[[1]][,keep.loc]) < mincount)
  dge.locs[['dge']][[j]] <- DGEList(counts=all.samps$dge[[1]][rkeep.loc,keep.loc], group=sample.mat.loc$group, 
                                    genes=all.samps$dge[[3]][rkeep.loc,])
  dge.locs[['sm']][[j]] <- sample.mat.loc
  dge.locs[['design']][[j]] <- model.matrix( ~ 0 + GT + time,data=sample.mat.loc)
  # normalisation factors
  dge.locs[["dge"]][[j]] <- calcNormFactors(dge.locs[["dge"]][[j]])
  filename=paste("voom_qw_loc",j,"pdf",sep = ".")
  pdf(file=filename)
  # calculate with quality weights
  dge.locs[["dge"]][[j]] <- voomWithQualityWeights(dge.locs[["dge"]][[j]], design=dge.locs[['design']][[j]], normalization="none", plot=TRUE)
  dev.off()
}
#save.image("thu17122015.RData")

# define the model matrices for the different GTs and times (already correctly done in the function get.GTandTime.dges)
for(loc in names(dge.locs[["dge"]]) ){
  dge.locs[["design.GTs"]][[loc]] <- model.matrix( ~  0 + time, data = dge.locs[["sm"]][[loc]])
  dge.locs[["design.times"]][[loc]] <- model.matrix( ~  0 + GT, data = dge.locs[["sm"]][[loc]])
}


fits.loc=list(list())
for (loc in names(dge.locs[["design.times"]]) ) {
  conts <- makeContrasts(t0t2="time0-time2",t0t8="time0-time8",
                         t2t8="time2-time8",levels=dge.locs[["design.GTs"]][[loc]])
  fit.GT <- lmFit(dge.locs[["dge"]][[loc]], design=dge.locs[["design.GTs"]][[loc]])
  fit.GT <- contrasts.fit(fit.GT, conts)
  fit.GT <- eBayes(fit.GT)
  fits.loc[["GTs"]][[loc]] <- fit.GT
  conts <- makeContrasts(WT_FGF23="GTWT-GTFGF23_VDR",WT_Klotho="GTWT-GTKlotho_VDR",
                         FGF23_Klotho="GTFGF23_VDR-GTKlotho_VDR",levels=dge.locs[["design.times"]][[loc]])
  fit.time <- lmFit(dge.locs[["dge"]][[loc]], design=dge.locs[["design.times"]][[loc]])
  fit.time <- contrasts.fit(fit.time, conts)
  fit.time <- eBayes(fit.time)
  fits.loc[["times"]][[loc]] <- fit.time
}

loc.list=list()
for (loc in  names(fits.loc[["GTs"]]) ){
  outspec=paste0("loc_",loc,"_times")
  fit=fits.loc[["times"]][[loc]]
  conts=colnames(fit$contrasts)
  loc.list=append(loc.list,get_de_lists(fit,outspec))
}
xlsbname=paste("de_genes_alltimes.xls",sep="_")
WriteXLS(loc.list, ExcelFileName = paste0(path.name,"/",xlsbname), SheetNames = names(loc.list), AdjWidth = TRUE, BoldHeaderRow = TRUE,FreezeRow = 1, FreezeCol = 5)
loc.list=list()
for (loc in  names(fits.loc[["GTs"]]) ){
  outspec=paste0("loc_",loc,"_GTs")
  fit=fits.loc[["GTs"]][[loc]]
  conts=colnames(fit$contrasts)
  loc.list=append(loc.list,get_de_lists(fit,outspec))
}
xlsbname=paste("de_genes_allGTs.xls",sep="_")
WriteXLS(loc.list, ExcelFileName = paste0(path.name,"/",xlsbname), SheetNames = names(loc.list), AdjWidth = TRUE, BoldHeaderRow = TRUE,FreezeRow = 1, FreezeCol = 5)
library("gplots")
print.sig.heatmaps <- function(fit,loc,gt,colbr=75, p.max = 0.05, min.lfc = 0.0, min.Amean = -10){
  # create two heatmaps for each contrast with only the genes significant, with minimum log2 fold change and min Amean values given
  cols1=colorRampPalette(c("green",  "red"))(n = colbr)
  cols2=colorRampPalette(c("green",  "red"))(n = colbr)
  for(t.name in dimnames(fit$p.value)$Contrasts ){
    sig=dimnames(fit$p.value)[[1]][p.adjust(fit$p.value[,t.name]) < p.max & abs(fit$coefficients[,t.name]) >= min.lfc & fit$Amean >= min.Amean ]
    if (length(sig) == 0) {next()}
    exprs=GT.times.dges[["dge.GTs"]][[loc]][[gt]][["E"]][sig,]
    exprs.z=t(scale(t(exprs)))
    pheatmap(exprs.z, col=cols1, key=TRUE, symkey=FALSE, density.info="none", trace="none", scale="row", keysize=1.0, margin=c(6,6),
             clustering_distance_rows="correlation",clustering_distance_cols="correlation", fontsize_row = 10/(1+2*floor(nrow(exprs)/100)),
             filename=paste(bn,t.name,"scaled.pdf",sep="_"), width = 8 , height = 10)
    pheatmap(exprs, col=cols2,key=TRUE, symkey=FALSE, density.info="none", trace="none", scale="none", keysize=1.0, margin=c(6,6),
             clustering_distance_rows="correlation",clustering_distance_cols="correlation", fontsize_row = 10/(1+2*floor(nrow(exprs)/100)),
             filename=paste(bn,t.name,"raw.pdf",sep="_"), width = 8 , height = 10)
  }
}

print.sig.heatmaps.2 <- function(fit,loc,type,colbr=75,dir=".", p.max = 0.05, min.lfc = 0.0, min.Amean = -10){
  cols1=colorRampPalette(c("green",  "red"))(n = colbr)
  cols2=colorRampPalette(c("green",  "red"))(n = colbr)
  bn=paste0(dir,"/sig_genes_heatmap_",type,"_location_",loc,".pdf")
  for(t.name in dimnames(fit$p.value)$Contrasts ){
    sig=dimnames(fit$p.value)[[1]][p.adjust(fit$p.value[,t.name]) < p.max & abs(fit$coefficients[,t.name]) >= min.lfc & fit$Amean >= min.Amean ]
    if (length(sig) == 0) {next()}
    exprs=dge.locs[["dge"]][[loc]][["E"]][sig,]
    exprs.z=t(scale(t(exprs)))
    pheatmap(exprs.z, col=cols1, key=TRUE, symkey=FALSE, density.info="none", trace="none", scale="row", keysize=1.0, margin=c(6,6),
             clustering_distance_rows="correlation",clustering_distance_cols="correlation", fontsize_row = 10/(1+2*floor(nrow(exprs)/100)),
             filename=paste(bn,t.name,"scaled.pdf",sep="_"), width = 8 , height = 10)
    pheatmap(exprs, col=cols2,key=TRUE, symkey=FALSE, density.info="none", trace="none", scale="none", keysize=1.0, margin=c(6,6),
             clustering_distance_rows="correlation",clustering_distance_cols="correlation", fontsize_row = 10/(1+2*floor(nrow(exprs)/100)),
             filename=paste(bn,t.name,"raw.pdf",sep="_"), width = 8 , height = 10)
  }
}

for (loc in  names(fits.loc[["GTs"]]) ){
  fit=fits.loc[["GTs"]][[loc]]
  print.sig.heatmaps.2(fit,loc,"GTs",colbr=75,dir="Highest1000/")
  fit=fits.loc[["times"]][[loc]]
  print.sig.heatmaps.2(fit,loc,"times",colbr=75,dir="Highest1000/")
}


# Easy Go and Kegg enrichment analysis
#library("Go.db")
source("http://bioconductor.org/biocLite.R")
library(pathview)
library(gage)
library(gageData)
data(kegg.sets.mm)
data(go.sets.mm)
data(sigmet.idx.mm)
data(go.subs.mm)
kegg.mmu <- kegg.gsets(species = "mmu", id.type = "entrez")
data(korg)
data(bods)

conts=colnames(fit$contrasts)
# create directory for saving files
dir.create(gt,showWarnings ='F') 

get_ttmerge <- function(fit){
  tts=list()
  conts=colnames(fit$contrasts)
  for (cont in conts){
    tts[[cont]]= topTable(fit,coef=cont,adjust.method="BH",number=nrow(fit$genes),sort.by="none")
    names(tts[[cont]])[7:12]=paste(names(tts[[cont]])[7:12],cont,sep="_")
    tts[[cont]]=cbind(rownames(tts[[cont]]),tts[[cont]])
    colnames(tts[[cont]])[1]="ENSEMBLEID"      
  }
  ttmerge=merge(tts[[1]],tts[[2]][,-c(2:7)],by="ENSEMBLEID")    
  ttmerge=merge(ttmerge,tts[[3]][,-c(2:7)],by="ENSEMBLEID")
  rownames(ttmerge)=ttmerge$ENSEMBLEID
  return(ttmerge)
}
# get toptable for each contrast, merge and write to file

cont="WT_Klotho"
cont.fc=ttmerge[,paste("logFC",cont,sep="_")]
names(cont.fc)=ttmerge$entrezgene
fc.kegg.p <- gage(cont.fc, gsets = kegg.mmu$kg.sets, ref = NULL, samp = NULL, use.fold= TRUE,same.dir = FALSE)
sel <- fc.kegg.p$greater[, "q.val"] < 0.05 & !is.na(fc.kegg.p$greater[, "q.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]
path.ids2 <- sapply(path.ids, function(x) unlist(strsplit(x," "))[1])
head(path.ids2)
out.suffix="limma.test.pathview"
pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = cont.fc, pathway.id = pid,species = "mmu", out.suffix=out.suffix))
pathview(gene.data = cont.fc, pathway.id = "mmu04668",species = "mmu", out.suffix=out.suffix)

library("topGO")
cont.q=ttmerge[,paste("adj.P.Val",cont,sep="_")]
cont.q = factor(as.integer(cont.q < 0.05))
cont.fc=ttmerge[,paste("logFC",cont,sep="_")]
names(cont.fc)=rownames(ttmerge)

names(cont.q)=rownames(ttmerge)
GOdata <- new("topGOdata", ontology = "BP", allGenes = cont.q , annotationFun = annFUN.org, mapping = "org.Mm.eg.db", ID = "ensembl")
ensemble.go=split(Mm_GO$ensembl_gene_id,Mm_GO$go_id)
fc.gage.p <- gage(cont.fc, gsets = ensemble.go, ref = NULL, samp = NULL, use.fold= TRUE,same.dir = FALSE,set.size = c(5,300))
fc.gage.pv = fc.gage.p$greater[,"p.val"]
fc.gage.pv[is.na(fc.gage.pv)] = 1.0
gage.topGo = new("topGOresult", description = "test top go", score=fc.gage.pv , testName="gage", algorithm="gage")
score(gage.topGo)

resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
pvalFis <- score(resultFis)
hist(pvalFis)
allRes <- GenTable(GOdata, GAGE = gage.topGo, topNodes = 20)
allRes
showSigOfNodes(GOdata, score(gage.topGo), firstSigNodes = 25, useInfo = 'all')
printGraph(GOdata, gage.topGo , firstSigNodes =5, useInfo = 'def',fn.prefix = "sampleFile",pdfSW = TRUE)
data(results.tGO)
s <- score(resultFisher)


# kegg.code = "mmu" entrez ids used
#new("topGOresult", description, score, testName, algorithm, geneData)
fit.go <- goana(fit,species="Mm",geneid = "entrezgene",coef=cont)
topGO(fit.go,ontology = c("BP"),n=100)
fit.kegg <- kegga(fit,species="Mm",geneid = "entrezgene",coef=cont,plot=TRUE)
topKEGG(fit.kegg,n=100)

library(EnrichmentBrowser)

# obtain grn for mouse
pwys <- download.kegg.pathways("mmu")
mmu.grn <- compile.grn.from.kegg(pwys)

kegg.gs <- get.kegg.genesets("mmu")
go.gs <- get.go.genesets(org="mmu", onto="BP", mode="GO.db")
all.genes=all.samps$dge$genes$entrezgene
go.gs = sapply( go.gs, function(x) x[x %in% all.genes] )
go.gs = go.gs[sapply( go.gs, function(x) length(x) > 5 )]
go.gs = go.gs[sapply( go.gs, function(x) length(x) < 400 )]
kegg.gs = sapply( kegg.gs, function(x) x[x %in% all.genes] )
kegg.gs = kegg.gs[sapply( kegg.gs, function(x) length(x) > 5 )]
kegg.gs = kegg.gs[sapply( kegg.gs, function(x) length(x) < 300 )]
# getting the pathways only
# library(KEGGREST)
# kegg.pw = keggList("pathway","mmu")
# kegg.pw=gsub("path:","",names(kegg.pw))
# k.dis = keggList("disease")
# k.pws = sapply(names(kegg.gs), function (x) strsplit(x,"_")[[1]][1] %in%  kegg.pw )



# create reports for enrichment of KeGG and GO using enrichementbrowser
fit <- fits[["GTs"]][[loc]][[gt]]
dge <- GT.times.dges[["dge.GTs"]][[loc]][[gt]]
conts=colnames(fit$contrasts)
cont="t0t2"
ttmerge <- get_ttmerge(fit)
featureData <- ttmerge[,c("logFC_t0t2","adj.P.Val_t0t2","t_t0t2")]
colnames(featureData) = c("FC","ADJ.PVAL","limma.STAT")

featureData <- new("AnnotatedDataFrame", data=featureData)
exprs=as.matrix(dge)
rownames(exprs) = rownames(dge)
exprs=exprs[as.character(ttmerge$ENSEMBLEID),]
pdata=as.data.frame(dge$targets)
colnames(pdata)=c("GROUP","lib.size","norm.factors")
pdata = new("AnnotatedDataFrame", data=pdata)
eset<-new("ExpressionSet", exprs=exprs, featureData = featureData, phenoData = pdata)
annotation(eset) <- "org.Mm.eg"
eset <- map.ids(eset, org="mmu", from="ENSEMBL", to="ENTREZID")

sbea.res <- sbea( method = "ora", eset=eset, gs=kegg.gs , alpha = 0.05, perm = 0, padj.method = "BH")
sbea.go.ora <- sbea( method = "ora", eset=eset, gs=go.gs , alpha = 0.05, perm = 0, padj.method = "BH")
#sbea.go.gsea <- sbea( method = "gsea", eset=eset, gs=go.gs , alpha = 0.05, perm = 1000, padj.method = "BH")
gs.ranking(sbea.res)
gs.ranking(sbea.go)
config.ebrowser( "OUTDIR.DEFAULT", value = getwd() )
ea.browse(sbea.go, html.only = TRUE)
ea.browse(sbea.gsea, html.only = TRUE)
pData(eset)


sbea.ora <- sbea( method = "ora", eset=eset, gs=kegg.gs , alpha = 0.05, perm = 0, padj.method = "BH")
sbea.gsenr <- sbea( method = use.gsenr , eset=eset, gs=kegg.gs , alpha = 0.05, perm = 1000, padj.method = "BH")
sbea.gsenr$method="GSENR"
sbea.gsenr$res.tbl$NR.SIG.GENES = get.nsig(sbea.gsenr)
sbea.gsenr$res.tbl$NR.GENES = get.ngenes(sbea.gsenr)


# save.image("tue26012016.RData")

use.camera <- function(eset , gs, alpha = 0.05, perm = 0 ) {
  cont.list = list(t0t2=c(1,-1,0),t0t8=c(1,0,-1),t2t8=c(0,1,-1), WT_FGF23=c(1,-1,0),WT_Klotho=c(1,0,-1), FGF23_Klotho=c(0,1,-1))
  #translate into index vectors
  gs.ind = ids2indices(gs,dge$genes$entrezgene, remove.empty=TRUE)
  gStests=as.data.frame(camera(dge,gs.ind,dge$design,contrast=cont.list[[cont]],inter.gene.cor=0.05))
  #interGeneCorrelation(dge$E, dge$design)
  gStests.p=gStests$FDR
  names(gStests.p)=rownames(gStests)
  return(gStests.p)
}

use.mroast <- function(eset , gs, alpha = 0.05, perm = 0) {
  cont.list = list(t0t2=c(1,-1,0),t0t8=c(1,0,-1),t2t8=c(0,1,-1), WT_FGF23=c(1,-1,0),WT_Klotho=c(1,0,-1), FGF23_Klotho=c(0,1,-1))
  #translate into index vectors
  gs.ind = ids2indices(gs,dge$genes$entrezgene, remove.empty=TRUE)
  gStests=as.data.frame(mroast(dge,gs.ind,dge$design,contrast=cont.list[[cont]], nrot = perm ))
  gStests.p=gStests$FDR.Mixed
  names(gStests.p)=rownames(gStests)
  return(gStests.p)
  
}

use.fry <- function(eset , gs, alpha = 0.05, perm = 0) {
  cont.list = list(t0t2=c(1,-1,0),t0t8=c(1,0,-1),t2t8=c(0,1,-1), WT_FGF23=c(1,-1,0),WT_Klotho=c(1,0,-1), FGF23_Klotho=c(0,1,-1))
  #translate into index vectors
  gs.ind = ids2indices(gs,dge$genes$entrezgene, remove.empty=TRUE)
  gStests=as.data.frame(fry(dge,gs.ind,dge$design,contrast=cont.list[[cont]]))
  gStests.p=gStests$FDR.Mixed
  names(gStests.p)=rownames(gStests)
  return(gStests.p)
  
}

use.romer <- function(eset , gs, alpha = 0.05, perm = 0) {
  cont.list = list(t0t2=c(1,-1,0),t0t8=c(1,0,-1),t2t8=c(0,1,-1), WT_FGF23=c(1,-1,0),WT_Klotho=c(1,0,-1), FGF23_Klotho=c(0,1,-1))
  #translate into index vectors
  gs.ind = ids2indices(gs,dge$genes$entrezgene, remove.empty=TRUE)
  gStests=as.data.frame(romer(dge,gs.ind,dge$design,contrast=cont.list[[cont]], nrot = perm ))
  gStests.p=gStests$Mixed
  names(gStests.p)=rownames(gStests)
  return(gStests.p)
}


use.goana <- function(eset , gs, alpha = 0.05, perm = 0) {
  genes = rownames(fData(eset))
  p.values = fData(eset)$ADJ.PVAL
  names(p.values) = rownames(fData(eset))
  #gs = sapply( gs, function(x) x[ x %in% genes] )
  gs.melt=melt(gs)
  gs.melt[,1]=as.character(gs.melt[,1])
  gs.melt[,2]=as.factor(gs.melt[,2])
  gs.ids=unique(gs.melt[,2])
  gStests=goana(genes[p.values < 0.05],species="Mm",geneid = genes,gene.pathway = gs.melt,pathway.names=cbind(gs.ids,gs.ids))
  #gStestsBH=p.adjust(gStests,method = "BH")
  gStests.p=gStests$P.DE
  names(gStests.p)=rownames(gStests)
  return(gStests.p)
}

use.kegga <- function(eset , gs, alpha = 0.05, perm = 9999) {
  genes = rownames(fData(eset))
  p.values = fData(eset)$ADJ.PVAL
  names(p.values) = rownames(fData(eset))
  #gs = sapply( gs, function(x) x[ x %in% genes] )
  gs.melt=melt(gs)
  gs.melt[,1]=as.character(gs.melt[,1])
  gs.melt[,2]=as.factor(gs.melt[,2])
  gs.ids=unique(gs.melt[,2])
  gStests=kegga(genes[p.values < 0.05],species="Mm",geneid = genes,gene.pathway = gs.melt,pathway.names=cbind(gs.ids,gs.ids))
  #gStestsBH=p.adjust(gStests,method = "BH")
  gStests.p=gStests$P.DE
  names(gStests.p)=rownames(gStests)
  return(gStests.p)
}

use.gage <- function(eset , gs, alpha = 0.05, perm = 9999) {
  genes = rownames(fData(eset))
  fc.values = fData(eset)$FC
  names(fc.values) = rownames(fData(eset))
  fc.gage.p <- gage(fc.values, gsets = gs, ref = NULL, samp = NULL, use.fold= TRUE,same.dir = FALSE,set.size = c(5,400))
  fc.gage.p <- fc.gage.p$greater[,"q.val"]
  fc.gage.p[is.na(fc.gage.p)] = 1.0
  return(fc.gage.p)
}



path.name="Highest1000"
for (gt in names(GT.times.dges[["dge.GTs"]][["DT"]])){
  loc.list=list()
  for (loc in names(GT.times.dges[["dge.times"]]) ){
    print(paste("going through",loc,gt))
    outspec=paste0("loc_",loc,"_gt_",gt)
    fit=fits[["GTs"]][[loc]][[gt]]
    loc.list=append(loc.list,get_de_lists(fit,outspec,dge.sub = GT.times.dges[["dge.GTs"]][[loc]][[gt]]))
  }
  xlsbname=paste("de_genes_genotype",gt,".xls",sep="_")
  WriteXLS("loc.list", ExcelFileName = paste0(path.name,"/",gt,"/",xlsbname), SheetNames = names(loc.list), AdjWidth = TRUE, BoldHeaderRow = TRUE,FreezeRow = 1, FreezeCol = 5)
}
for (time in names(GT.times.dges[["dge.times"]][["DT"]]) ) {
  loc.list=list()
  timepoint=paste0("timepoint_",time)
  for (loc in names(GT.times.dges[["dge.times"]]) ){
    print(paste("going through",loc,time))
    outspec=paste0("loc_",loc,"_time_",time)
    fit=fits[["times"]][[loc]][[time]]
    conts=colnames(fit$contrasts)
    loc.list=append(loc.list,get_de_lists(fit,outspec,dge.sub = GT.times.dges[["dge.times"]][[loc]][[time]]))
  }
  xlsbname=paste("de_genes_timepoint_",time,".xls",sep="_")
  WriteXLS("loc.list", ExcelFileName = paste0(path.name,"/",timepoint,"/",xlsbname), SheetNames = names(loc.list), AdjWidth = TRUE, BoldHeaderRow = TRUE,FreezeRow = 1, FreezeCol = 5)
}

loc.list=list()
for (loc in  names(fits.loc[["GTs"]]) ){
  outspec=paste0("loc_",loc,"_times")
  fit=fits.loc[["times"]][[loc]]
  conts=colnames(fit$contrasts)
  loc.list=append(loc.list,get_de_lists(fit,outspec,dge.sub = dge.locs[["dge"]][[loc]]))
}
xlsbname=paste("de_genes_alltimes.xls",sep="_")
WriteXLS("loc.list", ExcelFileName = paste0(path.name,"/",xlsbname), SheetNames = names(loc.list), AdjWidth = TRUE, BoldHeaderRow = TRUE,FreezeRow = 1, FreezeCol = 5)
loc.list=list()
for (loc in  names(fits.loc[["GTs"]]) ){
  outspec=paste0("loc_",loc,"_GTs")
  fit=fits.loc[["GTs"]][[loc]]
  conts=colnames(fit$contrasts)
  loc.list=append(loc.list,get_de_lists(fit,outspec,dge.sub = dge.locs[["dge"]][[loc]]))
}
xlsbname=paste("de_genes_allGTs.xls",sep="_")
WriteXLS("loc.list", ExcelFileName = paste0(path.name,"/",xlsbname), SheetNames = names(loc.list), AdjWidth = TRUE, BoldHeaderRow = TRUE,FreezeRow = 1, FreezeCol = 5)








sbea.safe <- sbea( method = "safe" , eset=eset, gs=kegg.gs , alpha = 0.05, perm = 0, padj.method = "BH")
sbea.gage <- sbea( method = use.gage , eset=eset, gs=kegg.gs , alpha = 0.05, perm = 0, padj.method = "none")
sbea.gage$method="GAGE"
sbea.camera <- sbea( method = use.camera , eset=eset, gs=kegg.gs , alpha = 0.05, perm = 0, padj.method = "none")
sbea.camera$method="CAMERA"
# sbea.romer <- sbea( method = use.romer , eset=eset, gs=kegg.gs , alpha = 0.05, perm = 10000, padj.method = "none")
# sbea.romer$method="ROMER"

go.ora <- sbea( method = "ora", eset=eset, gs=go.gs , alpha = 0.05, perm = 0, padj.method = "BH")
#go.safe <- sbea( method = "safe" , eset=eset, gs=go.gs , alpha = 0.05, perm = 0, padj.method = "BH")
go.camera <- sbea( method = use.camera , eset=eset, gs=go.gs , alpha = 0.05, perm = 0, padj.method = "none")
go.camera$method="CAMERA"
# go.romer <- sbea( method = use.romer , eset=eset, gs=go.gs , alpha = 0.05, perm = 10000, padj.method = "none")
# go.romer$method="ROMER"
#go.roast <- sbea( method = use.mroast , eset=eset, gs=go.gs , alpha = 0.05, perm = 10000, padj.method = "none")
#go.roast$method="ROAST"
start.time <- Sys.time()
go.gsenr <- sbea( method = use.gsenr , eset=eset, gs=go.gs , alpha = 0.05, perm = 10000, padj.method = "BH")
end.time <- Sys.time()
time.taken <- end.time - start.time
go.gsenr$method="GSENR"
go.gage <- sbea( method = use.gage , eset=eset, gs=go.gs , alpha = 0.05, perm = 0, padj.method = "none")
go.gage$method="GAGE"
# go.fry <- sbea( method = use.fry , eset=eset, gs=go.gs , alpha = 0.05, perm = 0, padj.method = "none")
# go.fry$method="fry"
gs.ranking(sbea.gsenr)
ea.browse(sbea.gsenr, html.only = TRUE)

res <- list(sbea.gsenr,sbea.ora)
comb.res <- comb.ea.results(res)
cor(as.matrix(comb.res$res.tbl[,grep(pattern = "PVAL",x = names(comb.res$res.tbl))]),method = "spearman")
#library(car)
scatterplotMatrix(comb.res$res.tbl[,grep(pattern = "PVAL",x = names(comb.res$res.tbl))], diagonal="histogram",nclass=20)
intersect(gs.ranking(sbea.gage)$GENE.SET, gs.ranking(sbea.ora)$GENE.SET)
setdiff(gs.ranking(sbea.gsenr)$GENE.SET, gs.ranking(sbea.ora)$GENE.SET)
gs.ranking(comb.res)

res.go <- list(go.camera,go.gsenr,go.ora,go.gage)
res.go <- list(go.gsenr,go.ora)
comb.res.go <- comb.ea.results(res.go)
cor(as.matrix(comb.res.go$res.tbl[,grep(pattern = "PVAL",x = names(comb.res.go$res.tbl))]),method = "spearman")
#library(car)
scatterplotMatrix(comb.res.go$res.tbl[,grep(pattern = "PVAL",x = names(comb.res.go$res.tbl))], diagonal="histogram",nclass=20)
intersect(gs.ranking(go.ora)$GENE.SET, gs.ranking(go.gsenr)$GENE.SET)
union(gs.ranking(go.ora)$GENE.SET, gs.ranking(go.gsenr)$GENE.SET)

nbea.res <- nbea(method="ggea", eset=eset, gs=kegg.gs, grn=mmu.grn, perm = 10000, padj.method = "BH")
gs.ranking(nbea.res)

res.nr <- list(sbea.gsenr,nbea.res)
comb.nr <- comb.ea.results(res.nr)
all = gs.ranking(comb.nr, signif.only = FALSE)
sig.sets=all[rowMin(as.matrix(all[,grep("PVAL",names(all))])) < 0.05  , ]
annotation(comb.nr$eset) = "mmu"
ea.browse(comb.nr,graph.view=mmu.grn, html.only = TRUE, nr.show = nrow(sig.sets))

# use other gene sets
gs.tf = Mm.Migs[["TF"]]
gs.react = Mm.Migs[["REACT"]]
gs.mir = Mm.Migs[["MIR"]]
names(gs.tf) = sapply(names(gs.tf), function (x) gsub("_","-",x) )
names(gs.react) = sapply(names(gs.react), function (x) gsub("_","-",x) )
names(gs.tf) = sapply(names(gs.tf), function (x) gsub("_","-",x) )
hh.ora <- sbea( method = "ora", eset=eset, gs=gs.react , alpha = 0.05, perm = 0, padj.method = "BH")
gs.ranking(hh.ora, signif.only = TRUE)
annotation(hh.ora$eset) = "mmu"
ea.browse(hh.ora, html.only = TRUE)


gs.ranking(comb.res, signif.only = TRUE)
gs.ranking(comb.res.go, signif.only = TRUE)
ea.browse(comb.res, html.only = TRUE)
ea.browse(comb.res,graph.view=mmu.grn, html.only = TRUE)
cor(as.matrix(comb.res.go$res.tbl[comb.res.go$res.tbl$ROMER.PVAL < 1,grep(pattern = "PVAL",x = names(comb.res.go$res.tbl))]),method = "spearman")
cor(comb.res.go$res.tbl$CAMERA.PVAL,comb.res.go$res.tbl$ROAST.PVAL, method = "spearman")

sbea.gage <- sbea( method = use.gage , eset=eset, gs=kegg.gs , alpha = 0.05, perm = 0, padj.method = "none")
sbea.gage$method="GAGE"
sbea.gsea <- sbea( method = "gsea" , eset=eset, gs=kegg.gs , alpha = 0.05, perm = 9999, padj.method = "BH")
sbea.samgs <- sbea( method = "samgs" , eset=eset, gs=kegg.gs , alpha = 0.05, perm = 1000, padj.method = "BH")

sb.gage=gs.ranking(sbea.gage)
sbea.ora <- sbea( method = "ora" , eset=eset, gs=kegg.gs , alpha = 0.05, perm = 1000, padj.method = "BH")
nbea.res <- nbea(method="ggea", eset=eset, gs=kegg.gs, grn=mmu.grn)
gs.ranking(nbea.res)
sb.ora=gs.ranking(sbea.ora)
merge()
res <- list(sbea.gage,nbea.res)
comb.res <- comb.ea.results(res)
gs.ranking(comb.res)
ea.browse(comb.res,graph.view=mmu.grn, html.only = TRUE)

save.image("tue26012016.RData")

undebug(ggea.graph)
ggea.graph(gs=kegg.gs[["mmu04512_ECM-receptor_interaction"]],grn=mmu.grn, eset=eset)


sbea.go <- sbea( method = use.gsenr , eset=eset, gs=go.gs , alpha = 0.05, perm = 9999, padj.method = "BH")
gs.ranking(sbea.go)
sbea.methods()
# transforming tab delimited files into 
for(i in )
b = goana(fit,species="Mm",geneid = "entrezgene",coef=cont,plot=TRUE)
a= goana(genes[p.values < 0.05],species="Mm",geneid = genes,plot=TRUE)
a[p.adjust(a$P.DE,method = "BH") < 0.05,]


# clusterProfiler and Dose
library(DOSE)
library(clusterProfiler)

genes = rownames(fData(eset))
fc.values = fData(eset)$FC
names(fc.values) = rownames(fData(eset))

kk2 <- gseKEGG(geneList     = sort(fc.values,decreasing=T),
               organism     = "mouse",
               nPerm        = 1000,
               minGSSize    = 5,
               pvalueCutoff = 0.05,
               verbose      = TRUE,
               use_internal_data = FALSE)

ego2 <- gseGO(geneList     = sort(fc.values,decreasing=T),
              organism     = "mouse",
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 5,
              pvalueCutoff = 0.01,
              verbose      = FALSE)


enrichMap(kk2)
cnetplot(kk2, categorySize="pvalue", foldChange=fc.values)

eset <- make.example.data(what="eset")
eset <- de.ana(eset)

library("HTSanalyzeR")
library("snow")
library("GO.db")
library("KEGG.db")
geneList=ttmerge$logFC_WT_Klotho
names(geneList)=ttmerge$ENSEMBLEID
hits <- geneList[which(ttmerge$adj.P.Val_WT_Klotho < 0.05)]
geneList = annotationConvertor(geneList, species="Mm", initialIDs="Ensembl.gene",finalIDs="Entrez.gene", keepMultipleMappings=TRUE, verbose=TRUE)
hits = annotationConvertor(hits, species="Mm", initialIDs="Ensembl.gene",finalIDs="Entrez.gene", keepMultipleMappings=TRUE, verbose=TRUE)

GO_BP <- GOGeneSets(species="Mm", ontologies=c("BP"))
PW_KEGG <- kegg.gs
ListGSC <- list(GO_BP=GO_BP, PW_KEGG=PW_KEGG)
gsca <- new("GSCA", listOfGeneSetCollections=ListGSC,geneList=geneList, hits=names(hits))
options(cluster=makeCluster(8, "SOCK"))
gsca<-analyze(gsca, para=list(pValueCutoff=0.05, pAdjustMethod ="BH", nPermutations=1000, minGeneSetSize=5, exponent=1))
if(is(getOption("cluster"), "cluster")) {
  stopCluster(getOption("cluster"))
  options(cluster=NULL) 
}
summarize(gsca)
getTopGeneSets(gsca,"GSEA.results", c("GO_BP","PW_KEGG"), allSig=TRUE)


viewEnrichMap(gsca, resultName="HyperGeo.results",
              gscs=c("GO_BP"),  allSig=FALSE, ntop=30, gsNameType="term",
              displayEdgeLabel=FALSE, layout="layout.fruchterman.reingold")

gsca<-appendGSTerms(gsca, goGSCs=c("GO_BP"), keggGSCs=c("PW_KEGG"))
report(gsca,allSig=TRUE,reportDir="report_gsca",goGSCs=c("GO_BP"),keggGSCs="PW_KEGG",species="Mm")

data("KcViab_GSCA")
viewEnrichMap(KcViab_GSCA, resultName="HyperGeo.results",
              gscs=c("PW_KEGG"), allSig=FALSE, ntop=30, gsNameType="id",
              displayEdgeLabel=FALSE, layout="layout.fruchterman.reingold")


library(EnrichmentBrowser)
# clean the Mm_Migs set:
for (name in names(Mm.Migs) ) {
  names(Mm.Migs[[name]]) = sapply(names(Mm.Migs[[name]]), function (x) gsub("_","-",x) )
} 
# obtain grn for mouse
#pwys <- download.kegg.pathways("mmu")
#mmu.grn <- compile.grn.from.kegg(pwys)
# add Kegg and GO-Bp
Mm.Migs$KEGG <- get.kegg.genesets("mmu")
Mm.Migs$GO <- get.go.genesets(org="mmu", onto="BP", mode="GO.db")
# clean sets to only use annotated genes
all.genes=all.samps$dge$genes$entrezgene
Mm.Migs$GO = sapply( Mm.Migs$GO, function(x) x[x %in% all.genes] )
Mm.Migs$GO = Mm.Migs$GO[sapply( Mm.Migs$GO, function(x) length(x) > 5 )]
Mm.Migs$GO = Mm.Migs$GO[sapply( Mm.Migs$GO, function(x) length(x) < 400 )]
Mm.Migs$KEGG = sapply( Mm.Migs$KEGG, function(x) x[x %in% all.genes] )
Mm.Migs$KEGG = Mm.Migs$KEGG[sapply( Mm.Migs$KEGG, function(x) length(x) > 5 )]
Mm.Migs$KEGG = Mm.Migs$KEGG[sapply( Mm.Migs$KEGG, function(x) length(x) < 300 )]

for (time in names(GT.times.dges[["dge.times"]][["DT"]]) ) {
  timepoint=paste0("timepoint_",time)
  for (loc in names(GT.times.dges[["dge.times"]]) ){
    print(paste("going through",loc,time))
    fit=fits[["times"]][[loc]][[time]]
    dge = GT.times.dges[["dge.times"]][[loc]][[time]]
    do_gs_ebrowser(fit,dge,time,loc,perm=10000)
  }
}

