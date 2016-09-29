#setwd("/Volumes/Temp/Lukas/Data/Test/")
setwd("/Volumes/vetgrid01/Data/Test/")
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
library("WriteXLS")

#library(mixOmics)

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
    bname=paste("Enrichments/enr",loc,gt.name,cont,sep="_")
    for( db in names(Mm.Migs)){
      print(paste("going through",gt,loc,db,sep=" "))
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
        ea.browse(sbea.gsenr, html.only = TRUE)
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


get_de_lists <- function (fit,outspec,dge.sub=FALSE) {
  conts=colnames(fit$contrasts)
  de.list=list()
  # create directory for saving files
  #dir.create(gt,showWarnings ='F') 
  # get toptable for each contrast, merge and write to file
  conts=colnames(fit$contrasts)
  ttmerge = get_ttmerge(fit)
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
  if(dge.sub) {
    ttmerge = cbind(ttmerge, get.median.rpkm(dge.sub, genes = as.character(ttmerge$ENSEMBLEID)) )
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
  names(t.values) = genes
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


get.median.rpkm <- function (dge.sub, genes = c()){
  # gets a dge for a certain model and only returns the median rpkms for the groups of samples used in the model from all.samps$dge
  # create list of groups containing sample names
  #dge.sub=GT.times.dges$dge.GTs$DT$WT
  if (! "rpkm" %in% names(all.samps$dge) ) {
    all.samps$dge$rpkm=as.data.frame(rpkm(all.samps$dge))
  }
  if (length(genes) == 0) {
    genes=rownames(dge.sub$genes)
  }
  smp_names = split(rownames(dge.sub$targets),dge.sub$targets$group)
  return( 
    as.data.frame( sapply(smp_names, 
                          function(x) apply(all.samps$dge$rpkm[genes,x],1,median, na.rm = TRUE ) ) ) )
}

load("tue26012016.RData")
#library(EnrichmentBrowser)
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
# getting the pathways only
# library(KEGGREST)
# kegg.pw = keggList("pathway","mmu")
# kegg.pw=gsub("path:","",names(kegg.pw))
# k.dis = keggList("disease")
# k.pws = sapply(names(kegg.gs), function (x) strsplit(x,"_")[[1]][1] %in%  kegg.pw )



for (time in names(GT.times.dges[["dge.times"]][["DT"]]) ) {
  timepoint=paste0("timepoint_",time)
  for (loc in names(GT.times.dges[["dge.times"]]) ){
    print(paste("going through",loc,time))
    fit=fits[["times"]][[loc]][[time]]
    dge = GT.times.dges[["dge.times"]][[loc]][[time]]
    do_gs_ebrowser(fit,dge,time,loc,perm=10000)
  }
}

