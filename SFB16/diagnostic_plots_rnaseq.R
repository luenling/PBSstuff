library("limma")
library("edgeR")
library("DESeq2")
library("ggplot2")
library("RColorBrewer")
library("reshape2")
library("ggrepel")
library("factoextra")
library("pheatmap")
library("corrplot")
library("GGally")
library("org.Mm.eg.db")
require(hexbin)
require(gridExtra)
library("gPCA")

#library("Mus.musculus")
# set your working directory
mywd = "~/path/to/working/directory/"     
mywd = "~/Data/SFB61/" 
# execute commands with ctrl+enter in Rstudio
#setwd("C:/Users/Lenovo/Desktop/SFB/")
setwd(mywd); mywd                         # the wd is now set. ';' allows to execute multiple commands in one row


# read the Annotation data to an object. 
# in my case, the file Annotation.csv is in my working directory
annotation = read.csv("Annotation.csv")

# create copy of annotation table for subsetting using only RNA experiments 
ad = subset(annotation,grepl("RNA",library))
gene_annos = read.table("RNA.exome.Probes.csv",header=TRUE, sep = "\t")
# read the rnaseq counts table (rows for genes, columns for samples) 
rnaseq_counts= read.csv("RNA.exome.Counts.csv")

# get sample names to use
samples=ad$sample_name[grep('RNA_MMAK00[2-4].*T8',ad$sample_name,perl=T)]
# only use sample names in rnaseq_counts
samples=intersect(samples,names(rnaseq_counts))

samples2=ad$sample_name[grepl('STAT3FLOX',ad$genotype,perl=T) & ad$library == "RNA" & ad$cell_type != "M"]
stat3dge=create_DGE_object(samples2,grouping = c("cell_type"))
keep=check_expression_threshold(stat3dge,min_cpm = 1.0,min_count = 10,n=2)
table(keep)
stat3dge=stat3dge[keep,,keep.lib.sizes=FALSE]
colnames(stat3dge) = gsub('RNA_CBMF|H_H_SmartRNA_','',colnames(stat3dge),perl=T)
rownames(stat3dge$samples) = colnames(stat3dge)
stat3dge=calcNormFactors(stat3dge,method="TMM")
plotBCV_groups(stat3dge)
stat3_disp=estimateDisp(stat3dge)
plotBCV(stat3_disp)
plot_overall_counts(stat3dge)
plot_MA_matrix(stat3dge)
pca_stat3=nicePCA(stat3dge,shape_by ="cell_number",topvar = 1500)
#mds_stat3=niceMDS(stat3dge,shape_by ="cell_number")
pca_stat3_rlog=nicePCA(rlog(stat3dge$counts),annos = stat3dge$samples,shape_by ="cell_number",topvar = 1500)
# only look at cell number == 400
stat3_400=stat3dge$samples$cell_number == "400"
pca_stat3_400=nicePCA(stat3dge[,stat3_400],shape_by ="cell_number",topvar = 1000)
pca_stat3_100=nicePCA(stat3dge[,! stat3_400],shape_by ="cell_number",topvar = 1000)

get_r2_corPC(pca_stat3,stat3dge$samples[,c("group","cell_number")])
get_r2_corPC(pca_stat3_rlog,stat3dge$samples[,c("group","cell_number")])
pca_stat3=nicePCA(stat3dge,color_by ="lib.size",shape_by ="cell_number",topvar = 0)

perform_gPCA(stat3dge,"date_of_collection")

# create an EdgeR DGE object to work with
mmdge=create_DGE_object(samples,grouping = c("genotype"))
# plot the fraction of genes filtered out
plot_filtered_genes(mmdge,per_lib = T)
# plot the fraction of genes with zero counts
plot_filtered_genes(mmdge,per_lib = F, n=3)


# remove all genes underneath a certain count per million and minimal count threshold
keep=check_expression_threshold(mmdge,min_cpm=1.0,min_count=10)

# look at how many genes are to be kept
table(keep)
mmdge=mmdge[keep,  , keep.lib.sizes = FALSE]
colnames(mmdge) = gsub('RNA_MMA|H_H_SmartRNA_','',colnames(mmdge),perl=T)
rownames(mmdge$samples) = colnames(mmdge)


# the readcounts need to be normalised (other than only by library size) to be comparable between samples
# for different methods see Maza, E. et al. (2013). Comparison of normalization methods for differential gene expression analysis in RNA-Seq experiments: A matter of relative size of studied transcriptomes. Communicative & Integrative Biology, 6(6), e25849.
# normalise data using TMM amd look at density plot (normally similar to raw CPM plot, the densities should be quite similar between samples) 
mmdge <- calcNormFactors(mmdge, method = "TMM")


# look at biological coefficient of variation (BCV)
# explained in https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf chapt.:2.8.2
# the BCV should correspond to the variation between replicates that is not just due to poisson sampling - that is additional biological and technical variation
# a value of eg. 0.25 means something like that the estimated true abundance of a gene varies by 25% between replicates
# from Chen et al., F1000Research 2016, 5:1438: http://dx.doi.org/10.12688/f1000research.8987.2
# For RNA-seq studies, the NB dispersions tend to be higher for genes with very low counts. The dispersion trend tends to decrease smoothly with abundance and to asymptotic to a constant value for genes with larger counts. From our past  experience, the asymptotic value for the BCV tends to be in range from 0.05 to 0.2 for genetically identical mice or cell lines, whereas somewhat larger values (> 0.3) are observed for human subjects.
mmdge_disp=estimateDisp(mmdge)
plotBCV(mmdge_disp)
plotBCV_groups(mmdge)
plot_overall_counts(mmdge)

# for comparing technical replicates create a factor called techrep that distinguishes groups of technical replicates
# in some cases just group, experiment_id and date_of_collection
mmdge$samples$techrep = factor(paste(mmdge$samples$date_of_collection,mmdge$samples$group,mmdge$samples$experiment_id))

# MA plots between all replicates of a group
# each point corresponds to a gene with the log expression averaged over a group of samples (A) on the x axis and the log2 fold expression change of a sample against the average as the y axis (M)
# this plot can reveal systematic biases (especially comparing a sample to other replicates) or changes in expression of large number of genes (against the overall average)

#plot_MA_matrix(mmdge,averageover=c("all","group"))
plot_MA_matrix(mmdge,averageover=c("group"))
plot_MA_matrix(mmdge,averageover = c("techrep"))

# plot Multidimensional scaling (MDS) plots for samples using plotMDS from the limma package
# Like Principle component analysis, MDS allows to represent high dimensional data in lower dimensions
# from Chen et al., F1000Research 2016, 5:1438, http://dx.doi.org/10.12688/f1000research.8987.2 :
# In the MDS plot, the distance between each pair of samples can be interpreted as the leading log-fold change between the samples for the genes that best distinguish that pair of samples. By default, leading fold-change is defined as the root-mean-square of the largest 500 log2-fold changes between that pair of samples.
# in this plots, replicates should be closer to each other, than unrelated samples
# you can color/label samples by different annotations, to see whether other influences lead to clustering of samples 
# the MDS dimensions/components are ranked by the amount of difference the explain between samples. you can also look at higher change dimensions, apart from the first 2
# the direct plots produced are alright, but do not look as nice and are a bit harder to configure, so we just extract the result data and plot it in ggplot
mds_gg=niceMDS(exobj = mmdge,color_by="genotype",shape_by = "replicates")


# we can use DESeq2 to obtain a variance stabilized approximately log2 expression matrix, in which variances for lowly expressed genes are shrunk more than those of highly expressed ones
# the variance stabilizing transformation is faster than the rlog, but more sentsitive to big differences in size factors
# the vst function is even faster, through only subsampling gense to estimate dispersion trends, but normally not blidned to the experimental design (see ?vst)
# create DESeq2 Data set and get either the vst or the rlog matrices
mm_rlog=rlog(mmdge$counts)
#mm_vst=vst(mmdge$counts)
mm_vst=varianceStabilizingTransformation(mmdge$counts)
q1=qplot(rowMeans(cpm(mmdge,log=T)),rowVars(mm_rlog)) + stat_binhex() + #geom_smooth(color="red",se=FALSE) + 
  stat_summary_bin(fun.y = median,color="green",geom="line") + ggtitle("rlog")
q2=qplot(rowMeans(cpm(mmdge,log=T)),rowVars(mm_vst)) + stat_binhex() + #geom_smooth(color="red",se=FALSE) + 
  stat_summary_bin(fun.y = median,color="green",geom="line") + ggtitle("varianceStabilizingTransformation")
q3=qplot(rowMeans(cpm(mmdge,log=T)),rowVars(cpm(mmdge,log=T))) + stat_binhex() + #geom_smooth(color="red",se=FALSE) + 
  stat_summary_bin(fun.y = median,color="green",geom="line") + ggtitle("logCPM, TMM")
grid.arrange(q1,q2,q3,nrow=3)

pca_mm=nicePCA(mmdge,shape_by ="date_of_collection",topvar = 1500)
get_r2_corPC(pca_mm,annos = mmdge$samples[,c("group","date_of_collection","cell_number")])


pca_mm_rlog=nicePCA(exobj = mm_rlog,annos = mmdge$samples,shape_by ="date_of_collection",topvar = 1500)


get_r2_corPC(pca_mm,mmdge$samples[,c("group","cell_number","date_of_collection")])
get_r2_corPC(pca_mm_rlog,mmdge$samples[,c("group","cell_number","date_of_collection")])

plot_MA_matrix <- function(mmdge, averageover=c("all","group","techrep")){
  # calculates average A and M values per all samples, group and techrep (if defined)
  # MA plots between all replicates of a group
  # each point corresponds to a gene with the log expression averaged over a group of samples (A) on the x axis and the log2 fold expression change of a sample against the average as the y axis (M)
  # this plot can reveal systematic biases (especially comparing a sample to other replicates) or changes in expression of large number of genes (against the overall average)
  # need to first create matrix of A (average log expression) and M (log fold expression difference) values
  # all is plotted in one big matrix to plot, use ggplot with hexagonal binning to reduce the number of points
  # define techrep factor for group membership of technical replicates
  AMm=get_A_M_values(mmdge,averageover)
  # plot MA of each sample against overall mean of all other samples
  # the red trend line indicates the mean M value and should be close to y=0
  if ("all" %in% averageover) {
    p = ggplot(AMm,aes(x=A,y=M)) + stat_binhex() + geom_smooth(color="red",se=FALSE) + 
      facet_grid(rows=vars(group),cols=vars(replicates)) + 
      labs(y = "log-ratio (sample against average)", x = "average log CPM (all samples)" )
    plot(p)
  }
  # plot MA of each replicate against mean of its sample group, these should be more similar and can reveal problems between replicates
  if ("group" %in% averageover) {
    p =  ggplot(AMm,aes(x=Ag,y=Mg)) + stat_binhex() + geom_smooth(color="red",se=FALSE) + 
      facet_grid(rows=vars(group),cols=vars(replicates)) + 
      labs(y = "log-ratio (sample against group average)", x = "average log CPM (replicate samples in group)" )  
    plot(p)
  }
  # plot MA of each technical replicate against mean of all technical replicates , these should be very similar 
  if ("techrep" %in% averageover & "techrep" %in% colnames(mmdge$samples) ) {
    p =   ggplot(AMm,aes(x=Ami,y=Mmi)) + stat_binhex() + geom_smooth(color="red",se=FALSE) + 
      facet_grid(rows=vars(techrep),cols=vars(replicates)) + 
      labs(y = "log-ratio (sample against group average)", x = "average log CPM (replicate samples in techrep)" ) 
    plot(p)
  }
}

plot_filtered_genes <- function(mmdge,minval=c(0.5,1,2.5,5,7.5,10,15),n=3,use.cpm=TRUE,palette="Blues",per_lib=FALSE) {
  # n: number of samples that have to have expression above threshold
  # per_lib: only look at genes filtered per library
  
  grp_samp_cols=get_group_sample_cols(mmdge)
  counts=if(use.cpm) cpm(mmdge) else mmdge$counts
  if (per_lib){
    fracs= sapply(minval,function (x) colSums(counts < x) )/nrow(counts)
  } else {
    fracs = sapply(minval,function (x) { keep=rowSums(counts >= x) > n;  colSums(counts[keep,] == 0)/sum(keep) }  )
  }
  colnames(fracs)=minval
  old.par=par()
  par(mar=c(10,3,4,5))
  leg_title=if(use.cpm) "min. cpm" else "min. count"
  main_title=if(per_lib) "Fraction of genes filtered with threshold per library" else paste("Fraction of unexpressed genes after filtering above threshold in at least",as.character(n),"samples") 
  barplot(t(fracs),beside=TRUE,legend.text=TRUE,args.legend=list(title=leg_title,x ="right",bty="n",inset=c(-0.10,0),xpd = TRUE),
          col=brewer.pal(n=ncol(fracs),name=palette),cex.names=0.5,las=2,main=main_title)
  par(old.par)
}


plotBCV_groups <- function(dgelist){
  # plot the biological coefficient of variation per group
  groups <- levels(dgelist$samples$group)
  if (length(groups) == 3) {
    par(mfrow = c(2,2))
  } else {
    par(mfrow = n2mfrow(length(groups)))
  }
  for(i in groups) {
    x <- estimateDisp(dgelist[,dgelist$samples$group == i])
    plotBCV(x,main=i)
  }
  # fill up page so not to leave plots empty
  while ( ! par()$page) { plot.new() }
  par(mfrow=c(1,1))
}


plot_overall_counts <- function(mmdge,normalized=TRUE,prior.count=0.25) {
  ## prior.count : added for getting logcpm not to be bothered by 0 vals   
  grp_samp_cols=get_group_sample_cols(mmdge)
  pcpm=prior.count/10^6
  norm=if(normalized) "TMM normalized" else ("library size normalized")
  # boxplot of log(CPM) with CPM: counts per million reads in library
  par(mfrow=c(2,2))
  names=paste(mmdge$samples$group,mmdge$samples$replicates,sep="_")
  boxplot(log10(cpm(mmdge,normalized.lib.sizes = normalized)+pcpm)[1:100,],col=c(grp_samp_cols$sample_cols),las=2,cex.axis=0.75, ylab="log10(CPM)",names=names,main=paste(norm,"CPM values of genes"),las=2)
  
  # plot library sizes (reads mapped to filtered genes)
  barplot(mmdge$samples$lib.size/10^6,col=grp_samp_cols$sample_cols,main="Library sizes",names.arg=names, las=2, cex.names =0.75, ylab="million reads")
  ## plot fraction of zero counts per sample
  barplot(apply(mmdge$counts,2,function (x) length(x[ x == 0])/length(x)),col=grp_samp_cols$sample_cols,
          names.arg=names, las=2, cex.names =0.75, ylab="fraction 0 count genes",main="genes with 0 counts")
  plot.new()
  par(xpd=TRUE)
  legend("topleft",legend=paste(names,rownames(mmdge$samples),sep=" : "),cex=0.75,fill=grp_samp_cols$sample_cols)
  par(xpd=FALSE)
  par(mfrow=c(1,1))
  # plot density of log cpm
  # change data frame to be used with ggplot
  lcpm_m=melt(log10(cpm(mmdge,normalized.lib.sizes = normalized)+pcpm),value.name = "log10CPM")
  colnames(lcpm_m)[1:2]=c("gene","sample")
  # density plot using ggplot2
  p = ggplot(lcpm_m, aes(x=log10CPM,color=sample)) + stat_density(geom="line",position = "identity") + 
    scale_color_manual(values = grp_samp_cols$sample_cols) + guides(col=guide_legend(ncol=1)) + labs(x="log10(CPM)", title=paste(norm,"CPM values of genes") ) +
    theme(legend.title = element_text(), legend.text = element_text(size=7.0))
  plot(p)
}


create_DGE_object <- function(samples, anno_table=ad, count_table=rnaseq_counts, gene_annos=NULL, grouping=c("genotype")){
  # create a DGE List object with only the samples in the vector samples
  # you need to say which columns of the annoation sheet are used for creating groups to be compared
  # eg. grouping=c("genotype","treatment")
  # create a sample list to subset with ( in my case not the automaticly generated sample name but the submitted sample_name was used for the sample columns)
  # you can also give the gene annotation table to obtain gene lengths and Entrez IDs
  anno_table=subset(anno_table, sample_name %in% colnames(count_table))
  # only use sample names in rnaseq_counts
  samples=intersect(samples,names(count_table))
  mmdata=subset(anno_table, sample_name %in% samples)
  # remove unused factor levels
  mmdata=droplevels(mmdata)
  #sample_list=mmdata$sample_name  
  rnaseq_counts=count_table[,colnames(count_table) %in% samples]
  # just to make sure, mmdata and the count matrix have the same order
  mmdata = mmdata[which(mmdata$sample_name %in%  colnames(rnaseq_counts)),]
  # define the grouping variable - in this case just treatment, otherwise combine multiple columns like so:
  # mmdata$grouping = factor(with(mmdata, paste(treatment,genotype)))
  rownames(mmdata)=mmdata$sample_name
  mmdata$grouping = if (length(grouping)>1) factor(do.call(paste,mmdata[,grouping])) else mmdata[,grouping]
  # if you want to have a certain order of grouping levles for analysis, eg. WT first
  mmdata$grouping = factor(mmdata$grouping,sort(levels(mmdata$grouping),decreasing=TRUE) ) 
  # assign replicate number just sequenctially from 1 to n for each samlpe in each grouping factor
  mmdata$replicates=sapply(1:length(mmdata$grouping), function(x) sum(mmdata$grouping[1:x] == mmdata$grouping[x]))
  # sort by groupings
  mmdata=mmdata[order(as.numeric(mmdata$grouping)),]
  rnaseq_counts=rnaseq_counts[,rownames(mmdata)]
  mmdge=DGEList(rnaseq_counts, samples = mmdata, group =  mmdata$grouping)
  # add information on genes to dge object
  # if there are multiple annotations for an ENSEMBL ID we just take the first, for non unique Entrez IDs we also only keep the first (could do better by checking whether chromosome or gene symbol NA)
  mmdge$genes=data.frame("ENSEMBL"=rownames(mmdge$counts))
 
    #mmdge$genes=select(Mus.musculus,keys = rownames(mmdge),keytype = "ENSEMBL",columns=c("ENSEMBL","SYMBOL","GENEID","GENENAME","CDSCHROM"))
    mmdge$genes$ENTREZID=mapIds(org.Mm.eg.db,keys = rownames(mmdge),keytype = "ENSEMBL",column="ENTREZID",
                                multiVals = "first")
    mmdge$genes$SYMBOL=mapIds(org.Mm.eg.db,keys = rownames(mmdge),keytype = "ENSEMBL",column="SYMBOL",
                              multiVals = "first")
    mmdge$genes$GENENAME=mapIds(org.Mm.eg.db,keys = rownames(mmdge),keytype = "ENSEMBL",column="GENENAME",
                                multiVals = "first")
    if (is.null(gene_annos)) {
      mmdge$genes$Length = gene_annos$width
      mmdge$genes$SYMBOL=gene_annos$gene
    }
  
  #mmdge$genes$CDSCHROM=mapIds(Mus.musculus,keys = rownames(mmdge),keytype = "ENSEMBL",column="CDSCHROM",multiVals = "first")
  mmdge$genes$UniqueEntrez=mmdge$genes$ENTREZID
  mmdge$genes$UniqueEntrez[duplicated(mmdge$genes$UniqueEntrez)] = NA
  rownames(mmdge$genes)=rownames(mmdge$counts)
  return(mmdge)
}
# calculate how many replicates in each group are above 1 cpm expression
# you might have to alter that threshold depending on sequencing depth and other considerations, for a library size of 12 million reads, that means around 12 counts
# see suggestions in Filtering to remove low counts in Chen et al., F1000Research 2016, 5:1438, http://dx.doi.org/10.12688/f1000research.8987.2

check_expression_threshold <- function(exobj,n=NULL,min_cpm=1.0,min_count=NULL,by_group=TRUE,groups=NULL){
  # this is similar to filterByExpr in the edgeR package
  # check ?filterByExpr and maybe use that if you feel more comfortable
  # returns a logical vector with length of genes with TRUE if expression in n replicates of a group is higher than min_cpm and min_count
  # n defaults to int(3/4) of the smalles group size. if by_group=FALSE no grouping is used, groups can be a factor vector indicating group membership, if a DGEList object is used and group is NULL, groups=dgeobject$samples$group
  # if min_count or min_cpm are set to NULL, either is ignored
  if (by_group){
    if (is.null(groups)) {
      if ("samples" %in% names(exobj) && ! is.null(exobj$samples$group)) {
        groups <- exobj$samples$group
      }
      else {
        groups=as.factor(rep(1,ncol(exobj)))
      }
    }
    else {
      groups = as.factor(groups)
    }
  } 
  else {
    groups=as.factor(rep(1,ncol(exobj)))
  }
  if (is.null(n)){
    n <- as.integer(min(tabulate(groups))*3/4)
  }
  keep=rep(TRUE,nrow(exobj))
  if (! is.null(min_cpm)) {
    suf_by_group=tapply(colnames(exobj),groups, function (x) rowSums(cpm(exobj)[,x] >= min_cpm ))
    suf_by_group=sapply(names(suf_by_group), function (x) suf_by_group[[x]])
    keep=keep & (rowSums(suf_by_group >= n) > 0)
  }
  if (! is.null(min_count)) {
    suf_by_group=tapply(colnames(exobj),groups, function (x) rowSums(as.matrix(exobj[,x]) >= min_count ))
    suf_by_group=sapply(names(suf_by_group), function (x) suf_by_group[[x]])
    keep=keep & (rowSums(suf_by_group >= n) > 0)
  }
  return(keep)
}

#set up colors for samples (similar colors for replicates, can be used later for plots)
get_group_sample_cols <- function(groups=mmdge$samples$group,replicates=mmdge$samples$replicates,palette="Dark2"){
  #takes two factor vectors for group and replicates and returns a vector for group and sample colors according to group and replicate vectors as a list
  # if the first element is a DGElist object, it extracts groups and replicates from that
  if(class(groups) == "DGEList") {
    replicates = groups$samples$replicates
    groups = groups$samples$group
  }
  group_cols=brewer.pal(n=length(levels(groups)),name=palette)
  names(group_cols)=levels(groups)
  sample_cols=group_cols[groups]
  # change the colors a tiny bit for each replicate
  sample_cols=sapply(1:length(sample_cols), function (x) 
    adjustcolor(sample_cols[x],transform=diag(c(1, 1, 1, 1) - (replicates[x]-mean(replicates))*c(0.20,0.20,0.20,0) ) ))
  names(sample_cols)=rownames(replicates)
  return(list("group_cols"=group_cols,"sample_cols"=sample_cols))
}


get_A_M_values <- function(mmdge,avagainst=c("all","group")){
  # this function returns M and A values averaged over all samples, by group, and by techrep
  # gives table with columns:
  # A, M : averaged over all samples
  # Ag,Mg : averaged by group
  # Ami,Mmi : averaged over technical replicates
  # averages over all samples
  lcpm_n=cpm(mmdge,log=TRUE)
  As=rowMeans(lcpm_n) 
  # Ms over all sample averages
  # melt data frame for use in ggplot2
  Mm=melt(lcpm_n)
  colnames(Mm)=c("gene","sample","lCPM")
  AMm=cbind(Mm,As[Mm$gene])
  colnames(AMm)=c("gene","sample","lCPM","A")
  AMm$M=AMm$lCPM-AMm$A
  if ("techrep" %in% avagainst & "techrep" %in% colnames(mmdge$samples)){
    AMm=cbind(AMm,mmdge$samples[Mm$sample,c("group","replicates","techrep")])
  } else {
    AMm=cbind(AMm,mmdge$samples[Mm$sample,c("group","replicates")])
  }
  # calculate A values per group (not overall)
  # this might be more appropriate, to see whether replicates deviate from each other
  if ("group" %in% avagainst){
    Ag=tapply(AMm$lCPM,AMm[,c("gene","group")],mean,simplify = TRUE)
    AMm$Ag=sapply(seq(1,nrow(AMm)), function (x) Ag[AMm$gene[x],AMm$group[x]] ) 
    AMm$Mg=AMm$lCPM-AMm$Ag
  }
  # this might be more appropriate, to see whether technical replicates deviate from each other for each indvidual
  if ("techrep" %in% avagainst & "techrep" %in% colnames(AMm)) {
    Ami=tapply(AMm$lCPM,AMm[,c("gene","techrep")],mean,simplify = TRUE)
    AMm$Ami=sapply(seq(1,nrow(AMm)), function (x) Ami[AMm$gene[x],AMm$techrep[x]] ) 
    AMm$Mmi=AMm$lCPM-AMm$Ami
  }
  return(AMm)
}





niceMDS <- function(exobj,annos=NULL,color_by="group", shape_by="replicates",
                    lables=NULL,ndims=5,xdim="Dim 1",ydim="Dim 2", title=NULL,top=500,gene.selection="pairwise"){
  ### takes an object which plotMDS from the limma package can take and an annotations sheet and plots the dimensions given by xdim and ydim ("Dim n")
  ### for coloring, shapes and labelling give the name of the factor column in the annotation sheet
  ### returns a data frame with the MDS dimensions and the annotation sheet concatenated for easy plotting with ggplot2
  # plot Multidimensional scaling (MDS) plots for samples using plotMDS from the limma package
  # Like Principle component analysis, MDS allows to represent high dimensional data in lower dimensions
  # from Chen et al., F1000Research 2016, 5:1438, http://dx.doi.org/10.12688/f1000research.8987.2 :
  # In the MDS plot, the distance between each pair of samples can be interpreted as the leading log-fold change between the samples for the genes that best distinguish that pair of samples. By default, leading fold-change is defined as the root-mean-square of the largest 500 log2-fold changes between that pair of samples.
  # in this plots, replicates should be closer to each other, than unrelated samples
  # you can color/label samples by different annotations, to see whether other influences lead to clustering of samples 
  # the MDS dimensions/components are ranked by the amount of difference the explain between samples. you can also look at higher change dimensions, apart from the first 2
  # the direct plots produced are alright, but do not look as nice and are a bit harder to configure, so we just extract the result data and plot it in ggplot
  mdsdata=plotMDS(exobj, top=top, gene.selection=gene.selection, ndim=ndims,plot=FALSE)
  if (is.null(annos) & is(exobj,"DGEList")) {
      annos = exobj$samples
  }
  mds_gg=cbind(mdsdata$cmdscale.out,annos[rownames(mdsdata$cmdscale.out),])
  mds_gg[c(color_by,shape_by)] <- lapply(mds_gg[c(color_by,shape_by)] , factor)
  #mds_gg$replicates=as.factor(mds_gg$replicates)
  dim_names = paste("Dim",1:ndims)
  colnames(mds_gg)[1:ndims]=dim_names
  #xdim=dim_names[1]; ydim=dim_names[2];
  gg = ggplot(mds_gg,aes(x = !!ensym(xdim), y = !!ensym(ydim),colour=!!ensym(color_by),shape=!!ensym(shape_by))) + geom_point() +
    labs(x = paste("leading logFC",xdim), y = paste("leading logFC",ydim) ) + ggtitle(title) +
    geom_hline(color="grey", linetype="dashed", yintercept = 0) + geom_vline(color="grey", linetype="dashed", xintercept=0)
  if (!is.null(lables)) {
    gg = gg + geom_label_repel(aes(label = mds_gg[,lables]))
  }
  plot(gg)
  return(mds_gg)
}

## helper function to only get the topvar rows of a data frame with highest variation
get_top_var<- function(exprs_mat,topvar=1000){
  ## only returns the topvar rows with the highest variation from a dataframe/matrix exprs_mat
  row_vars <- apply(exprs_mat,1,var)
  if ( ! is.finite(topvar) || topvar < 1 || topvar > length(row_vars)) {
    topvar = length(row_vars)
  }
  topvar<-as.integer(topvar)
  return(exprs_mat[order(row_vars,decreasing = TRUE)[1:topvar],])
}

get_topvar_by_group <- function(exprs_mat,groups,topvar=100) {
  # performs variance estimation per level in factor vector group for all rows in exprs_mat
  # returns the topvar genes for each group with the highest within group variance
  row_vars <- as.data.frame(sapply(levels(groups), function(x) rowVars(exprs_mat[],cols=which(groups == x))))
  rownames(row_vars) <- rownames(exprs_mat)
  if ( ! is.finite(topvar) || topvar < 1 || topvar > nrow(row_vars)) {
    topvar = nrow(row_vars)
  }
  topvar<-as.integer(topvar)
  top_row_var <- lapply(colnames(row_vars),function(x) { y<-row_vars[order(row_vars[,x],decreasing = TRUE)[1:topvar],x,drop=FALSE]; return(setNames(y[[1]],rownames(y))) } )
  names(top_row_var) <- colnames(row_vars)
  top_row_var[["COMBINED"]]=unique(as.vector(sapply(top_row_var,names)))
  names(top_row_var[["COMBINED"]])=top_row_var[["COMBINED"]]
  return(top_row_var)
  
}


# Performing Principle component analysis (PCA) and plotting Principle Components for the n genes varying most over all samples
# Here we only take the 1000 most varying genes for doing a PCA, similar to what DESeq2 does with their plotPCA function
nicePCA <- function(exobj=cpm(mmdge,log=TRUE),annos=mmdge$samples,topvar=1000,color_by="group", shape_by="replicates",lables=NULL,xdim="PC1",ydim="PC2",center=TRUE,scale=FALSE,title=NULL,plot_skree=TRUE){
  ### takes a numerical matrix with expression values of genes x sample (eg. cpm(mmdge,log=TRUE)) and an annotations sheet, performes centering over genes and PCA on the transposed matrix and plots a skree plot and the dimensions given by xdim and ydim ("Dim n")
  ### the PCA is performed on the topvar most varying genes (default 5000), if topvar = 0,NA,NULL or similar all genes are considered
  ### for coloring, shapes and labelling give the name of the factor column in the annotation sheet
  ### returns a prcomp output with a data frame added under pca$tab with the principle components and the annotation sheet concatenated for easy plotting with ggplot2
  ### if given a prcomp object, only does the plotting
  if (! class(exobj) == "prcomp") {
    if (is(exobj,"DGEList")) {
       annos=exobj$samples;
       exobj=cpm(exobj,log=TRUE);
    }
    exobj_top_var=get_top_var(exobj,topvar)
    # perform PCA (transposed exobj matrix to do it in sample space and to center over each gene (columns), not samples)
    pca <- prcomp(t(exobj_top_var), center = center, scale. = scale) 
    #skree plot (Variance explained by the different components)
    if (plot_skree) {
      barplot(pca$sdev^2/sum(pca$sdev^2),ylab="% variance explained",names.arg = colnames(pca$x),las=2)
    }
  } else {
    pca <- exobj
  }
  if (is.data.frame(annos)){
    pca$tab=cbind(pca$x,annos[rownames(pca$x),])
    pca$tab[c(color_by,shape_by)] <- lapply(pca$tab[c(color_by,shape_by)] , factor)
  }
  # get nice strings for variance explained for PCs
  var_expl=paste0("(",round(100*(pca$sdev^2/sum(pca$sdev^2)),digits=1),"%)")
  names(var_expl)=colnames(pca$x)
  # do actual plotting
  gg = ggplot(pca$tab,aes(x = !!ensym(xdim), y = !!ensym(ydim),colour=!!ensym(color_by),shape=!!ensym(shape_by))) + 
    geom_point(size=3) + labs(x = paste(xdim,var_expl[xdim]), y = paste(ydim,var_expl[ydim]) ) + ggtitle(title) +
    geom_hline(color="grey", linetype="dashed", yintercept = 0) + geom_vline(color="grey", linetype="dashed", xintercept=0)
  if (!is.null(lables)) {
    gg = gg + geom_label_repel(aes(label = pca$tab[,lables]))
  }
  plot(gg)
  return(pca)
}


get_r2_corPC <- function(pca,annos,dims=c(1:4)){
  res=list()
  for (x in names(annos)){
    res[[x]]=data.frame(p_lm=rep(NA,length(dims)), adj.r2=rep(NA,length(dims)),p_krusk=rep(NA,length(dims)))
    rownames(res[[x]])=as.character(dims)
    for(i in dims){
      res[[x]][as.character(i),"adj.r2"] <- summary(lm(pca$x[,i] ~ annos[,x]))$adj.r.squared
      #f <- a$fstatistic
      #res[[x]][as.character(i),"p_lm"] <- pf(f[1],f[2],f[3],lower.tail=F)
      res[[x]][as.character(i),"p_lm"] <- anova(lm(pca$x[,i] ~ annos[,x]))$`Pr(>F)`[1]
      res[[x]][as.character(i),"p_krusk"] <-kruskal.test(pca$x[,i]~annos[,x])$p.value
    }
  }
  return(res)
}










get_var
# we can use DESeq2 to obtain a variance stabilized approximately log2 expression matrix, in which variances for lowly expressed genes are shrunk more
# the variance stabilizing transformation is faster than the rlog, but more sentsitive to big differences in size factors
# the vst function is even faster, through only subsampling gense to estimate dispersion trends, but normally not blidned to the experimental design (see ?vst)
# create DESeq2 Data set and get either the vst or the rlog matrices
mmdeseq = DESeqDataSetFromMatrix(countData = mmdge$counts,
                                 colData = mmdge$samples,
                                 design = ~ group)
mmdeseq <- estimateSizeFactors(mmdeseq)
#mmdeseq <- estimateDispersions(mmdeseq)
#plotDispEsts(mmdeseq)
mmdeseq_vst <- varianceStabilizingTransformation(mmdeseq, blind = TRUE)
#mmdeseq_rlog <- rlog(mmdeseq,blind=TRUE)
#pca=nicePCA(assay(mmdeseq_rlog))
rlog_dist=dist(t(assay(mmdeseq)))

# you can access the log2 expression data of a deseq objet using eg.: assay(mmdeseq_vst)

pca=nicePCA(exobj=assay(mmdeseq_vst),topvar=1000,color="group",shape="replicate")
pca=nicePCA(exobj=cpm(mmdge,log=TRUE),topvar=1000,color="group",shape="replicates")
nicePCA(pca,color="genotype",shape="date_of_collection",lables="mouse_id")
mmdge_nobatch=removeBatchEffect(cpm(mmdge,log=TRUE),batch=mmdge$samples$cell_number,design=model.matrix(~1+group,data=mmdge$samples))
pca_nobatch=nicePCA(mmdge_nobatch,shape="cell_number",lables = "date_of_collection")

# plot first n=5 PC against each other 

#get sample symbols and colors
pchars=mmdge$samples[row.names(pca$x),c("replicates")]
#pcols=grp_samp_cols$group_cols[mmdge$samples[row.names(pca$x),c("group")]]

# for legend creation
reps=as.character(unique(pchars))
groups=unique(names(grp_samp_cols$group_cols))
# padding for legend
pad_r=max(0,length(groups)-length(reps))
pad_g=max(0,length(reps)-length(groups))
# pairs plot for the first 5 PCs
# add lines for zero axis
zero_lines <- function(x,y,...){
  points(x,y,...)
  abline(h=0,col = "lightgray", lty = 3)
  abline(v=0, col = "lightgray", lty = 3)
}
# plot pairwise comparisons
pairs(pca$x[,1:5],col=pcols,pch=pchars,lower.panel = NULL,upper.panel = zero_lines)
# add legend
par(xpd=NA)
legend("bottomleft", c(groups,rep.int(" ",pad_g),paste("rep",reps),rep.int(" ",pad_r)),ncol=2,
       col=c(grp_samp_cols$group_cols[groups],rep("black",pad_g+pad_r+length(reps))) , 
       pch=c(rep.int(0,length(groups)),rep.int(NA,pad_g),unique(pchars), rep.int(NA,pad_r))  )
par(xpd=FALSE)


# plot hierarchical clustering and heatmap with distances
# the clustering is done on the Pearson gene-gene, r, correlation of samples against each other
# 
# as we only want to cluster by distance, we need to create a distance matrix
# We use sqrt(2*(1-r)) to obtain a distance for clustering from the correlation r
cor_mat=cor(lcpm_n, method="pearson")
#lcmp_cor_dist=as.dist(sqrt(2*(1-cor_mat)))
lcmp_cor_dist=as.dist(1-cor_mat)
plot(hclust(lcmp_cor_dist))

pca=nicePCA(mmdge)


# Perform a DE analysis for genotypes without batch correction
design_mat=model.matrix(~ mmdge$samples$group)
#mmdge=voomWithQualityWeights(mmdge,design = design_mat, plot=TRUE)
v1=voom(mmdge,design = design_mat, plot=TRUE)
fit <- lmFit(v1, design = design_mat)
fit <- eBayes(fit)
summary(decideTests(fit))

design_mat2=model.matrix(~ mmdge$samples$group + mmdge$samples$date_of_collection)
v2=voom(mmdge,design = design_mat2, plot=TRUE)
fit2 <- lmFit(v2, design = design_mat2)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))

fit2 <- lmFit(v2, design = design_mat2) comobj

dupcor <- duplicateCorrelation(v2$E, design=design_mat2, block=mmdge$samples$mouse_id) 
dupcor$consensus.correlation
fit3 <- lmFit(v2, design=design_mat2, block=mmdge$samples$mouse_id, correlation=dupcor$consensus.correlation)
# the block indicates the technical replicates.
fit3 <- eBayes(fit3)
summary(decideTests(fit3))


# remove batch effect for data of collection to see if PCA gets any better
lcpm_batch=removeBatchEffect(lcpm_n,batch=mmdge$samples$date_of_collection,design=design_mat)
# calculate the consensus correlation for technical replicates
corfit <- duplicateCorrelation(lcpm_n,design=NULL,block=mmdge$samples$mouse_id)
corfit$consensus.correlation
#lcpm_batch=removeBatchEffect(lcpm,batch=NULL,design=design_mat) 
lcmp_batch_cor_dist=1-cor(lcpm_batch, method="pearson")
plot(hclust(as.dist(lcmp_batch_cor_dist)))


# to color dendrograms we have to use a different package
#install.packages("dendextend")
library(dendextend)
dd= as.dendrogram(hclust(as.dist(lcmp_cor_dist)))
# get the right colors per group
dd_lab_col=grp_samp_cols$group_cols[ mmdge$samples[names(labels_colors(dd)),"group"]]
names(dd_lab_col)=names(labels_colors(dd))
# set colors for dendrogram
dd=set(dd,"labels_color",dd_lab_col)
dd=set(dd, "labels_cex", 0.75)
# plot dendrogram
plot(dd)

# to use a different distance eg. use dist for euclidean:
lcmp_eucl_dist=dist(t(lcpm_n), method = "euclidean")
# plotting just a normal dendrogram
plot(hclust(lcmp_eucl_dist))



# define annotation table
annos=mmdge$samples[,c("group","replicate")]
annos$replicate=as.factor(annos$replicate)
pheatmap(as.dist(lcmp_cor_dist), annotation = annos)

# for more customisable but slower plotting
library("ComplexHeatmap")
ha=HeatmapAnnotation(annos)
Heatmap(lcmp_cor_dist,name="1-r",show_row_dend = FALSE, column_names_gp = gpar(fontsize=9), 
        row_names_gp = gpar(fontsize=9),top_annotation = ha)


# if you want, you use a refined method for getting normalise, moderated log2 values for clustering with the rlog function implemented in DESeq2
library("DESeq2")
# create DESeq2 Data set and get the rlog matrix
mmdeseq = DESeqDataSetFromMatrix(countData = mmdge$counts,
                                 colData = mmdge$samples,
                                 design = ~ group)
mmdeseq_rlog=rlog(mmdeseq,blind=TRUE)
nicePCA(assay(mmdeseq))
rlog_dist=dist(t(assay(mmdeseq)))

plot(hclust(rlog_dist))
pheatmap(rlog_dist, annotation = annos)




# let's get an overview about the possible values in each column
levels(ad$group_ID)
levels(ad$Experiment_ID)
levels(ad$genotype)
levels(ad$cell_type)

# subsetting of samples with factors
mmdata = ad[ad$Experiment_ID == "MM",]; length(mmdata$group_ID)
# drop unused factor levels
mmdata = lapply(mmdata, function(x) if(is.factor(x)) factor(x) else x)
mmdata

#subsetting with a list of sample names
mysamples = c("TDTD2_IRF9KO_BMDM__IFNg_1h30min_ATAC_2",
              "TDTD2_IRF9KO_BMDM___UT_ATAC_2",
              "TDTD2_WT_BMDM__IFNb_1h30min_ATAC_2") 

mysamples_ad = ad[ad$sample_name..automated. %in% mysamples,];length(mysamplesdata$group_ID); levels(mysamples_ad$genotype)
# drop unused factor levels
mysamples_ad = lapply(mysamples_ad, function(x) if(is.factor(x)) factor(x) else x)







# trying to do pvca to see influence of batch effects
library(rafalib)
library(pvca)
library(lme4)
# mmdge to Eset
mmeset=ExpressionSet(assayData = cpm(mmdge,log=TRUE), 
                     phenoData = AnnotatedDataFrame(mmdge$samples))
a <- pvcaBatchAssess_2(mmeset,batch.factors = c("genotype","date_of_collection"),threshold = 0.01)
mypar()
bp <- barplot(a$dat,  xlab = "Effects",
              ylab = "Weighted average proportion variance", ylim= c(0,1.1),
              col = c("blue"), las=2, main="PVCA estimation bar chart")
axis(1, at = bp, labels = a$label, xlab = "Effects", cex.axis = 0.5, las=2)
values = a$dat
new_values = round(values , 3)
text(bp,a$dat,labels = new_values, pos=3, cex = 0.8)

pvcaBatchAssess_2 <- function (abatch, batch.factors, threshold) {
  theDataMatrix <- exprs(abatch)
  dataRowN <- nrow(theDataMatrix)
  dataColN <- ncol(theDataMatrix)
  theDataMatrixCentered <- matrix(data = 0, nrow = dataRowN, 
                                  ncol = dataColN)
  theDataMatrixCentered_transposed = apply(theDataMatrix, 1, 
                                           scale, center = TRUE, scale = FALSE)
  theDataMatrixCentered = t(theDataMatrixCentered_transposed)
  theDataCor <- cor(theDataMatrixCentered)
  eigenData <- eigen(theDataCor)
  eigenValues = eigenData$values
  ev_n <- length(eigenValues)
  eigenVectorsMatrix = eigenData$vectors
  eigenValuesSum = sum(eigenValues)
  percents_PCs = eigenValues/eigenValuesSum
  expInfo <- pData(abatch)[, batch.factors]
  exp_design <- as.data.frame(expInfo)
  expDesignRowN <- nrow(exp_design)
  expDesignColN <- ncol(exp_design)
  my_counter_2 = 0
  my_sum_2 = 1
  for (i in ev_n:1) {
    my_sum_2 = my_sum_2 - percents_PCs[i]
    if ((my_sum_2) <= threshold) {
      my_counter_2 = my_counter_2 + 1
    }
  }
  if (my_counter_2 < 3) {
    pc_n = 3
  } else {
    pc_n = my_counter_2
  }
  pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN * 
                                               pc_n), ncol = 1)
  mycounter = 0
  for (i in 1:pc_n) {
    for (j in 1:expDesignRowN) {
      mycounter <- mycounter + 1
      pc_data_matrix[mycounter, 1] = eigenVectorsMatrix[j, 
                                                        i]
    }
  }
  AAA <- exp_design[rep(1:expDesignRowN, pc_n), ]
  Data <- cbind(AAA, pc_data_matrix)
  variables <- c(colnames(exp_design))
  for (i in 1:length(variables)) {
    Data$variables[i] <- as.factor(Data$variables[i])
  }
  op <- options(warn = (-1))
  effects_n = expDesignColN + choose(expDesignColN, 2) + 1
  randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  model.func <- c()
  index <- 1
  for (i in 1:length(variables)) {
    mod = paste("(1|", variables[i], ")", sep = "")
    model.func[index] = mod
    index = index + 1
  }
  for (i in 1:(length(variables) - 1)) {
    for (j in (i + 1):length(variables)) {
      mod = paste("(1|", variables[i], ":", variables[j], 
                  ")", sep = "")
      model.func[index] = mod
      index = index + 1
    }
  }
  function.mods <- paste(model.func, collapse = " + ")
  for (i in 1:pc_n) {
    y = (((i - 1) * expDesignRowN) + 1)
    funct <- paste("pc_data_matrix", function.mods, sep = " ~ ")
    Rm1ML <- lmer(funct, Data[y:(((i - 1) * expDesignRowN) + 
                                   expDesignRowN), ], REML = TRUE, verbose = FALSE, 
                  na.action = na.omit)
    randomEffects <- Rm1ML
    randomEffectsMatrix[i, ] <- c(unlist(VarCorr(Rm1ML)), 
                                  resid = sigma(Rm1ML)^2)
  }
  effectsNames <- c(names(getME(Rm1ML, "cnms")), "resid")
  randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, 
                                     ncol = effects_n)
  for (i in 1:pc_n) {
    mySum = sum(randomEffectsMatrix[i, ])
    for (j in 1:effects_n) {
      randomEffectsMatrixStdze[i, j] = randomEffectsMatrix[i, 
                                                           j]/mySum
    }
  }
  randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, 
                                      ncol = effects_n)
  for (i in 1:pc_n) {
    weight = eigenValues[i]/eigenValuesSum
    for (j in 1:effects_n) {
      randomEffectsMatrixWtProp[i, j] = randomEffectsMatrixStdze[i, 
                                                                 j] * weight
    }
  }
  randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
  randomEffectsSums <- colSums(randomEffectsMatrixWtProp)
  totalSum = sum(randomEffectsSums)
  randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, 
                                         ncol = effects_n)
  for (j in 1:effects_n) {
    randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum
  }
  return(list(dat = randomEffectsMatrixWtAveProp, label = effectsNames))
}


pheatmap(cpm(mmdge))
#cor_mat=cor(cpm(mmdge,log=F))
cor_mat=cor(cpm(stat3dge,log=TRUE))
#colnames(cor_mat)=gsub('RNA_|_SPN_H_H_RNA_',"",colnames(cor_mat),perl=T)
rownames(cor_mat)=colnames(cor_mat) 
library(corrplot)

corrplot(cor_mat)

corrplot.mixed(cor_mat,tl.pos = "lt",tl.col="black",tl.cex=0.75,cl.lim=c(-1,1.0))
# correlation plots
#transforming pearson correlation into distance: d=sqrt(2*(1-cor)), could also use d=(1-cor)/2
#hc = hclust(as.dist(sqrt(2*(1-cor_mat))))
hc = hclust(as.dist(1-cor_mat))
pheatmap(cor_mat,cluster_rows=hc,cluster_cols=hc,scale = "none")
m.title <- "CPM Comparisons"
aes.xlab <- bquote("log"["10"] ~ "(CPM)")
aes.ylab <- bquote("log"["10"] ~ "(CPM)")
ggpairs( log10(dat + 1),
         title = m.title,
         xlab = aes.xlab,
         ylab = aes.ylab,
         columns = seq_len(ncol(dat)),
         lower = list(continuous = "points") ) + theme_bw()
ggpairs( as.data.frame(log10(cpm(mmdge) + 1/10^6)),
         title = m.title,
         xlab = aes.xlab,
         ylab = aes.ylab,
         columns = seq_len(ncol(mmdge)),
         lower = list(continuous = my_bin )) + theme_bw()

my_bin <- function(data, mapping, ..., low = "#132B43", high = "#56B1F7") {
  ggplot(data = data, mapping = mapping) +
    geom_hex(...) +
    scale_fill_gradient(low = low, high = high)
}


#dist_mat <- sqrt(2*(1-cor_mat)) # magic step, not necessary - aequivalence between sclaed euclidian and correlation distancce?
dist_mat <- (1-cor_mat)
coords <- cmdscale(dist_mat, k=2)
ggplot(as.data.frame(coords),aes(x = coords[,1], y = coords[,2],colour=stat3dge$samples$group,shape=stat3dge$samples$date_of_collection)) + geom_point() 



pca
a<- get_pca_var(pca)


# looking at batch effects using gPCA (guided PCA)

perform_gPCA <- function(dgeobj,batch="group",topvar=1000){
  require("gPCA")
  gPCA.out <- gPCA.batchdetect(t(get_top_var(cpm(dgeobj,log=TRUE),topvar)),batch=as.numeric(dgeobj$samples[,batch]))
  old.par=par()
  gDist(gPCA.out)
  PCplot(gPCA.out,ug="guided",type="1v2")
  PCplot(gPCA.out,ug="unguided",type="1v2")
  #PCplot(gPCA.out,ug="guided",type="comp",npcs=3)
  par(old.par)
  CumulativeVarPlot(gPCA.out,ug="unguided",col="blue")
  CumulativeVarPlot(gPCA.out,ug="guided",col="blue")
  return(gPCA.out)
}

library(topGO)
# perform GO analysis
# top_var genes
library(org.Mm.eg.db)
top_row_var <- get_topvar_by_group(cpm(mmdge,log=T),mmdge$samples$group,topvar = 100)
allgenes=mmdge$genes$ENSEMBL
allgenes=as.factor(as.integer(allgenes %in% names(top_row_var$COMBINED)))
names(allgenes)=mmdge$genes$ENSEMBL
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = allgenes,
              nodeSize = 10,
              annot = annFUN.org, 
              mapping="org.Mm.eg.db",
              ID = "ensembl") 
resultTopGO.elim <- runTest(GOdata, algorithm = "elim", statistic = "Fisher" )
resultTopGO.parentchild <- runTest(GOdata, algorithm = "parentchild", statistic = "Fisher" )
resultTopGO.weight <- runTest(GOdata, algorithm = "weight", statistic = "Fisher" )
resultTopGO.classic <- runTest(GOdata, algorithm = "classic", statistic = "Fisher" )
GenTable(GOdata,elim=resultTopGO.elim,fis=resultTopGO.classic)
GenTable(GOdata,parenchil=resultTopGO.parentchild)
GenTable(GOdata,wei=resultTopGO.weight)
showSigOfNodes(GOdata, score(resultTopGO.classic), firstSigNodes = 5, useInfo = 'def')

res<-goana(as.character(mmdge$genes[top_row_var$COMBINED,c("UniqueEntrez")]),universe=as.character(mmdge$genes$UniqueEntrez),species="Mm")
topGO(res)

library("GOsummaries")
# annotation table not allowed to be to big/have NAs and all entries need to be character not factor - stupid tool
a=c("group","cell_number")
gosums=gosummaries(pca_stat3,organism="mmusculus",components = 1:2,
                   annotation=data.frame(lapply(stat3dge$samples[a],as.character),
                                         row.names = rownames(stat3dge$samples),stringsAsFactors =FALSE))
plot(gosums,classes="group")





# to do:
# try batch removal ruvseq svaseq combat?
# try remove batch effect
library("RUVSeq")
design <- model.matrix(~ mmdge$samples$group)
y <- DGEList(counts=mmdge$counts, group=mmdge$samples$group)
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")
#set <- newSeqExpressionSet(mmdge$counts)
set4 <- RUVr(mmdge$counts, rownames(mmdge), k=1, res)
set4$W
stripchart(set4$W ~ mmdge$samples$cell_number, vertical = TRUE)
boxplot(log(set4$normalizedCounts+0.11))
pca_ruvr=nicePCA(cpm(set4$normalizedCounts,log=T),annos = mmdge$samples,shape_by="cell_number")
pca_mmdge=nicePCA(cpm(mmdge,log=T),shape_by="cell_number")
library("sva")
?svaseq
mod  <- model.matrix(~ mmdge$samples$group,data=mmdge$samples)
mod0 <- model.matrix(~   1,data=mmdge$samples)
svseq <- svaseq(mmdge$counts, mod, mod0, n.sv = 1)
stripchart(svseq$sv ~ mmdge$samples$cell_number, vertical = TRUE)


require("BatchQC")
batch=mmdge$samples$date_of_collection
condition=mmdge$samples$group
batchQC(mmdge$counts, batch=batch, condition=condition, 
        report_file="batchqc_report.html", report_dir=".", 
        report_option_binary="111111111",
        view_report=FALSE, interactive=TRUE, batchqc_output=TRUE)


### fit an lm with batch to all genes simultaniously
fit_lm_over_all <- function(dgeobj,batch,group){
  expr=cpm(dgeobj,log=T)
  batches=rep(dgeobj$samples[,batch],each=nrow(expr))
  groups=rep(dgeobj$samples[,group],each=nrow(expr))
  y=as.vector(expr)
  lmdf=data.frame(y=y,batch=batches,group=groups)
  #batches=dgeobj$samples[,batch]
  #groups=dgeobj$samples[,group]
  #mod=model.matrix(~ batches + groups)
  #lm_fita=lm.fit(mod,t(expr))
  lm_fit=lm(y~batch+group,data = lmdf)
  return(lm_fit)
}


lmfit_overall_mmdge=fit_lm_over_all(mmdge,"cell_number","group")
### fit a gene wise model and extract some stuff

fit_lm_gene_wise <- function(dgeobj,batch,group){
  expr=cpm(dgeobj,log=T)
  batches=dgeobj$samples[,batch]
  groups=dgeobj$samples[,group]
  lm_fits=lapply(1:nrow(expr), function(x) return(lm(expr[x,]  ~ batches + groups)) )
  return(lm_fits)
}
lmfits_gw_mmdge=fit_lm_gene_wise(mmdge,"cell_number","group")

coef_ests <- sapply(lmfits_gw_mmdge,function(i) summary(i)$coefficients[2,1] )
var(coef_ests)



