# explorative plots of QC data for RNAseq mapper comparison
# QC data from N. Fortelny, 29.08.19
library(ggplot2)
library("pheatmap")


setwd("~/Data/SFB61/")
# load table
qc_data <- read.table("red_QC_Merged_Align_comp_29082019.txt",header=TRUE,sep = "\t")
rownames(qc_data)<-qc_data$sample_name

# just get relevant field names
grep("HISAT|RNASTAR",names(qc_data),value=T)
# Plots to do
# always STAR no MM agains Hisat all and Hisat no MM?

# log counts STAR / log counts Hisat
ggplot(data=qc_data, aes(x=HISAT_NoMultimappersreads.exons,y=RNASTAR.Exome_reads)) + geom_point() + 
  scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") 


ggplot(data=qc_data, aes(x=HISAT_NoMultimappersreads.fullTranscripts,y=RNASTAR.Transcriptome_reads)) + geom_point() + 
  scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") 

# scatter plots relative to all reads (FASTQC counts)
ggplot(data=qc_data, 
       aes(x=HISAT_NoMultimappersreads.fullTranscripts/RNA.Trimmed_reads,y=RNASTAR.Transcriptome_reads/RNA.Trimmed_reads)) +
  geom_point(aes(colour = RNA.Trimmed_reads) ) + scale_colour_gradient2(low = "red", mid="green", high = "darkblue",midpoint=3e+07) +
  geom_abline(intercept=0,slope=1,color="grey",linetype="dashed")
ggplot(data=qc_data, 
       aes(x=HISAT_NoMultimappersreads.exons/RNA.Trimmed_reads,y=RNASTAR.Exome_reads/RNA.Trimmed_reads)) +
  geom_point(aes(colour = RNA.Trimmed_reads) ) + scale_colour_gradient2(low = "red", mid="green", high = "darkblue",midpoint=3e+07) +
  geom_abline(intercept=0,slope=1,color="grey",linetype="dashed")

# histogram of percent differences
mean_counts_up_down <-function(x,y,z) {
  # takes three vectors of counts of reads with the same length (mapped by mapper a, b and total reads used)
  # returns list with : a>b=c( #a>b, mean((a-b)/total) for (a>b)) ,  b>a=c( #, mean((b-a)/total) for (b>a), a=b = c( #a=b, 0)
  fxy <- (x-y)/z
  xgb <- x > y
  xly <- x < y
  return(list(agb=c(sum(xgb),mean(fxy[xgb])),bga=c(sum(xly),mean(fxy[xly])) , aeb=c(length(x)-sum(xgb)-sum(xly),0)))
}

a <- with(qc_data,wmean_counts_up_down(HISAT_NoMultimappersreads.exons,RNASTAR.Exome_reads,RNA.Trimmed_reads))
ggplot(data=qc_data, 
       aes(x=(HISAT_NoMultimappersreads.exons-RNASTAR.Exome_reads)/RNA.Trimmed_reads*100)) +
  geom_histogram(bins=25,fill="pink",color="grey")+ geom_vline(xintercept = 0, linetype="dashed",color="grey") + xlab(label="Hisat2 (noMM) - STAR [ % total reads mapped to exons]") + ylab(label="libraries") + 
  geom_text(x=-3,y=22.5,label=paste("STAR better","\n","samples:",a$bga[1],"\n","avg(diff):",round(a$bga[2]*100,2),"%"),hjust = 0) +
  geom_text(x=0.25,y=22.5,label=paste("Hisat2 better","\n","samples:",a$agb[1],"\n","avg(diff):",round(a$agb[2]*100,2),"%"),hjust = 0) 


a <- with(qc_data,wmean_counts_up_down(HISAT_NoMultimappersreads.fullTranscripts,RNASTAR.Transcriptome_reads,RNA.Trimmed_reads))
ggplot(data=qc_data, 
       aes(x=(HISAT_NoMultimappersreads.fullTranscripts-RNASTAR.Transcriptome_reads)/RNA.Trimmed_reads*100)) +
  geom_histogram(bins=25,fill="pink",color="grey")+ geom_vline(xintercept = 0, linetype="dashed",color="grey") + xlab(label="Hisat2 (noMM) - STAR [ % total reads mapped to transcripts]") + ylab(label="libraries") + 
  geom_text(x=-4.5,y=22.5,label=paste("STAR better","\n","samples:",a$bga[1],"\n","avg(diff):",round(a$bga[2]*100,2),"%"),hjust = 0) +
  geom_text(x=0.25,y=22.5,label=paste("Hisat2 better","\n","samples:",a$agb[1],"\n","avg(diff):",round(a$agb[2]*100,2),"%"),hjust = 0) 


# read counts on gene level
# read matrices
hisat_counts <- read.csv("RNA.exome.Counts.HISAT.noMultiMappers.csv", header=TRUE)
star_counts <- read.csv("RNA.exome.Counts.STAR.noMultiMappers.csv", header=TRUE)
gene_annos <- read.table("RNA.exome.Probes.STAR.noMultiMappers.csv", header = TRUE)
rownames(gene_annos) <- gene_annos$ensG
gene_annos$probe<-NULL
gene_annos$gene_version<-NULL
gene_annos$gene_source<-NULL

# more samples in star_counts
star_counts <- star_counts[rownames(hisat_counts),colnames(hisat_counts)]
# get rid of zero counts
h_zero <- rowSums(hisat_counts) == 0 ; sum(h_zero)
s_zero <- rowSums(star_counts) == 0 ; sum(s_zero)
hs_zero <- s_zero & h_zero ; sum(hs_zero)
star_counts <- star_counts[!hs_zero,]
hisat_counts <- hisat_counts[!hs_zero,]

# check whether counts are identical
all(colSums(star_counts)[rownames(qc_data)]-qc_data$RNASTAR.Exome_reads == 0) # equals zero
all(colSums(hisat_counts)[rownames(qc_data)]-qc_data$HISAT_NoMultimappersreads.exons == 0) # equals zero

library(edgeR)
star_cpms <- cpm(star_counts) 
hisat_cpms <- cpm(hisat_counts) 
# or do star_rpms <- t(t(star_counts)/col)
#star_rpms <- t(t(star_counts)/colSums(star_counts))*1E06 
#hisat_rpms <- t(t(hisat_counts)/colSums(hisat_counts))*1E06 
#all.equal(star_rpms,cpm(star_counts), na.rm = T)
# remove all rows with median below 15 counts in both hisat and star



med_gt_25 <- (apply(hisat_counts,1,median) >=25) | (apply(star_counts,1,median) >=25)
med_gt_50 <- (apply(hisat_counts,1,median) >=50) | (apply(star_counts,1,median) >=50)
med_gt_100 <- (apply(hisat_counts,1,median) >=100) | (apply(star_counts,1,median) >=100)
sum(med_gt_25) # 9662
sum(med_gt_50) # 8728
sum(med_gt_100) # 7611

star_counts_50 <- star_counts[med_gt_50,]
hisat_counts_50 <- hisat_counts[med_gt_50,]
sh_rel_diff_counts_50 <- (hisat_counts_50 - star_counts_50)/pmax(star_counts_50,hisat_counts_50)
sh_rel_diff_counts_50[is.na(sh_rel_diff_counts_50)] <- 0
most_diff <- order(abs(rowSums(sh_rel_diff_counts_50)),decreasing = TRUE)
# get maximal percent change
all_maxperc <- apply(abs(round(sh_rel_diff_counts_50[most_diff,],digits=3))*100,1,max)
all_maxminperc <- apply(round(sh_rel_diff_counts_50[most_diff,],digits=3)*100,1,
                        function(x) {if (max(x) >= abs(min(x)) ) max(x) else min(x) })
all_median_perc <- apply(round(sh_rel_diff_counts_50[most_diff,],digits=3)*100,1,median)
hist(all_maxminperc)
hist(all_median_perc)

ggplot(data=as.data.frame(all_median_perc), aes(x=all_median_perc)) +
  geom_histogram(bins=50,fill="pink",color="grey")+ geom_vline(xintercept = 0, linetype="dashed",color="grey") + labs(title="Median percent difference between counts per gene (Hisat2-STAR)",x="median(Hisat2 (noMM) - STAR [ count difference in % of higher count])", y = "number of genes")  
all_perc_diff <- unlist(sh_rel_diff_counts_50,use.names = F)
ggplot(data=as.data.frame(all_perc_diff), aes(x=all_perc_diff,y=stat(count/length(all_perc_diff)))) +
  geom_histogram(bins=50,fill="pink",color="grey")+ geom_vline(xintercept = 0, linetype="dashed",color="grey") + labs(title="Percent difference between counts per gene and sample (Hisat2-STAR)",x="Hisat2 (noMM) - STAR [ count difference in % of higher count]", y = "fraction of genes x sample")

rnames <- rep("",10000)
rnames[seq.int(50,10000,by = 50)] <- as.character(seq.int(50,10000,by = 50))
rnames[1] <- "1"
md <- 700
png(filename = "heatmap_700genes_highest_perc_difference_hisat_star",width = 900,height=700)
pheatmap(sh_rel_diff_counts_50[most_diff[1:md],],scale="none",cluster_rows = FALSE,cluster_cols = FALSE,labels_col=rnames[1:ncol(sh_rel_diff_counts_50)],labels_row = rnames[1:md],main="700 genes with highest overall difference (Hisat2 - Star) over all 236 samples")
dev.off()



sum(abs(all_median_perc)>50) # 397
sum(abs(all_median_perc)>20) # 541
sum(abs(all_median_perc)>10) # 664
#### 

sh_rel_diff_counts_50[most_diff[1:10],1:10]
star_counts_50[most_diff[1:10],1:10]
hisat_counts_50[most_diff[1:10],1:10]
cat(rownames(sh_rel_diff_counts_50)[most_diff[1:100]],sep="\n")
get_gene_biotypes <- function(median_perc,ga = gene_annos) {
  biotype = ga[names(median_perc),"gene_biotype"]
  biotype=droplevels(biotype)
  median_perc <- sign(median_perc)
  ct <- table(median_perc,biotype)
  return(ct[,order(colSums(ct),decreasing = TRUE)])
}

a <- get_gene_biotypes(all_median_perc[1:700])
op <- par(mar = c(8,4,4,2) + 0.1)
plt <- barplot(a,beside=T, xaxt="n",ylab="number of genes for biotype",col=c("blue","yellow","red"),main="700 genes with highest difference between Hisat2 and Star")
text(plt[2,], par("usr")[3], colnames(a), xpd=TRUE, srt=45,adj = c(1.1,1.1), cex=0.7)
legend("topright", legend = c("less", "same", "more"), title = "HISAT2 counts",
       fill = c("blue","yellow","red"), bty=F)

par(op) 
a <- get_gene_biotypes(all_median_perc[1:350])
op <- par(mar = c(8,4,4,2) + 0.1)
plt <- barplot(a,beside=T, xaxt="n",ylab="number of genes for biotype",col=c("blue","red"),main="350 genes with highest difference between Hisat2 and Star")
text(plt[2,], par("usr")[3], colnames(a), xpd=TRUE, srt=45,adj = c(1.1,1.1), cex=0.7)
legend("topright", legend = c("less", "more"), title = "HISAT2 counts",
       fill = c("blue","red"), bty=F)

par(op) 


sort(table(gene_annos[rownames(sh_rel_diff_counts_50)[most_diff[1:500]],4]))
pheatmap(sh_rel_diff_counts_50[most_diff[1:750],],scale="none",cluster_rows = FALSE,cluster_cols = FALSE,labels_col=1:ncol(sh_rel_diff_rpms),show_rownames= FALSE)

pca <- prcomp(t( sh_rel_diff_counts_50[most_diff[1:750],] ) , center = TRUE, scale. = FALSE) 
barplot(pca$sdev^2/sum(pca$sdev^2),ylab="% variance explained",names.arg = colnames(pca$x),las=2)
ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC2)) + geom_point()
ggplot(data=as.data.frame(pca$x),aes(x=PC2,y=PC3)) + geom_point()

pca_cl1 <-rownames(pca$x[with(as.data.frame(pca$x), PC1 < 0 & PC2 > 0),])
pca_cl2 <-rownames(pca$x[with(as.data.frame(pca$x), PC1 > 0 ),])
pca_cl3 <-rownames(pca$x[with(as.data.frame(pca$x), PC1 < 0 & PC2 < 0),])
pca_cls <- 1:nrow(pca$x); names(pca_cls) <- rownames(pca$x); 
pca_cls[pca_cl1] <- 1
pca_cls[pca_cl2] <- 2
pca_cls[pca_cl3] <- 3



ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC2,colour=as.factor(pca_cls))) + geom_point()
km3<-kmeans(scale(t(sh_rel_diff_counts_50[most_diff[1:750],]),scale = T),3)
pam3<-pam(scale(t(sh_rel_diff_counts_50[most_diff[1:750],]),scale = T),3)
pam3$isolation
ggplot(data=as.data.frame(pca$x),aes(x=PC1,y=PC2,colour=as.factor(km3$cluster))) + geom_point()

pheatmap(sh_rel_diff_counts_50[most_diff[1:750],],scale="none",cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = FALSE,show_rownames= FALSE ,annotation_col = data.frame(km=as.factor(km3$cluster),pam=as.factor(pam3$clustering),pcacl=as.factor(pca_cls)) )


rclust=hclust(as.dist(1-cor(sh_rel_diff_counts_50[most_diff[1:550],])))
pheatmap(sh_rel_diff_counts_50[most_diff[1:750],],scale="none",cluster_cols = TRUE,show_colnames = FALSE,show_rownames= FALSE ,annotation_col = data.frame(pam=as.factor(pam3$clustering),cell_type=ad$cell_type,ad$tissue, ad$genotype, ad$method_ip) )

library(fpc)
# has got very nice cluster applications
plotcluster(sh_rel_diff_counts_50[most_diff[1:750],],pam4_g$clustering)
plotcluster(sh_rel_diff_counts_50[most_diff[1:750],],pam2_g$clustering)

names(pam2_g$clustering[pam2_g$clustering==1])

library("GOsummaries")
km3_g<-kmeans(sh_rel_diff_counts_50[most_diff[1:750],],3)
gs <- gosummaries(km3_g,organism = "mmusculus")
plot(gs)
gs




pam3_g<-pam(sh_rel_diff_counts_50[most_diff[1:750],],3)
plot(silhouette(pam3_g,dist(sh_rel_diff_counts_50[most_diff[1:750],])))
pam2_g<-pam(sh_rel_diff_counts_50[most_diff[1:750],],2)
plot(silhouette(pam2_g,dist(sh_rel_diff_counts_50[most_diff[1:750],])))


#read sample annotation table
annotation = read.csv("Annotation.csv")
ad = subset(annotation,grepl("RNA",library))
samples=intersect(ad$sample_name,names(hisat_counts))
ad <- ad[ad$sample_name %in% samples,]
rownames(ad) <- ad$sample_name
ad <- ad[colnames(hisat_counts),]
ad<-droplevels.data.frame(ad)
to_red <- unlist(lapply(ad, function(x) { (length(unique(x)) > 1 )}) )
to_red_2 <- unlist(lapply(ad, function(x) { (length(unique(x)) > 1) & (length(unique(x)) < 30)} ) )
ad_red <- ad[to_red_2]

boxplot(ad_red[pam3$clustering == 3,],las=2,ylim=c(0,25),border=rgb(1, 0, 0, alpha = 0.5), col=rgb(1, 0, 0, alpha = 0.5))
par(new=T)
boxplot(ad_red[pam3$clustering == 2,],las=2,ylim=c(0,25),border=rgb(0,1, 0, alpha = 0.5), col=rgb(0,1, 0, alpha = 0.5))
par(new=T)
boxplot(ad_red[pam3$clustering == 1,],las=2,ylim=c(0,25),border=rgb(0,0,1, alpha = 0.5), col=rgb(0,0,1, alpha = 0.5))

table(droplevels(ad_red$genotype[pam3$clustering == 1])) # WTN
table(droplevels(ad_red$genotype[pam3$clustering == 2])) # STAT1KO   STAT1AA   STAT1BB   WTN       TYK2KE    TYK2CMV   STAT1FLOX STAT1VAV  STAT2KO
table(droplevels(ad_red$genotype[pam3$clustering == 3])) # WTJ         WTN         STAT5BN642H STAT5VAV    STAT5FLOX   STAT3VAV    STAT3FLOX   TYK2KE      TYK2CMV     TYK2KO 
comp=c("cell_type","genotype","tissue")
table(droplevels(ad_red[pam3$clustering == 3,comp]))
table(droplevels(ad_red[pam3$clustering == 2,comp]))
table(droplevels(ad_red[pam3$clustering == 1,comp]))

str(ad)
library(cluster)
library(NbClust)
nb_test<-NbClust(sh_rel_diff_counts_50[most_diff[1:750],], distance = "euclidean", min.nc=2, max.nc=8, 
        method = "kmeans", index = "all")

ag_cl<-agnes(pca$x[,c(1,2,3)])
#plot(ag_cl)
pltree(ag_cl,labels = pca_cls)


star_cpms <- star_cpms[med_gt_15,]
hisat_cpms <- hisat_cpms[med_gt_15,]
star_cpms <- star_cpms[med_gt_100,]
hisat_cpms <- hisat_cpms[med_gt_100,]

sh_rel_diff_cpms <- (star_cpms - hisat_cpms)/pmax(star_cpms,hisat_cpms)
sh_rel_diff_cpms[is.na(sh_rel_diff_cpms)] <- 0
most_diff <- order(abs(rowSums(sh_rel_diff_cpms)),decreasing = TRUE)
# get maximal percent change
all_maxperc <- apply(abs(round(sh_rel_diff_cpms[most_diff,],digits=3))*100,1,max)
all_maxminperc <- apply(round(sh_rel_diff_cpms[most_diff,],digits=3)*100,1,
                        function(x) {if (max(x) >= abs(min(x)) ) max(x) else min(x) })
max(all_maxperc)
hist(all_maxminperc)
gene_annos[names(all_maxperc[all_maxperc >= 1]),]
#pheatmap(star_rpms[most_diff[1:100],] - hisat_rpms[most_diff[1:100],],scale="row",cluster_rows = FALSE)
pheatmap(sh_rel_diff_cpms[most_diff[1:1000],],scale="none",cluster_rows = FALSE,cluster_cols = FALSE,labels_col=1:ncol(sh_rel_diff_rpms),show_rownames= FALSE)
# get only gene names, one per row , in console
cat(rownames(sh_rel_diff_cpms)[most_diff[1:100]],sep="\n")

gene_annos[rownames(sh_rel_diff_rpms)[most_diff[1:100]],]
pheatmap(sh_rel_diff_rpms[most_diff[1:50],],scale="none",cluster_rows = FALSE,labels_col=1:ncol(sh_rel_diff_rpms))
mart = biomaRt::useMart(biomart="ensembl",dataset="mmusculus_gene_ensembl")
biomaRt::getBM(attributes = c("ensembl_gene_id","gene_biotype", "description", "chromosome_name"), values=names(all_maxperc[all_maxperc >= 1]), filters = "ensembl_gene_id",mart=mart)

library("pheatmap")

most_rel_diff <- order(abs(rowSums(sh_diff_counts)),decreasing = TRUE)
heatmap.2(as.matrix(sh_rel_diff_counts[most_rel_diff[1:10],]),scale="row",Rowv = NA)
pheatmap(as.matrix(sh_diff_counts[most_diff[1:10],]),scale="none",cluster_rows = FALSE)
