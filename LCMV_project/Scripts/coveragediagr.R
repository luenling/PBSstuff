setwd("/Volumes/Temp/LCMV/Data/")
setwd("/Volumes//vetgrid01/LCMV_project//Data/BSF_0176_H325VBBXX_5_samples/Coverages/")
virus.covs=read.table("BSF_0176_all.coverages",header=T)
levels(virus.covs$CHROM)=c("S","L")
primers=read.table("/Volumes//vetgrid01/LCMV_project/References/primers.bed")
colnames(primers)=c("CHROM","Start","End","Primer","MQ","Strand","Flag","CIGAR","MateCHR")
# bed format 0 based open intervals:
primers$Start = primers$Start + 1
levels(primers$CHROM)
levels(primers$CHROM)=c("S","L")

l.covs=virus.covs[virus.covs$CHROM=="L",2:ncol(virus.covs)]
s.covs=virus.covs[virus.covs$CHROM=="S",2:ncol(virus.covs)]
l.covs$windows=cut(l.covs$BPS,seq(0,max(l.covs$BPS),by=50))
s.covs$windows=cut(s.covs$BPS,seq(0,max(s.covs$BPS),by=50))

l.covs.50bp=aggregate(l.covs[,1:(ncol(l.covs)-1)],list(l.covs$win),mean)
s.covs.50bp=aggregate(s.covs[,1:(ncol(s.covs)-1)],list(s.covs$win),mean)
matplot(l.covs.50bp$BPS,l.covs.50bp[,3:5],type="l",xlab="segment L",ylab="coverage (50bp average)")
matplot(l.covs.50bp$BPS,l.covs.50bp[,3:ncol(l.covs.50bp)],type="l",xlab="segment L",ylab="coverage (50bp average)")

pdf("l.cov.pdf",width=8,height=6)
matplot(l.covs.50bp$BPS,l.covs.50bp[,4:ncol(l.covs.50bp)]+0.01,type="l",xlab="base position",ylab="coverage (50bp average)", log="y", main="L segment")
dev.off()
pdf("s.cov.pdf",width=8,height=6)
matplot(s.covs.50bp$BPS,s.covs.50bp[,4:ncol(s.covs.50bp)]+0.01,type="l",xlab="base position",ylab="coverage (50bp average)", log="y", main="S segment")
dev.off()

colSums(l.covs[,2:(ncol(l.covs)-1)] > 100)/nrow(l.covs)
colSums(s.covs[,2:(ncol(s.covs)-1)] > 100)/nrow(s.covs)

library(ggplot2)
library(RColorBrewer)
library("reshape2")
l.covs.50bp.melted <- melt(l.covs.50bp[,c(2,4:ncol(l.covs.50bp))],id.vars="BPS",variable.name = "sample", value.name="coverage")
s.covs.50bp.melted <- melt(s.covs.50bp[,c(2,4:ncol(s.covs.50bp))],id.vars="BPS",variable.name = "sample", value.name="coverage")

bp = ggplot(data=l.covs.50bp.melted, aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 12), rep("dashed", 11))) +
  scale_color_manual(values = c(brewer.pal(9, "Set1"), brewer.pal(6,"Accent"),brewer.pal(8, "Dark2"))) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("L Segment") 

primers$Y=10
primers$YE=10

bp + geom_vline(xintercept = primers$Start[primers$CHR == "L" ],colour="black", linetype = "dotted", linewidth=5) +
 annotate("text", x =primers$Start[primers$CHR == "L"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "L"]), size=3)
# geom_text(data= primers[primers$CHR == "L",], aes(x = Start, y = Start, label= Primer))

bpS = ggplot(data=s.covs.50bp.melted, aes(x=BPS, y=coverage, group = sample, linetype=sample,colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 12), rep("dotted", 11))) +
  scale_color_manual(values = c(brewer.pal(9, "Set1"), brewer.pal(6,"Accent"),brewer.pal(8, "Dark2"))) + theme_bw() + 
  guides(linetype=guide_legend(ncol=2)) + ggtitle("S Segment") + ggsave("ssegment.pdf",width=8,height=6)
bpS + ggsave("ssegment.pdf",width=8,height=6)
bpS +  geom_vline(xintercept = primers$Start[primers$CHR == "S" ],colour="black", linetype = "dotted", linewidth=5) +
  annotate("text", x =primers$Start[primers$CHR == "S"], y=primers$Start[primers$CHR == "S"]*0+c(9*10^6,7*10^6), label=gsub("RT_","",primers$Primer[primers$CHR == "S"]), size=3) + ggsave("ssegment_primers.pdf",width=8,height=6)


bpL = ggplot(data=l.covs.50bp.melted, aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 12), rep("dotted", 11))) +
  scale_color_manual(values = c(brewer.pal(9, "Set1"), brewer.pal(6,"Accent"),brewer.pal(8, "Dark2"))) + theme_bw() + 
  guides(linetype=guide_legend(ncol =2)) + ggtitle("L Segment") 
bpL + ggsave("lsegment.pdf",width=8,height=6)
bpL +  geom_vline(xintercept = primers$Start[primers$CHR == "L" ],colour="black", linetype = "dotted", linewidth=5) +
  annotate("text", x =primers$Start[primers$CHR == "L"], y=primers$Start[primers$CHR == "L"]*0+c(9*10^6,7*10^6), label=gsub("RT_","",primers$Primer[primers$CHR == "L"]), size=3) + ggsave("Lsegment_primers.pdf",width=8,height=6)


setwd("/Volumes/vetgrid01/LCMV_project/Variants/Varseq2/")
dat = dcast(m2data, type + index ~ variable, value.var="position")



variants=read.table("ds1000_varscan_perc_head.tab",header=T,na.strings = c("."))
dos = read.table("ds1000_varscan_DP_head.tab",header=T,na.strings = c("."))

all.var=variants[ ! is.na(rowSums(variants[,5:27])),]
fit <- princomp(sqrt(all.var[,5:27]/100), cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)
plot(fit)
library(psych)
fit <- factor.pa(all.var[,5:27], nfactors=2, rotation="varimax")
fit # print results
str(fit)
library("FactoMineR", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
result <- PCA(all.var[,5:27])
library("PopGenome")
GENOME.class <- readData("ds1000_varscan.vcf", format="VCF", include.unknown=TRUE)

fsts.varscan=read.table("../ds1000_varscan_snps.fst.tab",header=TRUE)
fsts.median=sapply(fsts.varscan[,3:ncol(fsts.varscan)], median)
fsts.mean=sapply(fsts.varscan[,3:ncol(fsts.varscan)], median)

dm.fst.mean=matrix(0,nrow=23,ncol=23)
dimnames(dm.fst.mean)=list(colnames(all.var[,5:27]),colnames(all.var[,5:27]))
a=c()
for(i in 1:22){
  for(j in (i+1):23){
    a=append(a,( (22 - i / 2 ) * (i-1) + j -1 ) )
    dm.fst.mean[j,i] =  fsts.mean[ ( (22 - i / 2 ) * (i-1) + j -1 ) ]
  }
} 
dm.fst.mean=as.dist(dm.fst.mean,diag = FALSE, upper = FALSE)
hc <- hclust(dm.fst.mean)
plot(hc)

setwd("/Volumes/vetgrid01/LCMV_project/NewData/bwa_mem_mm10/Coverages/")
virus.covs=read.table("ND_bwa_mem_mm10_all.coverages",header=T)
levels(virus.covs$CHROM)=c("S","L")
primers=read.table("/Volumes//vetgrid01/LCMV_project/References/primers.bed")
colnames(primers)=c("CHROM","Start","End","Primer","MQ","Strand","Flag","CIGAR","MateCHR")
primers$Start = primers$Start + 1
levels(primers$CHROM)
levels(primers$CHROM)=c("S","L")
l.covs=virus.covs[virus.covs$CHROM=="L",2:ncol(virus.covs)]
s.covs=virus.covs[virus.covs$CHROM=="S",2:ncol(virus.covs)]
l.covs$windows=cut(l.covs$BPS,seq(0,max(l.covs$BPS),by=50))
s.covs$windows=cut(s.covs$BPS,seq(0,max(s.covs$BPS),by=50))

l.covs.50bp=aggregate(l.covs[,1:(ncol(l.covs)-1)],list(l.covs$win),mean)
s.covs.50bp=aggregate(s.covs[,1:(ncol(s.covs)-1)],list(s.covs$win),mean)

l.covs$windows=cut(l.covs$BPS,seq(0,max(l.covs$BPS),by=10))
s.covs$windows=cut(s.covs$BPS,seq(0,max(s.covs$BPS),by=10))

l.covs.10bp=aggregate(l.covs[,1:(ncol(l.covs)-1)],list(l.covs$win),mean)
s.covs.10bp=aggregate(s.covs[,1:(ncol(s.covs)-1)],list(s.covs$win),mean)

matplot(l.covs.50bp$BPS,l.covs.50bp[,3:5],type="l",xlab="segment L",ylab="coverage (50bp average)")
matplot(l.covs.50bp$BPS,l.covs.50bp[,3:ncol(l.covs.50bp)],type="l",xlab="segment L",ylab="coverage (50bp average)")


library(ggplot2)
library(RColorBrewer)
library("reshape2")


l.covs.50bp.melted <- melt(l.covs.50bp[,c(2,3:ncol(l.covs.50bp))],id.vars="BPS",variable.name = "sample", value.name="coverage")
s.covs.50bp.melted <- melt(s.covs.50bp[,c(2,3:ncol(s.covs.50bp))],id.vars="BPS",variable.name = "sample", value.name="coverage")
l.covs.10bp.melted <- melt(l.covs.10bp[,c(2,3:ncol(l.covs.10bp))],id.vars="BPS",variable.name = "sample", value.name="coverage")
s.covs.10bp.melted <- melt(s.covs.10bp[,c(2,3:ncol(s.covs.10bp))],id.vars="BPS",variable.name = "sample", value.name="coverage")
colors = colorRampPalette(c("black", "yellow", "blue", "green", "red"))( 14 ) 
colors=c(rev(brewer.pal(6,"Set1")),"darkgrey")

cl13p=grep("Cl13.*plasmid",levels(l.covs.50bp.melted$sample))
cl13bh=grep("Cl13.*BHK",levels(l.covs.50bp.melted$sample))
brain=grep("Brain",levels(l.covs.50bp.melted$sample))
first=seq(1,6)
serum=grep("Serum",levels(l.covs.50bp.melted$sample))

displ=brain
bp = ggplot(data=l.covs.10bp.melted[l.covs.10bp.melted$sample %in% levels(l.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 14), rep("dashed", 14))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("L Segment") 
bp + geom_vline(xintercept = primers$Start[primers$CHR == "L" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "L"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "L"]), size=3) + ggsave("lsegment_brains.pdf",width=8,height=6)

bps = ggplot(data=s.covs.10bp.melted[s.covs.10bp.melted$sample %in% levels(s.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 14), rep("dashed", 14))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("S Segment") 
bps + geom_vline(xintercept = primers$Start[primers$CHR == "S" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "S"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "S"]), size=3) + ggsave("ssegment_brains.pdf",width=8,height=6)

displ=cl13p
bp = ggplot(data=l.covs.10bp.melted[l.covs.10bp.melted$sample %in% levels(l.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 7), rep("dashed", 7))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("L Segment") 
bp + geom_vline(xintercept = primers$Start[primers$CHR == "L" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "L"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "L"]), size=3) + ggsave("lsegment_cl13_plasmid.pdf",width=8,height=6)

bps = ggplot(data=s.covs.10bp.melted[s.covs.10bp.melted$sample %in% levels(s.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 7), rep("dashed", 7))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("S Segment") 
bps + geom_vline(xintercept = primers$Start[primers$CHR == "S" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "S"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "S"]), size=3) + ggsave("ssegment_cl13_plasmid.pdf",width=8,height=6)


displ=cl13bh
bp = ggplot(data=l.covs.10bp.melted[l.covs.10bp.melted$sample %in% levels(l.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 7), rep("dashed", 7))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("L Segment") 
bp + geom_vline(xintercept = primers$Start[primers$CHR == "L" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "L"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "L"]), size=3) + ggsave("lsegment_cl13_bhk21.pdf",width=8,height=6)

bps = ggplot(data=s.covs.10bp.melted[s.covs.10bp.melted$sample %in% levels(s.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 7), rep("dashed", 7))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("S Segment") 
bps + geom_vline(xintercept = primers$Start[primers$CHR == "S" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "S"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "S"]), size=3) + ggsave("ssegment_cl13_bhk21.pdf",width=8,height=6)

displ=first

bp = ggplot(data=l.covs.10bp.melted[l.covs.10bp.melted$sample %in% levels(l.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 7), rep("dashed", 7))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("L Segment") 
bp + geom_vline(xintercept = primers$Start[primers$CHR == "L" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "L"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "L"]), size=3) + ggsave("lsegment_0_ARM_B16.pdf",width=8,height=6)

bps = ggplot(data=s.covs.10bp.melted[s.covs.10bp.melted$sample %in% levels(s.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 7), rep("dashed", 7))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("S Segment") 
bps + geom_vline(xintercept = primers$Start[primers$CHR == "S" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "S"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "S"]), size=3) + ggsave("ssegment_0_ARM_B16.pdf",width=8,height=6)

displ=serum

bp = ggplot(data=l.covs.10bp.melted[l.covs.10bp.melted$sample %in% levels(l.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 7), rep("dashed", 7))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("L Segment") 
bp + geom_vline(xintercept = primers$Start[primers$CHR == "L" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "L"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "L"]), size=3) + ggsave("lsegment_serum.pdf",width=8,height=6)

bps = ggplot(data=s.covs.10bp.melted[s.covs.10bp.melted$sample %in% levels(s.covs.10bp.melted$sample)[displ],], aes(x=BPS, y=coverage, group = sample, linetype=sample, colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 7), rep("dashed", 7))) +
  scale_color_manual(values = c(colors,colors)) + theme_bw() + 
  theme(legend.direction ="horizontal",legend.position = "bottom") + guides(col=guide_legend(nrow=3)) + ggtitle("S Segment") 
bps + geom_vline(xintercept = primers$Start[primers$CHR == "S" ], aes(colour="black", linetype = "dotted", linewidth=5)) +
  annotate("text", x =primers$Start[primers$CHR == "S"], y=9*10^6, label=gsub("RT_","",primers$Primer[primers$CHR == "S"]), size=3) + ggsave("ssegment_serum.pdf",width=8,height=6)



primers$Y=10
primers$YE=10

# geom_text(data= primers[primers$CHR == "L",], aes(x = Start, y = Start, label= Primer))

bpS = ggplot(data=s.covs.50bp.melted, aes(x=BPS, y=coverage, group = sample, linetype=sample,colour = sample)) +
  geom_line() + scale_y_log10() + scale_linetype_manual(values = c(rep("solid", 6), rep("dotted", 6))) +
  scale_color_manual(values = c(brewer.pal(6, "Set1"), brewer.pal(6,"Accent"),brewer.pal(8, "Dark2"))) + theme_bw() + 
  guides(linetype=guide_legend(ncol=2)) + ggtitle("S Segment") + ggsave("ssegment.pdf",width=8,height=6)
