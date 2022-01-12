#BiocManager::install("ggtree")
#install.packages("remotes")
#remotes::install_github("GuangchuangYu/ggtree")
#BiocManager::install("treeio")
install.packages("phytools")
library(phytools)
library("ggplot2")
library("ggtree")
library("treeio")
library(ape)
library("tidytree")
library(ggrepel)
library("reshape2")

setwd("~/Data/Merima")

#sel_ort_raxml@phylo=root(sel_ort_raxml@phylo,outgroup="AJ704540.1")
annos_16S <- read.table(file="16S_annos_all.txt",header=T,sep = "\t")
annos_rpoB <- read.table(file="rpoB_annos_all.txt",header=T,sep = "\t")

sel_ort_raxml <- read.raxml("RAxML_bipartitionsBranchLabels.all_16Srna_nogap.mafft.raxml")
sel_ort_raxml <- read.raxml("RAxML_bipartitionsBranchLabels.all_16Srna_nogap.mafft.raxml")
str(sel_ort_raxml)
print(sel_ort_raxml)
show(sel_ort_raxml)
# drop the R. anati and other outgroups
new.phylo <- treeio::drop.tip(sel_ort_raxml, "Rie_anati" )
new.phylo <- treeio::drop.tip(new.phylo, c("UNSD01.1","UNSC01.1" ) )
#new.phylo <-root(new.phylo,c("UNSD01.1", "UNSC01.1"))
# plot the whole crap
ts=3
nptree <- ggtree(new.phylo) %<+% annos_16S
nptree + geom_tiplab(aes(label=Strain),align=T, size=ts) + 
  geom_text2(aes(x=branch,label=bootstrap, subset = !is.na(as.numeric(bootstrap)) & as.numeric(bootstrap) >= 101),position = position_nudge(y =0.5),size=ts) + geom_rootedge(0.001) +
  geom_treescale() + geom_text(aes(x=0.025,y=y,label=Host),size=ts) + geom_text(aes(x=0.0275,y=y,label=serovar),size=ts) + geom_text(aes(x=0.03,y=y,label=ST),size=ts)

# without R. Anati
new.phylo <- treeio::drop.tip(sel_ort_raxml, "Rie_anati" )
ts=3
nptree <- ggtree(new.phylo) %<+% annos_16S
nptree + geom_tiplab(aes(label=Strain),align=T, size=ts) + 
  geom_text2(aes(x=branch,label=bootstrap, subset = !is.na(as.numeric(bootstrap)) & as.numeric(bootstrap) >= 101),position = position_nudge(y =0.5),size=ts) + geom_rootedge(0.001) +
  geom_treescale() + geom_text(aes(x=0.11,y=y,label=Host),size=ts) + geom_text(aes(x=0.12,y=y,label=serovar),size=ts) 


# all strains + R. Anatipestifer
ort_16S_t <- ggtree(sel_ort_raxml) %<+% annos_16S
ort_16S_t + geom_tiplab(aes(label=Strain),align=T,size=ts)   +  geom_treescale() + 
   geom_text(aes(x=0.21,y=y,label=serovar),size=ts)  + geom_text(aes(x=0.18,y=y,label=Host),size=ts) + 
  geom_text2(aes(x=branch,label=bootstrap, subset = !is.na(as.numeric(bootstrap)) & as.numeric(bootstrap) >= 40),position = position_nudge(y =0.5),size=ts)



# load tree with selected 16S RNAs and no duplicated sequences, FARPER
sel_ort_raxml <- read.raxml("RAxML_bipartitionsBranchLabels.selected_FARP_16Srna_nogap_reduced.mafft.raxml")
sel_ort_raxml <- read.raxml("RAxML_bipartitionsBranchLabels.selected16S_reduced.mafft.mfa_alt.raxml")
str(sel_ort_raxml)
print(sel_ort_raxml)
show(sel_ort_raxml)
# drop the R. anati and other outgroups
#new.phylo <- treeio::drop.tip(sel_ort_raxml, c("UNSD01.1","UNSC01.1" ) )
#new.phylo <-root(new.phylo,c("UNSD01.1", "UNSC01.1"))
# plot the whole crap
sel_tibl=as_tibble(sel_ort_raxml)
outpar=sel_tibl$parent[sel_tibl$label %in% "UNSC01.1"]
outgrp=sel_tibl$node[sel_tibl$label %in% "UNSC01.1"]
outsib=sel_tibl$node[sel_tibl$parent == outpar & ! sel_tibl$node  %in% c(outpar,outgrp)]
sel_tibl$group = 0
sel_tibl$group[sel_tibl$node %in% c(outstart,outgrp,outsib)]=1
sel_tibl$group=as.factor(sel_tibl$group)
tot_branch=sum(sel_tibl$branch.length[sel_tibl$node %in% c(outgrp,outsib)])
sel_tibl$branch.length[sel_tibl$node ==  outgrp] = tot_branch*0.975/10
sel_tibl$branch.length[sel_tibl$node ==  outsib] = tot_branch*0.025/10

#outgrp=getMRCA(sel_ort_raxml@phylo,c("UNSD01.1", "UNSC01.1"))
#el_ort_raxml@phylo <- groupClade(sel_ort_raxml@phylo,outgrp)

#lbr=child(sel_ort_raxml,parent(sel_ort_raxml,outgrp))
#sel_tibl=as_tibble(sel_ort_raxml)

#sel_tibl$branch.length[sel_tibl$node %in% lbr] = sel_tibl$branch.length[sel_tibl$node %in% lbr]/100
sel_ort_alt = as.treedata(sel_tibl)

ts=3
annos_16S$seq=as.factor(c(rep("seq",3),rep("not_seq",nrow(annos_16S)-3)))
annos_16S$formLab=paste0('"',as.character(annos_16S$Strain),'"')
annos_16S$formLab[annos_16S$seq == "seq"] = paste0('bold(', annos_16S$formLab[annos_16S$seq == "seq"], ')') 
nptree <- ggtree(sel_ort_alt, aes(linetype=group)) %<+% annos_16S
#nptree$data$x = nptree$data$x - mean(nptree$data$x)
#nptree$data[nptree$data$node %in% c(outgrp), "x"] = max(nptree$data$x)
#nptree$data[nptree$data$label %in% c("UNSD01.1", "UNSC01.1"), "x"] = max(nptree$data$x)

 nptree %>% ggtree::rotate(24) %>% ggtree::rotate(29) +
  geom_tiplab(aes(label=formLab),align=T,size=ts,parse=TRUE)  +  geom_treescale() + theme(legend.position = 'none') + #geom_text(aes(x=branch,label=node)) +
  geom_text(aes(x=0.02625,y=y,label=Host),size=ts) + geom_text(aes(x=0.0285,y=y,label=serovar),size=ts) + geom_text(aes(x=0.02975,y=y,label=ST),size=ts) +
  geom_text2(aes(x=branch,label=bootstrap, subset = !is.na(as.numeric(bootstrap)) & as.numeric(bootstrap) >= 50),position = position_nudge(y =0.5),size=ts)

# load tree with selected 1geom_text(aes(x=branch,label=node)) *6S RNAs and no duplicated sequences
sel_ort_raxml <- read.raxml("RAxML_bipartitionsBranchLabels.selected16S_reduced.mafft.mfa_alt.raxml")
str(sel_ort_raxml)
print(sel_ort_raxml)
show(sel_ort_raxml)
# drop the R. anati and other outgroups
new.phylo <- treeio::drop.tip(sel_ort_raxml, c("UNSD01.1","UNSC01.1" ) )
#new.phylo <-root(new.phylo,c("UNSD01.1", "UNSC01.1"))
# plot the whole crap
ts=3
annos_16S$seq=as.factor(c(rep("seq",3),rep("not_seq",nrow(annos_16S)-3)))
annos_16S$formLab=paste0('"',as.character(annos_16S$Strain),'"')
annos_16S$formLab[annos_16S$seq == "seq"] = paste0('bold(', annos_16S$formLab[annos_16S$seq == "seq"], ')') 
nptree <- ggtree(new.phylo) %<+% annos_16S
p <- ggtree::rotate(nptree,27) #  %>%  ggtree::rotate(44)
p + geom_tiplab(aes(label=formLab),align=T,size=ts,parse=TRUE)   +  geom_treescale() + geom_rootedge(rootedge=0.0002) + 
  geom_text(aes(x=0.02,y=y,label=Host),size=ts) + geom_text(aes(x=0.022,y=y,label=serovar),size=ts) + geom_text(aes(x=0.0235,y=y,label=ST),size=ts) +
  geom_text2(aes(x=branch,label=bootstrap, subset = !is.na(as.numeric(bootstrap)) & as.numeric(bootstrap) >= 50),position = position_nudge(y =0.5),size=ts)

# rpoB tree
#rpo_raxml <- read.raxml("RAxML_bipartitionsBranchLabels.rpoB_reduced_alt.mafft_raxml")
  rpo_raxml <- read.raxml("RAxML_bipartitionsBranchLabels.rpoB_reduced_1og.mafft_raxml")
#rpo_raxml <- read.raxml("RAxML_bipartitionsBranchLabels.farper_rpo_alt.mfa.reduced_raxml")

rpo_tibl=as_tibble(rpo_raxml)
outpar=rpo_tibl$parent[rpo_tibl$label %in% "rpoB_UNSC01.1"]
outgrp=rpo_tibl$node[rpo_tibl$label %in% "rpoB_UNSC01.1"]
outsib=rpo_tibl$node[rpo_tibl$parent == outpar & ! rpo_tibl$node  %in% c(outpar,outgrp)]
rpo_tibl$group = 0
rpo_tibl$group[rpo_tibl$node %in% c(outstart,outgrp,outsib)]=1
rpo_tibl$group=as.factor(rpo_tibl$group)
tot_branch=sum(rpo_tibl$branch.length[rpo_tibl$node %in% c(outgrp,outsib)])
rpo_tibl$branch.length[rpo_tibl$node ==  outgrp] = tot_branch*0.975/5
rpo_tibl$branch.length[rpo_tibl$node ==  outsib] = tot_branch*0.025/5


#show(rpo_raxml)
# drop the R. anati and other outgroups
#new.phylo <- treeio::drop.tip(rpo_raxml, "rpoB_Riemerella_anatipestifer" )
#new.phylo <- treeio::drop.tip(new.phylo, c("rpoB_UNSD01.1","rpoB_UNSC01.1" ) )
#new.phylo <-root(new.phylo,c("UNSD01.1", "UNSC01.1"))
# plot the whole crap
ts=3
annos_rpoB$formLab=paste0('"',as.character(annos_rpoB$Strain),'"')
annos_rpoB$formLab[grep("LMG1[89]",annos_rpoB$Species)] = paste0('bold(', annos_rpoB$formLab[grep("LMG1[89]",annos_rpoB$Species)], ')') 

#outgrp=getMRCA(rpo_raxml@phylo,c("rpoB_UNSD01.1", "rpoB_UNSC01.1"))
#rpo_raxml@phylo <- groupClade(rpo_raxml@phylo,outgrp)
#lbr=child(rpo_raxml,parent(rpo_raxml,outgrp))
#rpo_tibl=as_tibble(rpo_raxml)

#rpo_tibl$branch.length[rpo_tibl$node %in% lbr] = rpo_tibl$branch.length[rpo_tibl$node %in% lbr]/100
rpo_alt = as.treedata(rpo_tibl)


#nptree <- ggtree(new_phylo) %<+% annos_rpoB
nptree <- ggtree(rpo_alt, aes(linetype=group)) %<+% annos_rpoB
#nptree$data[nptree$data$node %in% c(outgrp), "x"] = max(nptree$data$x)
#nptree$data[nptree$data$label %in% c("rpoB_UNSD01.1", "rpoB_UNSC01.1"), "x"] = max(nptree$data$x)

nptree  + geom_tiplab(aes(label=formLab),align=T, size=ts,parse=TRUE) +  theme(legend.position = 'none') +
  geom_text2(aes(x=branch,label=bootstrap, subset = !is.na(as.numeric(bootstrap)) & as.numeric(bootstrap) >= 50 & as.numeric(bootstrap) < 100),position = position_nudge(y =0.5),size=ts) + 
  geom_treescale() + geom_text(aes(x=0.24,y=y,label=Host),size=ts) + 
  geom_text(aes(x=0.26,y=y,label=serovar),size=ts) + geom_text(aes(x=0.275,y=y,label=ST),size=ts)

# all strains + R. Anatipestifer
ort_rpo <- ggtree(rpo_raxml) %<+% annos_rpoB
ort_rpo + geom_tiplab(aes(label=Strain),align=T,size=ts)   +  geom_treescale() + 
  geom_text(aes(x=1,y=y,label=serovar),size=ts)  + geom_text(aes(x=0.9,y=y,label=Host),size=ts) + 
  geom_text2(aes(x=branch,label=bootstrap, subset = !is.na(as.numeric(bootstrap)) & as.numeric(bootstrap) >= 50),position = position_nudge(y =0.5),size=ts)

# mlst tree
mlst_annos <- read.table(file = "mlst_combinations_host_country_serotype.tsv",header=T,sep = "\t")
mlst_annos$host=sub("chicken|turkey|pheasant",'Phasanidae',mlst_annos$host,perl=T)
mlst_annos$host=sub("other",'Turkey vulture',mlst_annos$host,perl=T)
mlst_annos$host= sapply(mlst_annos$host,function(x) paste0(toupper(substring(x,1,1)),substring(x,2)) )
mlst_annos$serotype=sub("\\s+\\(.*\\)",'',mlst_annos$serotype,perl=T)
mlst_annos$serotype=gsub("/",',',mlst_annos$serotype,perl = TRUE)
#brackets are only slight cross reactions
#mlst_annos$serotype=gsub("[()]",'',mlst_annos$serotype,perl = TRUE)
mlst_annos$serotype=sub(",\\s+",",",mlst_annos$serotype,perl=T)
mlst_annos$serotype=gsub("\\s+",',',mlst_annos$serotype,perl = TRUE)
#mlst_annos$serotype=sub(",\\s+","/",mlst_annos$serotype,perl=T)
mlst_sh <- aggregate(mlst_annos[,c("serotype","host")] ,by=list(mlst_annos$ST), function(x) (c(as.character(x))))
mlst_sh$serotype <- sapply(mlst_sh$serotype,function(x) paste(unique(sort(unlist(strsplit(x,",",fixed=T)))),collapse = ","))
mlst_sh$host <- sapply(mlst_sh$host,function(x) paste(unique(sort(x)),collapse = ","))
colnames(mlst_sh)[1]<-"ST"
mlst_sh$ST<-as.character(mlst_sh$ST)
mlst_anno <- rbind(mlst_sh,data.frame(ST=c("LMG18856","LMG18861","LMG19032"),host=c("Phasanidae"),serotype=c("F","K","M")))
mlst_anno <- rbind(mlst_sh,data.frame(ST=c("LMG18856","LMG18861","LMG19032"),host=c("Turkey"),serotype=c("F","K","M")))
#mlst_tree <- read.nexus(file="mlst_bootstraps_rotated_rooted")
#mlst_tree <- treeio::read.nexus(file="mlst_bootstraps_rotated_rooted_alt")
mlst_tree <- treeio::read.beast(file="mlst_bootstraps_rotated_rooted_alt")
x <- as_tibble(mlst_tree)
x[,'!rotate'] <- NULL
x$label[x$label == "4/LMG18856"] = "LMG18856"
mlst_tree <- as.treedata(x)
#mlst_tree <- read.raxml(file="~/vetlinux01_temp2/Lukas/Merima/MLST/RAxML_bipartitionsBranchLabels.no4_1000bs_raxml")
MRCA(mlst_tree,"21","28")
#mlst_tree <-root(mlst_tree,node=63)
#mlst_tree <- new.phylo
ggtree(mlst_tree)  + geom_tiplab(aes(label=label),align=F,size=ts) + geom_text(aes(x=branch,label=node))
#mlst_tree@phylo <- root(mlst_tree@phylo, which(mlst_tree@phylo$tip.label == "LMG19032"),edgelabel = TRUE)

mlst_anno$species<-as.character(mlst_anno$ST)
#mlst_anno$serotype2 = gsub("/",",",mlst_anno$serotype)
#mlst_anno$serotype2=sapply(mlst_anno$serotype2, function(x) paste(unique(unlist(strsplit(x,",",fixed=T))),collapse=",") )

mlst_anno$species[grep('^\\d',mlst_anno$species,perl=T)] = paste0("ST",mlst_anno$species[grep('^\\d',mlst_anno$species,perl=T)])
mlst_anno$species[mlst_anno$species=="LMG18856"] <- "LMG18856, ST4"
#MRCA(mlst_tree,"21","28")

#mlst_tree <-root(mlst_tree,node=63)
#mlst_tree <-root(mlst_tree,node=39)
ort_mlst <- ggtree(mlst_tree) %<+% mlst_anno
ort_mlst + geom_tiplab(aes(label=species),align=T,size=ts)  +  geom_treescale() +   geom_text(aes(x=branch,label=node))

 
  geom_text(aes(x=0.3,y=y,label=serotype),size=ts)  + geom_text(aes(x=0.4,y=y,label=host),size=ts)


#ort_mlst + geom_tiplab(aes(label=species),align=F,size=ts)  +  geom_text2(aes(label=node))
#p <- ggtree::rotate(ort_mlst,43)  %>%  ggtree::rotate(44)  %>%  ggtree::rotate(45)  %>%  ggtree::rotate(46) %>%  ggtree::rotate(63) %>%  ggtree::rotate(70)  %>%  ggtree::rotate(71) %>%  ggtree::rotate(67)  %>%  ggtree::rotate(59)  %>%  ggtree::rotate(57) %>% ggtree::rotate(38) 
  p <- ggtree::rotate(ort_mlst,65) %>% ggtree::flip(59,40) %>%  ggtree::rotate(40) %>%  ggtree::flip(42,45) %>%  ggtree::rotate(45) %>%  
    ggtree::rotate(55) %>%  ggtree::rotate(53)
  ost=0.04
  p  + geom_tiplab(aes(label=species),align=T,size=ts)  +  geom_treescale()  +  #geom_text(aes(x=branch,label=node)) + 
    geom_cladelabel(node=40, label="B", align=TRUE,  offset = ost) + geom_cladelabel(node=59, label="A", align=TRUE,  offset = ost) + geom_cladelabel(node=65, label="C", align=TRUE,  offset = ost) + 
    geom_cladelabel(node=45, label="Ba", align=TRUE,  offset = ost+0.01) + geom_cladelabel(node=42, label="Bb", align=TRUE,  offset = ost+0.01) +
    geom_text(aes(x=0.3625,y=y,label=serotype),size=ts)  + geom_text(aes(x=0.32,y=y,label=host),size=ts)  + xlim_tree(0.37) +
    geom_text2(aes(x=branch,label=bootstrap, subset = !is.na(as.numeric(bootstrap)) & as.numeric(bootstrap) >= 50),position = position_nudge(y =0.5),size=ts)
# for long names
  p <- ggtree::rotate(ort_mlst,65) %>% ggtree::flip(59,40) %>%  ggtree::rotate(40) %>%  ggtree::flip(42,45) %>%  ggtree::rotate(45) %>%  
    ggtree::rotate(55) %>%  ggtree::rotate(53)
  ost=0.041
  p  + geom_tiplab(aes(label=species),align=T,size=ts)  +  geom_treescale()  +  #geom_text(aes(x=branch,label=node)) + 
    geom_cladelabel(node=40, label="B", align=TRUE,  offset = ost) + geom_cladelabel(node=59, label="A", align=TRUE,  offset = ost) + geom_cladelabel(node=65, label="C", align=TRUE,  offset = ost) + 
    geom_cladelabel(node=45, label="Ba", align=TRUE,  offset = ost+0.01) + geom_cladelabel(node=42, label="Bb", align=TRUE,  offset = ost+0.01) +
    geom_text(aes(x=0.3675,y=y,label=serotype),size=ts)  + geom_text(aes(x=0.32,y=y,label=host),size=ts)  + xlim_tree(0.374) +
    geom_text2(aes(x=branch,label=bootstrap, subset = !is.na(as.numeric(bootstrap)) & as.numeric(bootstrap) >= 50),position = position_nudge(y =0.5),size=ts)
  
  
ubcg_tree <- treeio::read.beast("~/vetlinux01_temp2/Lukas/Merima/UBCG/output/all_strains_raxml/all_strains_tree")  
x <- as_tibble(ubcg_tree)  
ubcg_names=as.data.frame(matrix(c("DSM_15997", "DSM-15997",
"H06-030791", "H06-030791",
"LMG18856", "LMG-18856",
"LMG18861", "LMG-18861",
"LMG19032", "LMG-19032",
"ORT-UMN88", "ORT-UMN 88",
"Riemerella_anatipestifer", "R. anatipestifer",
'UNSC01.1',"Or. sp. OH-22767",
'UNSD01.1',"Or. sp. OH-2280"),ncol = 2,byrow = T))
colnames(ubcg_names)=c("label","species")

ubcg_tree <- treeio::read.beast("~/vetlinux01_temp2/Lukas/Merima/UBCG/output/all_strains_raxml/all_strains_tree")
new.phylo <- treeio::root(ubcg_tree,"Riemerella_anatipestifer",resolve.root = T,edgelabel = T)
new.phylo <- phytools::reroot(as.phylo(ubcg_tree),6,position=0.15)
new.phylo <- treeio::drop.tip(new.phylo, "Riemerella_anatipestifer" )

xgsi = x[,c("branch.length","gsi")]
xgsi = merge(as_tibble(new.phylo),xgsi)

#as_tibble(ubcg_tree) %<+% x %<+% ubcg_names
ggtree(ubcg_tree) %<+% x %<+% ubcg_names + geom_tiplab(aes(label=species),align=T,size=ts,hjust=-0.5) +  geom_treescale() +  xlim_tree(1.25) +
  geom_text2(aes(x=branch,label=gsi),position = position_nudge(y =0.1),size=ts) + geom_text2(aes(label=node)) 
ggtree(new.phylo) %<+% xgsi %<+% ubcg_names + geom_tiplab(aes(label=species),align=T,size=ts,hjust=-0.01) +  geom_treescale()  + 
  xlim_tree(0.55) +
  geom_text2(aes(x=branch,label=gsi, subset = !(node == 10)),position = position_nudge(y =0.15),size=ts) +
  geom_text2(aes(x=branch,label=gsi, subset = (node == 10)),position = position_nudge(x =0.01),size=ts) # + geom_text2(aes(label=node)) 

ggtree(tree_subset(ubcg_tree,node ="LMG18856" ,  levels_back = 5)) %<+% x %<+% ubcg_names + geom_tiplab(aes(label=species),align=T,size=ts,hjust=-0.1) +  geom_treescale() +  xlim_tree(0.7) +
  geom_text2(aes(x=branch,label=gsi),position = position_nudge(y =0.1),size=ts) #+ geom_text2(aes(label=node)) 

#+ geom_text2(aes(x=branch,label=bootstrap, subset = !is.na(as.numeric(bootstrap)) & as.numeric(bootstrap) >= 101),position = position_nudge(y =0.5),size=ts) 

ort_16S_t +   geom_treescale() + geom_tiplab(aes(label=Strain,align=T)) +  geom_facet(panel="Serovar",data=annos_16S,geom=geom_text,mapping=aes(x=0,label=serovar,hjust="middle"))

geom_facet(panel="ST",data=annos_16S ,geom=geom_text,aes(label=serovar))

ort_16S_t + geom_tiplab(aes(label=paste0(Strain, Host, ST, serovar, sep="\t")),align=T) + xlim_tree(0.4) +  geom_treescale()

ort_16S_t + geom_text2(aes(subset=!isTip, label=bootstrap), hjust=-.3)


ggtree(sel_ort_raxml) + geom_text2(aes(subset=!isTip, label=bootstrap), hjust=-.3) + geom_tiplab()
ggtree(sel_ort_raxml) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()

of_og_overlap=read.table("/Volumes/vetlinux01/Merima/Orthofinder/all_ort/OrthoFinder/Results_Feb20/Comparative_Genomics_Statistics/Orthogroups_SpeciesOverlaps.tsv",header=T)

ogs_per_species=read.table("/Volumes/vetlinux01/Merima/Orthofinder/all_ort/OrthoFinder/Results_Feb20/Orthogroups/Orthogroups_n0_prefix.tsv",header=T,sep = "\t")

ogs_per_species$OG[which(ogs_per_species$LMG18856 != "" & ogs_per_species$LMG18861 != "" & ogs_per_species$LMG19032 != "" & ogs_per_species$ORT_DSM_15997 == "" & ogs_per_species$ORT_H06_030791 == "" & ogs_per_species$ORT_UMN_88 == "")]
lmg_spec_ogs=which(ogs_per_species$LMG18856 != "" & ogs_per_species$LMG18861 != "" & ogs_per_species$LMG19032 != "" & ogs_per_species$ORT_DSM_15997 == "" & ogs_per_species$ORT_H06_030791 == "" & ogs_per_species$ORT_UMN_88 == "")

a = ogs_per_species$LMG18861[lmg_spec_ogs]
b =  unlist(strsplit(as.character(a), ', ', fixed=T))
egg_LMG18861$COG[egg_LMG18861$query_name %in% b]
egg_LMG18861 = read.table("/Volumes/vetlinux01/Merima//eggnog/HMM/LMG18861.emapper.annotations_no_pref",sep = "\t")
colnames(egg_LMG18861) = c("query_name", "seed_eggNOG_ortholog","seed_ortholog_evalue","seed_ortholog_score","predicted_gene_name", "GO_terms","KEGG_KOs","BiGG_reactions","Annotation_tax_scope","OGs","bestOG|evalue|score","COG","eggNOG_annot")
table(egg_LMG18861$COG)






# orthoani 
library(pheatmap)
library(RColorBrewer)


anis <- read.csv(file="~/merima_rem/ORT/orthoani_orts.csv",header=TRUE,sep=",")
anis$X<-NULL
for(i in ncol(anis):2){
  anis[,i] <- anis[,i-1]
}
colnames(anis) <- sub('_rseed','',colnames(anis))
rownames(anis) <- colnames(anis)
anis<-as.matrix(anis)
diag(anis) <- 100
anis[lower.tri(anis)]
anis[lower.tri(anis)]=t(anis)[lower.tri(anis)]
hc = hclust(as.dist(1 - anis/100),method="average")
pheatmap(anis, cluster_col = hc, cluster_row = hc, display_numbers = TRUE,treeheight_row=200,treeheight_col = FALSE)
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(88, 95, 100), c("red", "green","blue"))
Heatmap(anis,cluster_rows = hc,cluster_columns = hc,show_column_dend = FALSE,cell_fun = function(j, i,x, y, width, height, fill){grid.text(sprintf("%.2f", anis[i,j] ),x, y)}, col = col_fun)

anis_emen <- read.csv(file="~/merima_rem/ORT/orthoani_emen_orts.csv",header=TRUE,sep=",")
anis_emen$X<-NULL
#for(i in ncol(anis_emen):2){
##  anis_emen[,i] <- anis_emen[,i-1]
#}
colnames(anis_emen) <- sub('_rseed','',colnames(anis_emen))
colnames(anis_emen) <- sub('Elizab.*septica','E. men.',colnames(anis_emen))
rownames(anis_emen) <- colnames(anis_emen)
anis_emen<-as.matrix(anis_emen)
diag(anis_emen) <- 100
anis_emen[lower.tri(anis_emen)]
anis_emen[lower.tri(anis_emen)]=t(anis_emen)[lower.tri(anis_emen)]
hc = hclust(as.dist(1 - anis_emen/100),method="average")
col_map=c(colorRampPalette(colors=c("red","orange"))((90-60)/0.1),colorRampPalette(colors=c("orange","lightblue"))((95-90)/0.1),colorRampPalette(colors=c("lightblue","blue"))((100-95)/0.1))
break_map=seq(60,100,by=0.1)
pheatmap(anis_emen, cluster_col = hc, cluster_row = hc, display_numbers = TRUE,treeheight_row=100,treeheight_col = FALSE,color = col_map,breaks = break_map)
col_map

library(seqinr)
help("seqinr")
sel_orts_mask = read.alignment("/Volumes/vetlinux01/Merima/16S-trees/all_16Srna_nogap.mafft.mfa",format = "fasta" )
sel_orts_mask_dist = dist.alignment(sel_orts_mask,matrix="identity") #,diag=TRUE)
id_mat = 1-as.matrix(sel_orts_mask_dist^2)
library(RColorBrewer)

library(corrplot)
vignette("corrplot-intro")

corrplot(id_mat,is.corr = FALSE,cl.lim=c(0.80,1.0),type="upper",diag=FALSE,tl.col="black", tl.srt=45,col=brewer.pal(n=10, name="RdBu"))

all_anis <- read.csv(file="~/Data/Merima/orthoani_ort_oh_ran_matrix_alldists.csv",header=TRUE,sep=",",skip=2)

all_anis$X<-NULL
for(i in ncol(all_anis):2){
  all_anis[,i] <- all_anis[,i-1]
}

anis <- list()
anis[["orig"]] <- all_anis[grep('(Original ANI)',rownames(all_anis),fixed=TRUE),]
anis[["ortho"]] <- all_anis[grep('(OrthoANI)',rownames(all_anis),fixed=TRUE),]
anis[["ggdc"]] <- all_anis[grep('(GGDC)',rownames(all_anis),fixed=TRUE),]
for(i in names(anis)){
  colnames(anis[[i]]) <- sub('_rseed','',colnames(anis[[i]]))
  rownames(anis[[i]]) <- colnames(anis[[i]])
  anis[[i]]<-as.matrix(anis[[i]])
  diag(anis[[i]]) <- 100
  anis[[i]][lower.tri(anis[[i]])]
  anis[[i]][lower.tri(anis[[i]])]=t(anis[[i]])[lower.tri(anis[[i]])]
  hc = hclust(as.dist(1 - anis[[i]]/100),method="average")
} 

farp_anis <- read.csv(file="~/Data/Merima/orthoani_ort_farper.csv",header=TRUE,sep=",",skip=3, row.names = NULL)
farp_anis <- farp_anis[,1:3]
all_labs=unique(c(as.character(farp_anis[,1]),as.character(farp_anis[,2])))
farp_mat <- as.data.frame(matrix(100,nrow=length(all_labs),ncol=length(all_labs)))
rownames(farp_mat)=all_labs
colnames(farp_mat)=all_labs
for(i in 1:nrow(farp_anis)) {
  #farp_mat[as.character(farp_anis[i,1]),as.character(farp_anis[i,1])] <- 100
  farp_mat[as.character(farp_anis[i,1]),as.character(farp_anis[i,2])] <- farp_anis[i,3]
  farp_mat[as.character(farp_anis[i,2]),as.character(farp_anis[i,1])] <- farp_anis[i,3]
}

new_labs <- c("LMG-18856","LMG-18861","LMG-19032","FARPER-174b","O. sp. OH-22767","O. sp. OH-2280", "ORT-UMN 88" ,"H06-030791","DSM-15997")
hc = hclust(as.dist(1 - farp_mat/100),method="average")
col_map=c(colorRampPalette(colors=c("red","orange"))((90-60)/0.1),colorRampPalette(colors=c("orange","lightblue"))((95-90)/0.1),colorRampPalette(colors=c("lightblue","blue"))((100-95)/0.1))
break_map=seq(60,100,by=0.1)
pheatmap(farp_mat, cluster_col = hc, cluster_row = hc, display_numbers = TRUE,treeheight_row=100, treeheight_col = FALSE,
         color = col_map,breaks = break_map,labels_row = new_labs,labels_col = new_labs,number_color="black",fontsize_number=10)

anis_i<- anis[["orig"]]
anis_i<- anis[["ortho"]]
new_labs <- c("DSM-15997" ,"H06-030791","ORT-UMN 88","R. anatipestifer","O. sp. OH-22767","O. sp. OH-2280","LMG-18856","LMG-18861","LMG-19032")
hc = hclust(as.dist(1 - anis_i/100),method="average")
col_map=c(colorRampPalette(colors=c("red","orange"))((90-60)/0.1),colorRampPalette(colors=c("orange","lightblue"))((95-90)/0.1),colorRampPalette(colors=c("lightblue","blue"))((100-95)/0.1))
break_map=seq(60,100,by=0.1)
pheatmap(anis_i, cluster_col = hc, cluster_row = hc, display_numbers = TRUE,treeheight_row=100,treeheight_col = FALSE,
         color = col_map,breaks = break_map,labels_row = new_labs,labels_col = new_labs,number_color="black",fontsize_number=10)
drop<-grep("ana",colnames(anis_i))
anis_i<-anis_i[- drop, - drop]
hc = hclust(as.dist(1 - anis_i/100),method="average")
pheatmap(anis_i, cluster_col = hc, cluster_row = hc, display_numbers = TRUE,treeheight_row=100,treeheight_col = FALSE,
         color = col_map,breaks = break_map,labels_row = new_labs[-drop],labels_col = new_labs[-drop],number_color="black",fontsize_number=10)

anis_i<- anis[["ggdc"]]

ggdc_tab=read.table(file="ggdc.txt",sep = "\t", header = T)

reshape(ggdc_tab[,1:3],idvar="Query.genome",timevar = "Reference.genome",v.names = "DDH",direction = "wide")
#reshape2
ggdc<-acast(ggdc_tab,Query.genome~Reference.genome,value.var = "DDH")
ggdc_int<-acast(ggdc_tab,Query.genome~Reference.genome,value.var = "Model.C.I.")
col_map=c(colorRampPalette(colors=c("red","orange"))((70-10)/0.1),colorRampPalette(colors=c("orange","lightblue"))((79-70)/0.1),colorRampPalette(colors=c("lightblue","blue"))((100-79)/0.1))
break_map=seq(10,100,by=0.1)
new_labs <- c("LMG-18856","LMG-18861","LMG-19032","DSM-15997" ,"H06-030791","ORT-UMN 88","O. sp. OH-22767","O. sp. OH-2280")
hc = hclust(as.dist(1 - ggdc/100),method="average")
pheatmap(ggdc, cluster_col = hc, cluster_row = hc, display_numbers = TRUE,treeheight_row=100,treeheight_col = FALSE,color = col_map,breaks = break_map,labels_row = new_labs,labels_col = new_labs)
# with intervals
pheatmap(ggdc, cluster_col = hc, cluster_row = hc, treeheight_row=100,treeheight_col = FALSE,color = col_map,breaks = break_map,labels_row = new_labs,labels_col = new_labs,display_numbers=ggdc_int)


         
