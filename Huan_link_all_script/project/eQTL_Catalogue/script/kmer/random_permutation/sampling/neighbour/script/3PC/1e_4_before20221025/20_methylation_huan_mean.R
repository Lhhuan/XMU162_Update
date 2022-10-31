library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(Seurat)
library(reshape2)
library(R.utils)
library(ggsci)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/methylation/")
methylation <-read.table("/home/huanhuan/project/eQTL_Catalogue/output/annotation/extend_18snp/Methylation_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(methylation)[1:7]<-c("h_chr","h_start","h_end","m_chr","m_start","m_end","tissue_score")
A <- strsplit(methylation$tissue_score,split=";")
tmp <-do.call(rbind,A)%>%as.data.frame()
tmp <-tmp[,1:2]
colnames(tmp) <-c("tissue","score")
dat <- bind_cols(methylation,tmp)
dat$hotspot <- paste0(dat$h_chr,":",dat$h_start,"-",dat$h_end)
dat$methylation_site<-paste(dat$m_chr,dat$m_start,dat$m_end,sep=":")
dat$h_length <- dat$h_end - dat$h_start
dat <- dat[,c("hotspot","methylation_site","tissue","h_length","score")]
dat$score <-as.numeric(dat$score)
Mscore <- dat%>%group_by(hotspot,tissue,h_length)%>%summarise(total_score=sum(score))%>%data.frame()
Mscore$adjust_score <- Mscore$total_score/ Mscore$h_length

Cluster <-read.table("../../11_whole_genome_leiden_pca3_k50_resolution1e-04.txt",header = T,sep = "\t") %>% as.data.frame()
Mscore1 <-left_join(Cluster[,c("cluster","hotspot")],Mscore,by="hotspot")
Mscore1$adjust_score[is.na(Mscore1$adjust_score)] <-0

p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5))


Mscore1$cluster <-factor(Mscore1$cluster,levels=c(6,4,1,2,5,3))
p1 <- ggplot(Mscore1, aes(x=as.factor(cluster), y=log(adjust_score),group=cluster,fill = as.factor(cluster))) + 
# geom_boxplot(outlier.colour = NA)+
geom_boxplot()+
# scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
scale_fill_manual(values=c("#8C554B","#C22324","#1E77B4","#FF7F0E","#9567BD","#2CA02C"))+
theme_bw()+
labs(x="Cluster",y="Log(score per bp)")+
ggtitle("DNA methylation")+
theme(legend.position ="none")+
p_theme

pdf("20_boxpolt_whole_genome_dna_methylation_60_cell.pdf",height=4.7,width=4.5)
print(p1)
dev.off()

p1 <- ggplot(Mscore1, aes(x=as.factor(cluster), y=log(adjust_score+0.00000001),group=cluster,fill = as.factor(cluster))) + 
geom_boxplot()+
# scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
scale_fill_manual(values=c("#8C554B","#C22324","#1E77B4","#FF7F0E","#9567BD","#2CA02C"))+
theme_bw()+
labs(x="Cluster",y="Log(score per bp)")+
ggtitle("DNA methylation(1e-8)")+
theme(legend.position ="none")+
p_theme

pdf("20_boxpolt_whole_genome_dna_methylation_60_cell_1e-08.pdf",height=4.7,width=4.5)
print(p1)
dev.off()
#---------------------------------------------------mean
Mscore2 <- Mscore1%>%group_by(hotspot,cluster)%>%summarise(mean_score=mean(adjust_score))%>%data.frame()

Mscore2$cluster <-factor(Mscore2$cluster,levels=c(6,4,1,2,5,3))
# p1 <- ggplot(Mscore2, aes(x=as.factor(cluster), y=log(mean_score),group=cluster,fill = as.factor(cluster))) + 
p1 <- ggplot(Mscore2, aes(x=as.factor(cluster), y=log(mean_score),group=cluster,fill = as.factor(cluster))) + 
# geom_boxplot(outlier.colour = NA)+
geom_boxplot()+
# scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
scale_fill_manual(values=c("#8C554B","#C22324","#1E77B4","#FF7F0E","#9567BD","#2CA02C"))+
theme_bw()+
labs(x="Cluster",y="Log(mean score per bp across tissue)")+
ggtitle("DNA methylation")+
theme(legend.position ="none")+
p_theme

pdf("20_boxpolt_whole_genome_dna_methylation_LOG_mean_score.pdf",height=4.7,width=4.5)
print(p1)
dev.off()
#----------------------------------
p1 <- ggplot(Mscore2, aes(x=as.factor(cluster), y=log(mean_score+0.00000001),group=cluster,fill = as.factor(cluster))) + 
# geom_boxplot(outlier.colour = NA)+
geom_boxplot()+
# scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
scale_fill_manual(values=c("#8C554B","#C22324","#1E77B4","#FF7F0E","#9567BD","#2CA02C"))+
theme_bw()+
labs(x="Cluster",y="Log(mean score per bp across tissue)")+
ggtitle("DNA methylation(1e-8)")+
theme(legend.position ="none")+
p_theme

pdf("20_boxpolt_whole_genome_dna_methylation_LOG_mean_score_1e-08.pdf",height=4.7,width=4.5)
print(p1)
dev.off()






Mscore3<-filter(Mscore2,mean_score!=0)
all_dat <-Cluster%>%group_by(cluster)%>%summarise(cluster_count=n())%>%data.frame()
all_dat$cluster <-as.factor(all_dat$cluster)
methy_dat <- Mscore3%>%group_by(cluster)%>%summarise(methy_cluster_count=n())%>%data.frame()

fdat <-left_join(methy_dat,all_dat,by="cluster")
fdat$proportion <- fdat$methy_cluster_count/fdat$cluster_count

fdat$cluster <-factor(fdat$cluster,levels=c(6,4,1,2,5,3))
p1 <- ggplot(fdat, aes(x=as.factor(cluster), y=proportion,group=cluster,fill = as.factor(cluster))) + 
# geom_boxplot(outlier.colour = NA)+
geom_bar(stat = 'identity', width=0.8)+
# scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
# scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B"))+
scale_fill_manual(values=c("#8C554B","#C22324","#1E77B4","#FF7F0E","#9567BD","#2CA02C"))+
theme_bw()+
labs(x="Cluster",y="Fraction of hotspots")+
ggtitle("DNA methylation")+
theme(legend.position ="none")+
p_theme

pdf("20_barpolt_whole_genome_DNA_methylation_PRO.pdf",height=4.7,width=4.5)
print(p1)
dev.off()


