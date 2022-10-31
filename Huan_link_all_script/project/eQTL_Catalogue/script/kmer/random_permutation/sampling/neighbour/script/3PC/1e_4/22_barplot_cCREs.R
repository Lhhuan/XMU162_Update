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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/Cis_Regulatory_Elements/")
org <-read.table("hotspot_cluster_cCREs.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(org)<-c("CHR","start","end","cluster","cre_chr","cre_start","cre_end","name","score","strand","thickStart","thickEnd","reserved","ccre","encodeLabel","zScore","ucscLabel","accessionLabel","description")
org$hotspot <-paste(org$CHR,org$start,org$end,sep=":")
org$cluster <-as.factor(org$cluster)
org <- org[,c("hotspot","cluster","ucscLabel","name")]
org1 <- unique(org[,c("hotspot","cluster","ucscLabel")])

dat <-org1%>%group_by(ucscLabel,cluster)%>%summarise(ucscLabel_count=n())%>%data.frame()

all_hotspot <- read.table("../GWAS/17_whole_genome_leiden_pca3_k50_resolution1e-04_sorted.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(all_hotspot) <-c("CHR","start","end","cluster")
all_hotspot$cluster <-as.factor(all_hotspot$cluster)
all_hotspot$hotspot <- paste(all_hotspot$CHR,all_hotspot$start,all_hotspot$end,sep=":")
all_hotspot <- all_hotspot[,c("hotspot","cluster")]
cluster_n <- all_hotspot %>%group_by(cluster)%>%summarise(hotspot_count=n())%>%data.frame()


dat <- left_join(dat,cluster_n,by=c("cluster"))
dat$proporation <-dat$ucscLabel_count/dat$hotspot_count

p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5))

plist<-list()
for(i in unique(dat$ucscLabel)){
    fdatcount <-filter(dat,ucscLabel==i)
    fdatcount$cluster <-factor(fdatcount$cluster,levels=c(4,2,1,5,6,3))
   plist[[i]] <- ggplot(fdatcount, aes(x=as.factor(cluster), y=proporation,group=cluster,fill = as.factor(cluster))) + 
    # geom_boxplot(outlier.colour = NA)+
    geom_bar(stat = 'identity', width=0.8)+
    scale_fill_manual(values=c("#C22324","#FF7F0E","#1E77B4","#9567BD","#8C554B","#2CA02C"))+
    theme_bw()+
    labs(x="Cluster",y="Fraction of hotspots")+
    ggtitle(i)+
    theme(legend.position ="none")+
    p_theme
    print(i)
    pdf(paste0("12_barplot_of_",i,".pdf"),height=4.7,width=4.5)
    print(plist[[i]])
    dev.off()
}

pdf("22_barplot_of_cCRE.pdf",width=12, height=2.5)
gridExtra::marrangeGrob(plist,ncol=5,nrow=1,top="Candidate Cis-Regulatory Elements")
dev.off()





names(pca_list)=c("rbfdot","tanhdot","anovadot","vanilladot","polydot")
# save(pca_list,file = "08_kpca.Rdata")
plot <-function(i=NULL){
    aa=as.data.frame(pca_list[[i]]@rotated)
    colnames(aa)[1:2] <- c("PC1","PC2")
    p =ggplot(aa) + geom_point(aes(x = PC1, y = PC2), size = 0.1, alpha = 0.3)+ggtitle(names(pca_list)[i])+theme_bw() + p_theme
    return(p)
}

plist <-lapply(c(1:5),plot)
pdf(paste0("./egene_cluster/13_1_kpca_sampling_result.pdf"),width=8.3, height=4)
CombinePlots(plist,ncol=4,nrow=2)
dev.off()


p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5))

dat1 <-filter(all_hotspot,cluster %in%c(1:6))
dat1$cluster <-factor(dat1$cluster,levels=c(4,2,1,5,6,3))
p1 <- ggplot(dat1, aes(x=as.factor(cluster), y=log(adjust_trait_count),group=cluster,fill = as.factor(cluster))) + 
# geom_boxplot(outlier.colour = NA)+
geom_boxplot()+
# scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
# scale_fill_manual(values=c("#8C554B","#C22324","#1E77B4","#FF7F0E","#9567BD","#2CA02C"))+
scale_fill_manual(values=c("#C22324","#FF7F0E","#1E77B4","#9567BD","#8C554B","#2CA02C"))+
theme_bw()+
labs(x="Cluster",y="Log(number of snp-trait pairs per kb)")+
ggtitle("SNP-trait pairs")+
theme(legend.position ="none")+
p_theme

pdf("18_boxpolt_whole_genome_gwas_trait.pdf",height=4.7,width=4.5)
print(p1)
dev.off()



p1 <- ggplot(dat1, aes(x=as.factor(cluster), y=log(adjust_trait_count+0.1),group=cluster,fill = as.factor(cluster))) + 
# geom_boxplot(outlier.colour = NA)+
geom_boxplot()+
# scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
# scale_fill_manual(values=c("#8C554B","#C22324","#1E77B4","#FF7F0E","#9567BD","#2CA02C"))+
scale_fill_manual(values=c("#C22324","#FF7F0E","#1E77B4","#9567BD","#8C554B","#2CA02C"))+
theme_bw()+
labs(x="Cluster",y="Log(number of snp-trait pairs per kb)")+
ggtitle("SNP-trait pairs")+
theme(legend.position ="none")+
p_theme

pdf("18_boxpolt_whole_genome_gwas_trait_1E-8.pdf",height=4.7,width=4.5)
print(p1)
dev.off()



p <- ggplot(dat1, aes(x=log(adjust_trait_count),color=cluster)) + 
geom_density(size=0.8)+
scale_color_d3("category20") + 
theme_bw()+
labs(x="Log(Number of snp-trait pairs per kb)",color="Cluster") +
ggtitle("SNP-trait pairs")+
p_theme

pdf("18_densityplot_whole_genome_gwas_trait.pdf",height=4.5,width=4.9)
print(p)
dev.off()

all_dat <-dat1%>%group_by(cluster)%>%summarise(cluster_count=n())%>%data.frame()
aa <-filter(dat1,adjust_trait_count>0)
aa_dat <-aa%>%group_by(cluster)%>%summarise(GWAS_anno_count=n())%>%data.frame()
fdatcount <- left_join(all_dat,aa_dat,by="cluster")
fdatcount$propotion <- fdatcount$GWAS_anno_count /fdatcount$cluster_count
#===================
fdatcount$cluster <-factor(fdatcount$cluster,levels=c(4,2,1,5,6,3))
p1 <- ggplot(fdatcount, aes(x=as.factor(cluster), y=propotion,group=cluster,fill = as.factor(cluster))) + 
# geom_boxplot(outlier.colour = NA)+
geom_bar(stat = 'identity', width=0.8)+
# scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
# scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B"))+
# scale_fill_manual(values=c("#8C554B","#C22324","#1E77B4","#FF7F0E","#9567BD","#2CA02C"))+
scale_fill_manual(values=c("#C22324","#FF7F0E","#1E77B4","#9567BD","#8C554B","#2CA02C"))+
theme_bw()+
labs(x="Cluster",y="Fraction of hotspots")+
ggtitle("GWAS")+
theme(legend.position ="none")+
p_theme

pdf("18_barpolt_whole_genome_GWAS.pdf",height=4.7,width=4.5)
print(p1)
dev.off()