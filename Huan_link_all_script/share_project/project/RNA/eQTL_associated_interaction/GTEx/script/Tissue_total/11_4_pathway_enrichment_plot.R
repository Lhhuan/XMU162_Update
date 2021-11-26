library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(patchwork)
library(Seurat)
# library(clusterProfiler)
# library(org.Hs.eg.db)
#-------------------------------------------------------
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_total/gene/david_analyse/")
#-------------normal
org <- read.table("11_3_49_union_david.txt",header = T,sep = "\t") %>% as.data.frame()
colnames(org)[4]="gene_ratio"
org$p.adjust <-p.adjust(org$PValue,method="fdr")

h_sig <-filter(org,p.adjust<0.05)%>%select("Category","Term","gene_ratio","Count","PValue","FDR","p.adjust")
# sig <-filter(org,FDR<=0.05)%>%select("Category","Term","gene_ratio","Count","PValue","FDR","p.adjust")
h_sig$Term <-gsub("^.*:","",h_sig$Term)
h_sig$Term <-gsub("^.*~","",h_sig$Term)
unique(h_sig$Term)
unique_Category <-unique(h_sig$Category)


# for(class in unique_Category){
#     print(class)
all <-function(class){
    a <- filter(h_sig,Category==class)   
    class <-gsub("_"," ",class)
    ggplot(a,aes(x=gene_ratio,y=reorder(Term,X=gene_ratio)))+
    geom_point(aes(size=Count*2,color=-1*log10(p.adjust)))+scale_color_gradient(low="blue",high ="red")+
    labs(color=expression(-log[10](p.adjust)),size="Gene",x="gene ratio",y=class)+theme_bw()+
    theme(axis.text.y=element_text(size=8.5),plot.title=element_text(hjust = 0.5))
}

plist = lapply(unique_Category,all)

pdf("all_gene_DAVID.pdf",width=20, height=10)
CombinePlots(plist,ncol=2,nrow=2)
dev.off()

all_sg <-function(class){
    a <- filter(h_sig,Category==class)   
    class <-gsub("_"," ",class)
    ggplot(a,aes(x=gene_ratio,y=reorder(Term,X=gene_ratio)))+
    geom_point(aes(size=Count*2,color=-1*log10(p.adjust)))+scale_color_gradient(low="blue",high ="red")+
    labs(color=expression(-log[10](p.adjust)),size="Gene",x="gene ratio",y=class)+theme_bw()+
    theme(axis.text.y=element_text(size=8.5),plot.title=element_text(hjust = 0.5))
}











colnames(a)[4]="gene_ratio"
ggplot(a,aes(x=gene_ratio,y=Term))+
geom_point(aes(size=Count*2,color=-1*log10(PValue)))+scale_color_gradient(low="blue",high ="red")+
  labs(color=expression(-log[10](PValue)),size="Gene",x="gene ratio",y="Pathway name",title="Pathway enrichment")+theme_bw()

