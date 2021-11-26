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


setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/chr1_6/6kmer/5_community/")

# aa <-c(1;5)
i=1
a <- read.table(paste0(i,"_GO_MF.txt"),header = T,sep = "\t") %>% as.data.frame()
colnames(org)[3]="gene_ratio"
# org$p.adjust <-p.adjust(org$PValue,method="fdr")

# h_sig <-filter(org,p.adjust<0.05)%>%select("Category","Term","gene_ratio","Count","PValue","FDR","p.adjust")
# sig <-filter(org,FDR<=0.05)%>%select("Category","Term","gene_ratio","Count","PValue","FDR","p.adjust")
# h_sig$Term <-gsub("^.*:","",h_sig$Term)
# h_sig$Term <-gsub("^.*~","",h_sig$Term)
# unique(h_sig$Term)
# unique_Category <-unique(h_sig$Category)


# for(class in unique_Category){
#     print(class)
all <-function(i){
    # a <- filter(h_sig,Category==class)   
    # class <-gsub("_"," ",class)
    a <- read.table(paste0(i,"_GO_MF.txt"),header = T,sep = "\t") %>% as.data.frame()
    colnames(a)[3]="gene_ratio"
    a$
    a$gene_ratio <-as.character(a$gene_ratio)
    b <- as.data.frame(apply(str_split(a$gene_ratio,"/",simplify = T),2,as.numeric))
    a$gene_ratio <- b[,1]/b[,2]  
    # a$gene_ratio <-as.numeric(as.character(a$gene_ratio))
    colnames(a)[2]="Term"
    a$Term <-str_replace(a$Term,"RNA polymerase II-specific","R II")
    p1 <-ggplot(a,aes(x=gene_ratio,y=reorder(Term,X=gene_ratio)))+
    geom_point(aes(size=Count*2,color=-1*log10(p.adjust)))+scale_color_gradient(low="blue",high ="red")+
    labs(color=expression(-log[10](p.adjust)),size="Gene",x="gene ratio",y="")+theme_bw()+ggtitle(i)+
    theme(axis.text.y=element_text(size=8.5,color="black"),plot.title=element_text(hjust = 0.5))
    pdf(paste0(i,"_GOMF_refine.pdf"))
    print(p1)
    dev.off()
}

plist = lapply(c(1:5),all)

pdf("./all_go_mf.pdf",width=15, height=15)
CombinePlots(plist,ncol=3,nrow=2)
dev.off()


    p1 <-ggplot(a,aes(x=gene_ratio,y=reorder(Term,X=gene_ratio)))+
    geom_point(aes(size=Count*2,color=-1*log10(p.adjust)))+scale_color_gradient(low="blue",high ="red")+
    labs(color=expression(-log[10](p.adjust)),size="Gene",x="gene ratio",y=i)+theme_bw()+ggtitle(i)+
    theme(axis.text.y=element_text(size=8),plot.title=element_text(hjust = 0.5))
    pdf(paste0(i,"_GOMF_refine.pdf"))
    print(p1)
    dev.off()