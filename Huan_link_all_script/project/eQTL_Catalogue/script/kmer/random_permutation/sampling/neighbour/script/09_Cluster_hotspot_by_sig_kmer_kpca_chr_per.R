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
library(mclust)
library(kernlab)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure")
load("08_permutation_wilocx_overlap_sig_kmer_0_1000.Rdata")
Sorg_all <-Sorg
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) 



chrs <-paste0("chr",1:22,":")
# for(chr in chrs){
kpca_chr <-function(chr=NULL){
    Sorg <-Sorg_all[grep(chr,rownames(Sorg_all)),]
    chr_new<-chr
    chr_new <-gsub(":","",chr_new)
    # pca=prcomp(dat1)
    options(future.globals.maxSize = 50* 1024^3)
    pca1<- kpca(~.,data=Sorg, kernel="rbfdot",kpar = list(sigma = 0.1),features = 0, th = 1e-4)
    save(pca1,file=paste0("./09_chr/09_kpca_rbfdot_",chr_new,".Rdata"))
    pca=pca1@pcv
    aa=as.data.frame(pca)
    colnames(aa)[1:2] <- c("PC1","PC2")
    p =ggplot(aa) + geom_point(aes(x = PC1, y = PC2), size = 0.1, alpha = 0.3)+ggtitle(paste0(chr," rbfdot"))+theme_bw() + p_theme
    pdf(paste0("./09_chr/09_kpca_rbfdot_",chr_new,".pdf"),width=4.1, height=4)    
    print(p)
    dev.off()
    return(p)
    print(chr)
}
plist <- lapply(chrs,kpca_chr)

p =gridExtra::marrangeGrob(plist,ncol=4,nrow=2)
# ggsave("09_kpca_rbfdot_all_per_chr.png",p,dpi=120,width=12,height=12)
pdf("09_kpca_rbfdot_all_per_chr.pdf",width=9, height=5)    
print(p)
dev.off()