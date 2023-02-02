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
library(umap)
library(uwot)
library(kernlab)
library(data.table)
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/cancer/output/")

f1 <- read.table("02_cancer_cellLine_marker.txt",header = T,sep = "\t") %>% as.data.frame()
tissues <-sort(unique(f1$tmp_tissue))
p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    legend.position ="none",
    plot.title = element_text(hjust = 0.5,size=11))

plist <- list()
for (i in 1:length(tissues)){
    tissue = tissues[i]
    new_dir <- paste0("/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/cancer/cell_line/",tissue,"/GC_content")
    if(dir.exists(new_dir)){
        print(c(new_dir, "exists"))
    }else{
        dir.create(new_dir, recursive = TRUE)
    }
    setwd(new_dir)
    org <- fread("../predicted_matrix_0_mean_allFeatute.txt.gz",header =T,sep = "\t") %>% as.data.frame()
    org1 <- org[,c("hotspot","GC_content","predict_class")]
    org1$cluster <-factor(org1$predict_class,levels=c(0,1,2))
    tissue <-capitalize(tissue)
    #===========================================
    plist[[i]] <- ggplot(org1, aes(x=cluster, y=GC_content,fill =cluster)) + 
        geom_boxplot(outlier.shape = NA)+
        # scale_fill_manual(values=c("#1E77B4","#FF7F0E"))+
        scale_fill_manual(values=c("#A593E0","#84B1ED","#F6B352"))+
        theme_bw()+
        ggtitle(tissue)+
            stat_compare_means(method = "wilcox.test",comparisons =list(c("0","1"),c("1","2")),label="p.signif",method.args = list(alternative = "less"))+
        p_theme +
        labs(x="Cluster",y="GC content")+
        # ylim(0,1.21)
        scale_y_continuous(breaks=seq(0.0, 1.2, 0.25), limits=c(0, 1.2))
        # labs(y="GC content",x="")
    pdf("08_gc_content_boxplot.pdf",height=3.1,width=3)
    print(plist[[i]])
    dev.off()
    print(c(i,tissue))
}

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/normal_tissue/output/figure/")
pdf("08_ALL_tissue_boxplot_gc_content.pdf",width=12, height=9.5)
p2<-gridExtra::marrangeGrob(plist,nrow=3,ncol=4,title="GC_content")
print(p2)
dev.off()
print("finish")
