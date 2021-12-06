#---------------------
library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(Seurat)
library(circlize)
# cancers <-c("BRCA","KIRP")

setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/prediction/output/")
trans <- fread("05_trans_predict.txt.gz",header = T,sep = "\t") %>% as.data.frame()
trans1 <-filter(trans,gene_chr %in% unique(trans$h_chr))

cis <- fread("05_cis_predict.txt.gz",header = T,sep = "\t") %>% as.data.frame()


hotspot <-trans1%>%dplyr::select(h_chr,h_start,h_end)
gene <- trans1%>%dplyr::select(gene_chr,gene_start,gene_end)
#_-------------
circos.clear()
set.seed(80)

pdf("./figure/06_trans_circos_link_prediction.pdf")
par(mar = c(1, 1, 1, 1))
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species= "hg19",chromosome.index = paste0("chr", 1:22))
circos.genomicLink(hotspot, gene, col = sample(nrow(hotspot), nrow(gene), replace = TRUE))
# circos.genomicLink(trans_qtl, trans_egene)
dev.off()
circos.clear()
#-------------------
    #--------BRCA
# cancer =