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
library(ComplexHeatmap)
library(circlize)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/Cis_Regulatory_Elements/")
enrichment <-read.table("22_3_merge_cCREs_enrichment.txt.gz",header = T,sep = "\t") %>% as.data.frame()
# sample_id <- read_excel("chromHmm_1581.xlsx", sheet = "IDs")%>% as.data.frame()

enrichment$EnrichFDR <- p.adjust(enrichment$EnrichPValue,method="fdr")
enrichment$DepleteFDR <- p.adjust(enrichment$DepletePValue,method="fdr")
enrichment$SignificantEnrichment <-NA
for(i in 1:nrow(enrichment)){
    if(enrichment[i,"EnrichFDR"]<0.05 |enrichment[i,"DepleteFDR"]<0.05 ){
        enrichment[i,"SignificantEnrichment"] <-enrichment[i,"Enrichment"]
    }
}

dat <- dcast(enrichment[,c("cluster","type","SignificantEnrichment")],cluster~type)
rownames(dat) <-dat$cluster
dat <-dat[,-1]

p1 <- Heatmap(as.matrix(dat), name = "Overlap enrichment", row_names_side = "left",
    # column_names_side = "top",
    column_names_rot = 45,
    cluster_rows = FALSE, cluster_columns = FALSE,
    width = ncol(dat)*unit(10, "mm"),height = nrow(dat)*unit(10, "mm"),
    row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),
    col = colorRamp2(c(0, 1, 5), c("#0000FF", "white", "#FF0000")) ,
    # Legend(labels_gp = gpar(fontsize = 5, fontface = "bold")))
    heatmap_legend_param = list(title_gp = gpar(fontsize = 6.5)))
    # col = colorRamp2(c(0, 1, 4), c("#1E295C", "white", "#C02025")))
pdf("22_4_heatmap_cCREs_enrichement.pdf",width=3.5,height=3.5)
print(p1)
dev.off()

# color2 = colorRampPalette(c('#c6dbef','#6baed6','#2171b5'))(50)

# library(pheatmap)
# pdf("22_2_cCRE_propotion_heatmap.pdf",width=2.7,height=2.8)
# pheatmap(fdat[,2:6],cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = TRUE,show_rownames = T,show_colnames = T)
# dev.off()
