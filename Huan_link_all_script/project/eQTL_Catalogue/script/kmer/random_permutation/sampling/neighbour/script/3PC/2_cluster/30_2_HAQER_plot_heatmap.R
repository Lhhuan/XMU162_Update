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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/HAQER/")
enrichment <-read.table("30_1_merge_haqer_enrichment.txt.gz",header = T,sep = "\t") %>% as.data.frame()
# sample_id <- read_excel("chromHmm_1581.xlsx", sheet = "IDs")%>% as.data.frame()

enrichment$EnrichFDR <- p.adjust(enrichment$EnrichPValue,method="fdr")
enrichment$DepleteFDR <- p.adjust(enrichment$DepletePValue,method="fdr")
enrichment$SignificantEnrichment <-NA
for(i in 1:nrow(enrichment)){
    if(enrichment[i,"EnrichFDR"]<0.05 |enrichment[i,"DepleteFDR"]<0.05 ){
        enrichment[i,"SignificantEnrichment"] <-enrichment[i,"Enrichment"]
    }
}

write.table(enrichment,"30_2_Enrichment.txt",quote = FALSE, sep = "\t")

M_enrichment_C1 <- enrichment[,c("cluster","SignificantEnrichment")]
rownames(M_enrichment_C1)<-M_enrichment_C1[,1]
M_enrichment_C1<-M_enrichment_C1[,-1]

M_enrichment_C1 <-as.matrix(M_enrichment_C1)
colnames(M_enrichment_C1) <- "HAQER"
rownames(M_enrichment_C1) <-paste0("C",1:2)

p1 <- Heatmap(M_enrichment_C1, name = "Overlap enrichment", row_names_side = "left", column_names_side = "top",column_names_rot=0,
    cluster_rows = FALSE, cluster_columns = FALSE,
    width = ncol(M_enrichment_C1)*unit(1, "cm"),height = nrow(M_enrichment_C1)*unit(1, "cm"),
    row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 8),
    col = colorRamp2(c(0, 1, 4), c("#0000FF", "white", "#FF0000")) ,
    # Legend(labels_gp = gpar(fontsize = 5, fontface = "bold")))
    heatmap_legend_param = list(title_gp = gpar(fontsize = 6.5)))
    # col = colorRamp2(c(0, 1, 4), c("#1E295C", "white", "#C02025")))
pdf("HAQER_enrichment.pdf",width=5,height=5)
print(p1)
dev.off()

# ccc = colorRamp2(c(0, 1, 2), c("#0000FF", "white", "#FF0000"))
# library(pheatmap)
# pdf("HAQER_enrichment1.pdf",width=5,height=5)
# pheatmap(M_enrichment_C1,cluster_rows = FALSE, cluster_cols = FALSE,color= ccc,angle_col = 0,display_numbers = TRUE,show_rownames = T,cellwidth= 60)
# dev.off()

# color2 = colorRampPalette(c("#0000FF", "white", "#FF0000"))(30)

# pdf("HAQER_enrichment1.pdf",width=5,height=5)
# pheatmap(M_enrichment_C1,cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = TRUE,cellwidth= 60)
# dev.off()
