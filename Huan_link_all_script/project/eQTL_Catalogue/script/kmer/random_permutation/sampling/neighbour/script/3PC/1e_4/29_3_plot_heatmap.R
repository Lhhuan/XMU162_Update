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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/ChromHMM/")
enrichment <-read.table("/share/data0/QTLbase/huan/eQTL_Catalogue/ChromHMM/output/merge_15_state_enrichment.txt.gz",header = T,sep = "\t") %>% as.data.frame()
sample_id <- read_excel("chromHmm_1581.xlsx", sheet = "IDs")%>% as.data.frame()

enrichment$EnrichFDR <- p.adjust(enrichment$EnrichPValue,method="fdr")
enrichment$DepleteFDR <- p.adjust(enrichment$DepletePValue,method="fdr")
enrichment$SignificantEnrichment <-NA
for(i in 1:nrow(enrichment)){
    if(enrichment[i,"EnrichFDR"]<0.05 |enrichment[i,"DepleteFDR"]<0.05 ){
        enrichment[i,"SignificantEnrichment"] <-enrichment[i,"Enrichment"]
    }
}

for (cc in c("C1","C2","C3","C4","C5","C6")){
    enrichment_C1 <- filter(enrichment,cluster==cc)
    M_enrichment_C1 <- dcast(enrichment_C1[,c("state","sample","SignificantEnrichment")],state~sample)
    rownames(M_enrichment_C1)<-M_enrichment_C1[,1]
    M_enrichment_C1<-M_enrichment_C1[,-1]
    M_enrichment_C1 <- M_enrichment_C1[c("1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","8_ZNF_Rpts","9_Het","10_TssBiv","11_BivFlnk","12_EnhBiv","13_ReprPC","14_ReprPCWk","15_Quies"),]

    
    p1 <- Heatmap(as.matrix(M_enrichment_C1), name = "Overlap enrichment", row_names_side = "left", column_names_side = "top",
        cluster_rows = FALSE, cluster_columns = FALSE,
        width = ncol(M_enrichment_C1)*unit(2, "mm"),height = nrow(M_enrichment_C1)*unit(2, "mm"),
        row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 5),
        col = colorRamp2(c(0, 1, 4), c("#0000FF", "white", "#FF0000")) ,
        # Legend(labels_gp = gpar(fontsize = 5, fontface = "bold")))
        heatmap_legend_param = list(title_gp = gpar(fontsize = 6.5)))
        # col = colorRamp2(c(0, 1, 4), c("#1E295C", "white", "#C02025")))
    pdf(paste0(cc,"_enrichment.pdf"),width=11.5,height=2)
    print(p1)
    dev.off()

}


#======================================
for (cc in c("C1","C2","C3","C4","C5","C6")){
    enrichment_C1 <- filter(enrichment,cluster==cc)
    M_enrichment_C1 <- dcast(enrichment_C1[,c("state","sample","SignificantEnrichment")],state~sample)
    rownames(M_enrichment_C1)<-M_enrichment_C1[,1]
    M_enrichment_C1$state <-gsub("1_TssA|2_TssAFlnk","TSS",M_enrichment_C1$state)
    M_enrichment_C1$state <-gsub("3_TxFlnk|4_Tx|5_TxWk","Tx",M_enrichment_C1$state)
    M_enrichment_C1$state <-gsub("6_EnhG|7_Enh","Enh",M_enrichment_C1$state)
    M_enrichment_C1$state <-gsub("10_TssBiv|11_BivFlnk|12_EnhBiv","Bivalent",M_enrichment_C1$state)
    M_enrichment_C1$state <-gsub("9_Het|13_ReprPC|14_ReprPCWk","Repr",M_enrichment_C1$state)
    M_enrichment_C1$state <-gsub("8_ZNF_Rpts|15_Quies","Other",M_enrichment_C1$state)
    # M_enrichment_C1<-M_enrichment_C1[,-1]
    M_enrichment_C1 <- M_enrichment_C1[c("1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","10_TssBiv","11_BivFlnk","12_EnhBiv","9_Het","13_ReprPC","14_ReprPCWk","8_ZNF_Rpts","15_Quies"),]

    # col_fun = colorRamp2(c("TSS", "Tx", "Enh","Bivalent","Repr","Other"), 
    # c("#B4C6E7", "#F8CBAD", "#DBDBDB","#BDD7EE","#C6E0B4","#FFE699"))
    col_fun = list(Category=c("1_TssA"="#B4C6E7","2_TssAFlnk"="#B4C6E7", "3_TxFlnk"="#F8CBAD","4_Tx"="#F8CBAD","5_TxWk"="#F8CBAD", "6_EnhG"="#DBDBDB","7_Enh"="#DBDBDB","10_TssBiv"="#BDD7EE","11_BivFlnk"="#BDD7EE","12_EnhBiv"="#BDD7EE","9_Het"="#C6E0B4","13_ReprPC"="#C6E0B4","14_ReprPCWk"="#C6E0B4","8_ZNF_Rpts"="#FFE699","15_Quies"="#FFE699"))
    df=data.frame(Category=rownames(M_enrichment_C1))

    p1 = Heatmap(as.matrix(M_enrichment_C1[,2:ncol(M_enrichment_C1)]), name = "Overlap enrichment",                row_names_side = "right", column_names_side = "top",
        column_title = cc,column_title_gp = gpar(fontsize = 6.5),
        cluster_rows = FALSE, cluster_columns = FALSE,
        width = ncol(M_enrichment_C1)*unit(1, "mm"),height = nrow(M_enrichment_C1)*unit(2, "mm"),
        row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 3),
        col = colorRamp2(c(0, 1, 4), c("#0000FF", "white", "#FF0000")) ,
        # Legend(labels_gp = gpar(fontsize = 5, fontface = "bold")))
        heatmap_legend_param = list(title_gp = gpar(fontsize = 6.5),
        legend_direction = "horizontal",
        legend_width = unit(3, "cm") ))#, 
        # title_position = "lefttop"))
    p2=rowAnnotation(df = df, col = col_fun,show_legend = FALSE,width=unit(2,"mm"))
        # col = colorRamp2(c(0, 1, 4), c("#1E295C", "white", "#C02025")))
    p=p2+p1

    pdf(paste0(cc,"_enrichment.pdf"),width=6,height=3)
    draw(p,heatmap_legend_side = "bottom")
    # print(p1)
    dev.off()
}


#============合并样本
rs <-data.frame()
for (cc in c("C1","C2","C3","C4","C5","C6")){
    enrichment_C1 <- filter(enrichment,cluster==cc)
    M_enrichment_C1 <- dcast(enrichment_C1[,c("state","sample","SignificantEnrichment")],state~sample)
    rownames(M_enrichment_C1) <- M_enrichment_C1$state 
    M_enrichment_C1 <- M_enrichment_C1[,-1]
    dat <-data.frame(state=rownames(M_enrichment_C1),mean_enrichment= rowMeans(M_enrichment_C1,na.rm=TRUE),cluster=cc)
    rownames(dat) <-NULL
    rs <- bind_rows(rs,dat)
}

fdat <- dcast(rs[,c("state","cluster","mean_enrichment")],state~cluster)
rownames(fdat) <- fdat$state
fdat <- fdat[,-1]
fdat <- fdat[c("1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","10_TssBiv","11_BivFlnk","12_EnhBiv","9_Het","13_ReprPC","14_ReprPCWk","8_ZNF_Rpts","15_Quies"),]
col_fun = list(Category=c("1_TssA"="#B4C6E7","2_TssAFlnk"="#B4C6E7", "3_TxFlnk"="#F8CBAD","4_Tx"="#F8CBAD","5_TxWk"="#F8CBAD", "6_EnhG"="#DBDBDB","7_Enh"="#DBDBDB","10_TssBiv"="#BDD7EE","11_BivFlnk"="#BDD7EE","12_EnhBiv"="#BDD7EE","9_Het"="#C6E0B4","13_ReprPC"="#C6E0B4","14_ReprPCWk"="#C6E0B4","8_ZNF_Rpts"="#FFE699","15_Quies"="#FFE699"))
df=data.frame(Category=rownames(fdat))

p1 = Heatmap(as.matrix(fdat), name = "Overlap enrichment",                row_names_side = "right", 
# column_names_side = "top",
    column_names_rot=0,
    cluster_rows = FALSE, cluster_columns = FALSE,
    width = ncol(fdat)*unit(10, "mm"),height = nrow(fdat)*unit(5, "mm"),
    row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 10),
    col = colorRamp2(c(0, 1, 10), c("#0000FF", "white", "#FF0000")) ,
    # Legend(labels_gp = gpar(fontsize = 5, fontface = "bold")))
    heatmap_legend_param = list(title_gp = gpar(fontsize = 6.5),
    legend_direction = "horizontal",
    legend_width = unit(3, "cm") ))#, 
    # title_position = "lefttop"))
p2=rowAnnotation(df = df, col = col_fun,show_legend = FALSE,width=unit(2,"mm"))
    # col = colorRamp2(c(0, 1, 4), c("#1E295C", "white", "#C02025")))
p=p2+p1

pdf("29_3_mean_enrichment.pdf",width=5,height=5)
draw(p,heatmap_legend_side = "bottom")
# print(p1)
dev.off()

#==========================pheatmap 
# color2 = colorRampPalette(c('#c6dbef','#6baed6','#2171b5'))(50)
color2 = colorRamp2(c(0, 1, 10), c("#0000FF", "white", "#FF0000"))
library(pheatmap)
pdf("29_3_mean_enrichment_pheatmap.pdf",width=5,height=5)
pheatmap(as.matrix(fdat),cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = TRUE,show_rownames = T,show_colnames = T)
dev.off()
