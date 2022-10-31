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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/gene/")
all_gene <-read.table("/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_eQTL_5e_8.bed.gz",header = F,sep = "\t") %>% as.data.frame()
Cluster <-read.table("../../11_whole_genome_leiden_pca3_k50_resolution1e-04.txt",header = T,sep = "\t") %>% as.data.frame()
tmp1 <-strsplit(Cluster$hotspot,":")
tmp2 <- do.call(rbind,tmp1)%>%data.frame()
Cluster$CHR <- tmp2$X1
cluster_chr <-Cluster %>%group_by(cluster,CHR)%>%summarise(count=n())%>%data.frame()
cluster_n <- Cluster%>%group_by(cluster)%>%summarise(cluster_count=n())%>%data.frame()
cluster_chr <-left_join(cluster_chr,cluster_n,by="cluster")
cluster_chr$ratio <- cluster_chr$count/cluster_chr$cluster_count


all_gene$hotspot <-paste0(all_gene$V1,":",all_gene$V2,"-",all_gene$V3)
all_gene$eqtl <- paste0(all_gene$V4,":",all_gene$V5,"-",all_gene$V6)
all_gene$length <- all_gene$V3 - all_gene$V2
all_gene<-left_join(all_gene,Cluster[,c("cluster","hotspot")],by="hotspot")
colnames(all_gene)[7]<-"ENSEMBL"
all_gene<-all_gene[,c("cluster","hotspot","eqtl","length","ENSEMBL")]
all_gene_count <-all_gene%>%group_by(hotspot,cluster,length)%>%summarise(eqtl_gene_count=n())%>%data.frame()
all_gene_count$adjust_eqtl_gene_count <-all_gene_count$eqtl_gene_count/all_gene_count$length*1000
all_gene_count1 <-filter(all_gene_count,cluster %in%c(1:6))

p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5))

library(clusterProfiler)
library(org.Hs.eg.db)

data = all_gene[,"ENSEMBL"]%>%unique()
data = as.vector(data)
annots <- select(org.Hs.eg.db, keys=data, 
                columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")
annots <-annots[!duplicated(annots$ENSEMBL),]
all_gene <-left_join(all_gene,annots,by="ENSEMBL")
all_gene1 <-unique(all_gene[,c("cluster","ENSEMBL","SYMBOL","ENTREZID")])

#======
bbb<-list()
for(i in 1:6){
    itmp <-filter(all_gene1,cluster==i)
    bbb[[i]]<-unique(itmp$ENSEMBL)
    names(bbb)[i] <-paste0("C",i)
    print(i)
    # print(nrow(filter(all_gene1,cluster==i)%>%dplyr::select(ENSEMBL)%>%unique()))
}
library(venn)
pdf("13_egene_venn.pdf",width=3.1,height=3.1)
venn(bbb, ilab = TRUE, zcolor = "style")
dev.off()

rs <-data.frame()
for(j in c(1:length(bbb))){
    b <-bbb
    all <-bbb[[j]]
    CC <- names(bbb[j])
    CC <-gsub("C","",CC)
    b[[j]]=NULL
    names(b)=NULL
    other <-unlist(b)
    specific <- setdiff(all,other)%>%data.frame()
    colnames(specific) <-"ENSG"
    specific$cluster <-j 
    rs <-bind_rows(rs,specific)
}

gene_cluster_n <-all_gene1%>%group_by(ENSEMBL)%>%summarise(count=n())%>%data.frame()
share_gene <-filter(gene_cluster_n,count==6)
share_gene$cluster <-"Share"
colnames(share_gene)[1] <-"ENSG"
# random <- data.frame()
# for(i in 1:6){
#     c_specific <-filter(rs,cluster==i)
#     set.seed(123)
#     tmp3 <-c_specific[sample(1:nrow(c_specific), 20, replace = FALSE),]
#     random <- bind_rows(random,tmp3)
# }
# set.seed(123)
# tmp3 <-share_gene[sample(1:nrow(share_gene), 20, replace = FALSE),]
# random$cluster <- as.character(random$cluster)
# random <- bind_rows(random,tmp3[,c(1,3)])
#==============
gene_exp <-read.table("/share/data0/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",header = T,sep = "\t",skip=2) %>% as.data.frame()
# gene_exp$Name<-gsub("\\..*","",gene_exp$Name)
colnames(gene_exp) <-gsub("\\.","_",colnames(gene_exp))
colnames(gene_exp) <-gsub("___","_",colnames(gene_exp))
colnames(gene_exp) <-gsub("__","_",colnames(gene_exp))
colnames(gene_exp)[1] <-"ENSG"
gene_exp$sd <- apply(gene_exp[,3:ncol(gene_exp)],1,sd)
gene_exp <-gene_exp[order(-gene_exp$sd),]
rownames(gene_exp)<-gene_exp$ENSG
library(pheatmap)
color2 = colorRampPalette(c('#c6dbef','#6baed6','#2171b5'))(50)
pdf("13_1_gene_expression_tissue_heatmap_top200.pdf",width=13,height=20)
pheatmap(as.matrix(gene_exp[1:200,3:56]),cluster_rows = TRUE, cluster_cols = FALSE,angle_col = 90,color= color2,display_numbers = FALSE,show_rownames = T,cellwidth= 10)
dev.off()
log_gene_exp <- log(gene_exp[,3:56])
pdf("13_1_gene_expression_tissue_heatmap_top200_log.pdf",width=13,height=20)
pheatmap(as.matrix(log_gene_exp[1:200,]),cluster_rows = TRUE, cluster_cols = FALSE,angle_col = 90,color= color2,display_numbers = FALSE,show_rownames = T,cellwidth= 10)
dev.off()





random <-filter(random,ENSG %in%gene_exp$ENSG)
fgene <-left_join(random,gene_exp,by="ENSG")
fgene$Name <-paste(fgene$cluster,fgene$ENSG,sep="-")
rownames(fgene) <-fgene$Name


# fgene <-fgene[random1$ENSG,]
fgene <-fgene[,-c(1:3,58)]
library(pheatmap)
color2 = colorRampPalette(c('#c6dbef','#6baed6','#2171b5'))(50)
pdf("13_1_gene_expression_tissue_heatmap.pdf",width=13,height=20)
pheatmap(as.matrix(fgene),cluster_rows = TRUE, cluster_cols = FALSE,angle_col = 90,color= color2,display_numbers = FALSE,show_rownames = T,cellwidth= 10)
dev.off()