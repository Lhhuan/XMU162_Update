library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tibble)
library(R.utils)
library(reshape2)

markers = c("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3","CTCF","TFBS","CHROMATIN_Accessibility")

setwd("/home/huanhuan/project/eQTL_Catalogue/output/annotation/extend_18snp/signalValue/egene/")

gene_signal<-function(marker= NULL){
    file_name <-paste0(marker,"_max_mean_median_09_egene_pos_5e_8_sorted.bed.gz")
    org <- read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
    org <-org[,c(1:4,7)]
    org$marker <-marker
    return(org)
}

tmp <- lapply(markers,gene_signal)
rs<-do.call(rbind,tmp)
gene <-dcast(rs[,c("Chr","start","end","egene","marker","max_signalvalue")], Chr+start+end+egene~marker)
colnames(gene)[4] <-"Gene"

setwd("/home/huanhuan/project/eQTL_Catalogue/output/annotation/extend_18snp/signalValue/hotspot/")
hotspot_signal<-function(marker= NULL){
    file_name <-paste0(marker,"_max_mean_median_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz")
    org <- read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
    org <-org[,c(1:3,6)]
    org$marker <-marker
    return(org)
}
tmp2 <-lapply(markers,hotspot_signal)
rs2<-do.call(rbind,tmp2)
hotspot <-dcast(rs2[,c("Chr","start","end","marker","max_signalvalue")], Chr+start+end~marker)
hotspot$Hotspot <- paste(hotspot$Chr,hotspot$start,hotspot$end,sep="_")
# hotspot <-hotspot[,4:14]


setwd("/home/huanhuan/project/eQTL_Catalogue/script/prediction/YY/chr1/output/")

gene_idx <- read.table("01_gene_idx.txt.gz",header = T,sep = "\t") %>% as.data.frame()
hotspot_idx <- read.table("01_hotspot_idx.txt.gz",header = T,sep = "\t") %>% as.data.frame()

gene_idx <-left_join(gene_idx,gene,by="Gene")
hotspot_idx <-left_join(hotspot_idx,hotspot,by="Hotspot")

save(gene_idx,file="04_gene_anno_idx.Rdata")
save(hotspot_idx,file="04_hotspot_anno_idx.Rdata")
colnames(gene_idx)[1:2]<-c("Name","node_idx")
colnames(hotspot_idx)[1:2]<-c("Name","node_idx")
idx <-bind_rows(hotspot_idx,gene_idx)
save(idx,file="04_hotspot_egene_node_anno_idx.Rdata")
idx<-idx[,2:15]

interaction<- read.table("02_hotspot_egene_idx.txt.gz",header = T,sep = "\t") %>% as.data.frame()
set.seed(1)
# positive<- interaction[sample(1:nrow(interaction),10000, replace = TRUE),]
positive  <-interaction
set.seed(1)
random_h <- hotspot_idx[sample(1:nrow(hotspot_idx), 2*nrow(positive), replace = TRUE),"idx"]
random_e <- gene_idx[sample(1:nrow(gene_idx), 2*nrow(positive), replace = TRUE),"idx"]

rs3 <-data.frame(Hotspot_idx=random_h,Egene_idx=random_e)
potential_negative <- anti_join(rs3,interaction,by=c("Hotspot_idx","Egene_idx"))
set.seed(1)
negative <-potential_negative[sample(1:nrow(potential_negative),nrow(positive),replace = FALSE),]
negative$lable <-0
positive$lable <-1
finter<-bind_rows(positive,negative)
write.table(finter,"edges.csv",row.names = F, col.names = T,quote =F,sep=",")
write.table(finter,"edges.txt",row.names = F, col.names = T,quote =F,sep="\t")

set.seed(1)
node_features <-filter(idx,node_idx %in% unique(c(finter$Hotspot_idx,finter$Egene_idx)))
node_features <-node_features[sample(1:nrow(node_features),nrow(node_features),replace = FALSE),]
write.table(node_features,"node.csv",row.names = F, col.names = T,quote =F,sep=",")
write.table(node_features,"node.txt",row.names = F, col.names = T,quote =F,sep="\t")

