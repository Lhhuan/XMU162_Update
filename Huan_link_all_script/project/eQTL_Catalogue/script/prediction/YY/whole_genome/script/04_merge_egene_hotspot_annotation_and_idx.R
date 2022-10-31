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

kmer<-read.csv("/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/egene_pos_0.05_6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()
load("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/06_permutation_test_1000_sig_kmer.Rdata")
sigK <-fdat0

colnames(kmer)[1] <-"hotspot"
kmer$hotspot <-gsub(">","",kmer$hotspot)
rownames(kmer)=kmer$hotspot
kmer <-kmer[,-1]
egene_Skmer <-kmer[,which(colnames(kmer) %in% sigK$seq)]
egene_Skmer$pos <- rownames(egene_Skmer)

load("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/08_permutation_wilocx_overlap_sig_kmer_0_1000.Rdata")
hotspot_Skmer <- Sorg
hotspot_Skmer$name <- rownames(hotspot_Skmer)

#==============================================================================================================================
markers = c("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3","CTCF","TFBS","CHROMATIN_Accessibility")

setwd("/home/huanhuan/project/eQTL_Catalogue/output/annotation/extend_18snp/signalValue/egene/")

gene_signal<-function(marker= NULL){
    file_name <-paste0(marker,"_max_mean_median_egene0.05_pos_sorted.bed.gz")
    org <- read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
    org <-org[,c(1:4,7)]
    org$marker <-marker
    return(org)
}

tmp <- lapply(markers,gene_signal)
rs<-do.call(rbind,tmp)
gene <-dcast(rs[,c("Chr","start","end","egene","marker","max_signalvalue")], Chr+start+end+egene~marker)
colnames(gene)[4] <-"name"
gene$pos <- paste0(gene$Chr,":",gene$start,"-",gene$end)
gene <-left_join(gene,egene_Skmer,by ="pos")
gene <- gene%>%dplyr::select(-pos)


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

hotspot$name <- paste0(hotspot$Chr,":",hotspot$start,"-",hotspot$end)
hotspot <- left_join(hotspot,hotspot_Skmer,by="name")
hotspot$name <- paste(hotspot$Chr,hotspot$start,hotspot$end,sep="_")
# hotspot <-hotspot[,4:14]
#====================================================================================
setwd("/home/huanhuan/project/eQTL_Catalogue/script/prediction/YY/whole_genome/output/")

node_edge <-function(chr=NULL){
    all_idx <- read.table(paste0("01_",chr,"_hotspot_egene_idx.txt.gz"),header = T,sep = "\t") %>% as.data.frame()
    G_idx <-all_idx[grep("ENSG",all_idx$name),] 
    h_idx <-all_idx[grep("chr",all_idx$name),]
    gene_idx <-left_join(G_idx,gene,,by="name")
    # aaaa <-inner_join(h_idx,hotspot,by="name")
    hotspot_idx <-left_join(h_idx,hotspot,by="name")
    colnames(gene_idx)[1:2]<-c("Name","node_idx")
    colnames(hotspot_idx)[1:2]<-c("Name","node_idx")
    idx <-bind_rows(hotspot_idx,gene_idx)
    save(idx,file=paste0("04_",chr,"_hotspot_egene_node_anno_idx.Rdata"))
    idx<-idx[,c(2,1,3:411)]
    #---------------------------------------
    interaction<- read.table(paste0("02_",chr,"_hotspot_egene_idx.txt.gz"),header = T,sep = "\t") %>% as.data.frame()
    set.seed(1)
    # positive<- interaction[sample(1:nrow(interaction),10000, replace = TRUE),]
    positive  <-interaction
    set.seed(1)
    random_h <- hotspot_idx[sample(1:nrow(hotspot_idx), 2*nrow(positive), replace = TRUE),"node_idx"]
    random_e <- gene_idx[sample(1:nrow(gene_idx), 2*nrow(positive), replace = TRUE),"node_idx"]

    rs3 <-data.frame(Hotspot_idx=random_h,Egene_idx=random_e)
    potential_negative <- anti_join(rs3,interaction,by=c("Hotspot_idx","Egene_idx"))
    set.seed(1)
    negative <-potential_negative[sample(1:nrow(potential_negative),nrow(positive)*1.5,replace = FALSE),]
    negative$lable <-0
    positive$lable <-1
    finter<-bind_rows(positive,negative)
    write.table(finter,paste0("./raw/04_",chr,"_edges.csv"),row.names = F, col.names = T,quote =F,sep=",")
    write.table(finter,paste0("./raw/04_",chr,"_edges.txt"),row.names = F, col.names = T,quote =F,sep="\t")

    set.seed(1)
    node_features <-filter(idx,node_idx %in% unique(c(finter$Hotspot_idx,finter$Egene_idx)))
    node_features <-node_features[order(node_features$node_idx),]
    write.table(node_features,paste0("./raw/04_",chr,"_node.csv"),row.names = F, col.names = T,quote =F,sep=",")
    write.table(node_features,paste0("./raw/04_",chr,"_node.txt"),row.names = F, col.names = T,quote =F,sep="\t")
    print(chr)
}

lapply(paste0("chr",c(1:22)),node_edge)

