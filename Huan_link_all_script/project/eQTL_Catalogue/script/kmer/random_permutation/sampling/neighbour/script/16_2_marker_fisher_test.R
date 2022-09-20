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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/markers/")
Cluster <-read.table("../../11_chr1_louvain_pca5_k500.txt",header = T,sep = "\t") %>% as.data.frame()
cluster_n <- Cluster%>%group_by(cluster)%>%summarise(hotspot_num=n())%>%data.frame()
org <-read.table("/home/huanhuan/project/eQTL_Catalogue/output/enrichment/hotspot_cutoff_0.176_marker_jaccard_index_hotspot.txt.gz",header = T,sep = "\t") %>% as.data.frame()

chr1_org <- filter(org,hotspot_chr=="chr1" & jaacard_index>0)
chr1_org$hotspot <-paste0(chr1_org$hotspot_chr,":",chr1_org$hotspot_start,"-",chr1_org$hotspot_end)
chr1_org<-left_join(chr1_org,Cluster[,c("cluster","hotspot")],by="hotspot")
chr1_org_n <- chr1_org %>% group_by(cluster,Marker)%>%summarise(hotspot_hit_number=n())%>%data.frame()
chr1_org_n <-left_join(chr1_org_n,cluster_n,by="cluster")
chr1_org_n$hotspot_not_hit <- chr1_org_n$hotspot_num -chr1_org_n$hotspot_hit_number

anno_dir <-"/home/huanhuan/project/eQTL_Catalogue/output/annotation/extend_18snp"

rs <-data.frame()
histone=c("H3K27ac","H3K4me1","H3K4me3","H3K9ac","H3K36me3","H3K27me3","H3K9me3")
for (marker in histone){
    dat <- read.table(paste0(anno_dir,"/",marker,"_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz"),header = F,sep = "\t") %>% as.data.frame()
    dat <-filter(dat,V1=="chr1")
    colnames(dat)[1:6] <-c("h_chr","h_start","h_end","m_chr","m_start","m_end")
    dat$hotspot <-paste0(dat$h_chr,":",dat$h_start,"-",dat$h_end)
    dat$marker <-paste0(dat$m_chr,":",dat$m_start,"-",dat$m_end)
    dat <-left_join(dat,Cluster[,c("cluster","hotspot")],by="hotspot")
    dat_n <-dat %>%group_by(cluster)%>%summarise(anno_marker=n())%>%data.frame()
    fanno <- read.table(paste0("/share/data0/GTEx/annotation/ROADMAP/sample/merge/hg38/",marker,"_sorted_merge.bed.gz"),header = F,sep = "\t") %>% as.data.frame()
    fanno <-filter(fanno,V1=="chr1")
    dat_n$number_of_total_marker <-nrow(fanno)
    dat_n$Marker <-marker
    rs<-bind_rows(rs,dat_n)
    print(marker)
}

cistromeDB=c("Human_FACTOR","Human_CHROMATIN_Accessibility")
for (marker in cistromeDB){
    marker1 <-gsub("Human_CHROMATIN_Accessibility","CHROMATIN_Accessibility",marker)
    marker1 <-gsub("Human_FACTOR","TFBS",marker1)
    dat <- read.table(paste0(anno_dir,"/",marker1,"_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz"),header = F,sep = "\t") %>% as.data.frame()
    dat <-filter(dat,V1=="chr1")
    colnames(dat)[1:6] <-c("h_chr","h_start","h_end","m_chr","m_start","m_end")
    dat$hotspot <-paste0(dat$h_chr,":",dat$h_start,"-",dat$h_end)
    dat$marker <-paste0(dat$m_chr,":",dat$m_start,"-",dat$m_end)
    dat <-left_join(dat,Cluster[,c("cluster","hotspot")],by="hotspot")
    dat_n <-dat %>%group_by(cluster)%>%summarise(anno_marker=n())%>%data.frame()
    fanno <- read.table(paste0("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/",marker,"/merge_pos_info_narrow_peak_sort_union.bed.gz"),header = F,sep = "\t") %>% as.data.frame()
    fanno <-filter(fanno,V1=="chr1")
    dat_n$number_of_total_marker <-nrow(fanno)
    dat_n$Marker <-marker1
    rs<-bind_rows(rs,dat_n)
    print(marker)
}

CTCF=c("CTCF")
for (marker in CTCF){
    dat <- read.table(paste0(anno_dir,"/",marker,"_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz"),header = F,sep = "\t") %>% as.data.frame()
    dat <-filter(dat,V1=="chr1")
    colnames(dat)[1:6] <-c("h_chr","h_start","h_end","m_chr","m_start","m_end")
    dat$hotspot <-paste0(dat$h_chr,":",dat$h_start,"-",dat$h_end)
    dat$marker <-paste0(dat$m_chr,":",dat$m_start,"-",dat$m_end)
    dat <-left_join(dat,Cluster[,c("cluster","hotspot")],by="hotspot")
    dat_n <-dat %>%group_by(cluster)%>%summarise(anno_marker=n())%>%data.frame()
    fanno <- read.table("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/CTCF/normal_cell_line/hg38/05_normal_cell_line_ctcf_sort_union_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame()
    fanno <-filter(fanno,V1=="chr1")
    dat_n$number_of_total_marker <-nrow(fanno)
    dat_n$Marker <-marker
    rs<-bind_rows(rs,dat_n)
    print(marker)
}

fdat <-left_join(chr1_org_n,rs,by=c("cluster","Marker"))
fdat$number_of_not_anno_mark <-fdat$number_of_total_marker -fdat$anno_marker
write.table(fdat,"16_2_marker_anno_for_fisher_test.txt",row.names = F, col.names = T,quote =F,sep="\t")


all_fisher <-data.frame()
for(marker in unique(fdat$Marker)){
    for(i in unique(fdat$cluster)){
        dat <- filter(fdat,Marker==marker & cluster==i)
        pre_data <-matrix(c(dat$hotspot_hit_number,dat$anno_marker, dat$hotspot_not_hit, dat$number_of_not_anno_mark), nrow = 2)
        fisher_test_result <-fisher.test(pre_data)
        p_value =fisher_test_result$p.value
        conf_int1 <-fisher_test_result$conf.int[1]
        conf_int2 <-fisher_test_result$conf.int[2]
        od <-fisher_test_result$estimate
        names(od)=NULL
        if(p_value <= 0.0001){
            annotation = "****"
        }else if(p_value <= 0.001){
            annotation = "***"
        }else if(p_value <= 0.01){
            annotation = "**"
        }else if(p_value <= 0.05){
            annotation = "*"
        }else{
            annotation = "ns"
        }
        tmp <-data.frame(clutser= i, marker =marker,  p_value =p_value,conf_int_bottom=conf_int1,conf_int_up=conf_int2,OR = od,significant=annotation)
        all_fisher <- bind_rows(all_fisher,tmp)
        print(marker)
        print(i)
    }
}

write.table(all_fisher,"16_2_marker_cluster_fisher_test.txt",row.names = F, col.names = T,quote =F,sep="\t")

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 12),
                                                axis.title.x = element_text(size = 12),
                                                axis.line = element_line(color = "black"),
                                                axis.text.x = element_text(color = "black",size = 12),axis.text.y = element_text(color = "black",size = 12),
                                                plot.title = element_text(hjust = 0.5),legend.position="none")


for(Marker in unique(all_fisher$marker)){
    dat <-filter(all_fisher,marker==Marker)
    Marker <-gsub("CHROMATIN_Accessibility","CA",Marker)
    dat$clutser <-paste0("C",dat$clutser)
    dat$clutser <-factor(dat$clutser, levels=dat$clutser)
    dat <-dat%>%arrange(OR)
    dat$FDR <-p.adjust(dat$p_value,method="fdr")
    p1<-ggplot(data=dat, aes(x=reorder(clutser,OR),y=log(OR), ymin=log(conf_int_bottom) , ymax=log(conf_int_up),color=clutser))+
        geom_pointrange(size = 0.5)+
        scale_color_manual(values=c("#42B540FF","#0099B4FF","#A593E0","#fcbe32","#F16B6F","#F68657","#5A9367"))+
        # scale_y_continuous(limits= c(1.5,6), breaks= seq(1.5,6,1.5))+
        coord_flip()+
        ylab("Odds ratio (log scale)")+xlab("Cluster")+
        p_theme+
        geom_text(aes(x=clutser,y=6,label = significant))+
        ggtitle(Marker)
    pdf(paste0("16_2_",Marker,"_cluster_fisher_test.pdf"),height = 3.5,width = 3.8)
    print(p1)
    dev.off()
    print(Marker)
}


markers=c("TFBS")
for(Marker in markers){
    dat <-filter(all_fisher,marker==Marker)
    Marker <-gsub("CHROMATIN_Accessibility","CA",Marker)
    dat$clutser <-paste0("C",dat$clutser)
    dat$clutser <-factor(dat$clutser, levels=dat$clutser)
    dat <-dat%>%arrange(OR)
    dat$FDR <-p.adjust(dat$p_value,method="fdr")
    p1<-ggplot(data=dat, aes(x=reorder(clutser,OR),y=log(OR), ymin=log(conf_int_bottom) , ymax=log(conf_int_up),color=clutser))+
        geom_pointrange(size = 0.5)+
        scale_color_manual(values=c("#42B540FF","#0099B4FF","#A593E0","#fcbe32","#F16B6F","#F68657","#5A9367"))+
        scale_y_continuous(limits= c(5,7.5))+
        coord_flip()+
        ylab("Odds ratio (log scale)")+xlab("Cluster")+
        p_theme+
        geom_text(aes(x=clutser,y=7.5,label = significant))+
        ggtitle(Marker)
    pdf(paste0("16_2_",Marker,"_cluster_fisher_test.pdf"),height = 3.5,width = 3.8)
    print(p1)
    dev.off()
    print(Marker)
}