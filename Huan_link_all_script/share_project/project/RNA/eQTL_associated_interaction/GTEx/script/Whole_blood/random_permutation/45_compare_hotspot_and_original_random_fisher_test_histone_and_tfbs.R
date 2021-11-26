library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(Seurat)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/enrichment/interval_18/ALL/")
hotspot<-read.table("hotspot_cutoff_0.176_marker_jaccard_index.txt.gz",header = T,sep = "\t") %>% as.data.frame()
# random<-fread("0_0.176_jaccard_index_histone_marker_1000.txt.gz",header = T,sep = "\t")
random<-read.table("original_random_jaccard_index_histone_marker_1000.txt.gz",header = T,sep = "\t") %>% as.data.frame()

hotspot$Random_number <- 1
hotspot$Class <- "hotspot"
random$Class <- "random"

rs<- bind_rows(hotspot,random)
# marks <-c("H3K27ac","H3K9ac","H3K36me3","H3K4me1","H3K4me3","H3K27me3","H3K9me3")
marks <-c("H3K27ac","H3K9ac","H3K36me3","H3K4me1","H3K4me3","H3K27me3","H3K9me3","CHROMATIN_Accessibility","TFBS","CTCF")
#-------------two side
all_fisher <-data.frame()
for(state in marks){
    aa <-filter(rs,Marker==state)%>%filter(jaacard_index >0)%>%group_by(Class)%>%summarise(table(Class)%>%as.data.frame())%>%as.data.frame()
    a <-filter(aa,Class=="hotspot")$Freq
    r_c<-filter(aa,Class=="random")$Freq
    pre_data <-matrix(c(a,r_c, nrow(hotspot)/10-a, nrow(hotspot)/10*1000-r_c), nrow = 2)

    # fisher_test_result <-fisher.test(pre_data,alternative = "greater")
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
    # p_value<-as.character(p_value)
    tmp <-data.frame(marker =state,  p_value =p_value,conf_int_bottom=conf_int1,conf_int_up=conf_int2,OR = od,significant=annotation)
    all_fisher <- bind_rows(all_fisher,tmp)
    print(state)
}

write.table(all_fisher,"./figure/original_random/compare_original_random_jaacard_index_fisher_test_histone_tfbs_two_side.txt",row.names = F, col.names = T,quote =F,sep="\t")

#---------------alternative = "greater"
all_fisher <-data.frame()
for(state in marks){
    aa <-filter(rs,Marker==state)%>%filter(jaacard_index >0)%>%group_by(Class)%>%summarise(table(Class)%>%as.data.frame())%>%as.data.frame()
    a <-filter(aa,Class=="hotspot")$Freq
    r_c<-filter(aa,Class=="random")$Freq
    pre_data <-matrix(c(a,r_c, nrow(hotspot)/10-a, nrow(hotspot)/10*1000-r_c), nrow = 2)

    fisher_test_result <-fisher.test(pre_data,alternative = "greater")
    # fisher_test_result <-fisher.test(pre_data)

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
    # p_value<-as.character(p_value)
    tmp <-data.frame(marker =state,  p_value =p_value,conf_int_bottom=conf_int1,conf_int_up=conf_int2,OR = od,significant=annotation)
    all_fisher <- bind_rows(all_fisher,tmp)
    print(state)
}

write.table(all_fisher,"./figure/original_random/compare_original_random_jaacard_index_fisher_test_greater.txt",row.names = F, col.names = T,quote =F,sep="\t")

