library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(Seurat)

setwd("/home/huanhuan/project/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/enrichment")
hotspot<-read.table("hotspot_cutoff_0.176_marker_jaccard_index.txt.gz",header = T,sep = "\t") %>% as.data.frame()
# random<-fread("0_0.176_jaccard_index_histone_marker_1000.txt.gz",header = T,sep = "\t")
random<-read.table("0_0.176_jaccard_index_marker_1000.txt.gz",header = T,sep = "\t") %>% as.data.frame()


# active_mark <-c("H3K27ac","H3K9ac","H3K36me3","H3K4me1","H3K4me3")
# repressed_mark <-c("H3K27me3","H3K9me3")

hotspot$Random_number <- 1
hotspot$Class <- "hotspot"
random$Class <- "random"

rs<- bind_rows(hotspot,random)
marks <-c("H3K27ac","H3K9ac","H3K36me3","H3K4me1","H3K4me3","H3K27me3","H3K9me3","CHROMATIN_Accessibility","TFBS","CTCF")
#-------------two side
# all_fisher <-data.frame()
# for(state in marks){
fisher_two <-function(state){
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
    # all_fisher <- bind_rows(all_fisher,tmp)
    print(state)
    return(tmp)
}
re  <-lapply(marks,fisher_two)

result <-do.call(rbind,re)
write.table(result,"./figure/0_0.176/compare_0_0.176_jaacard_index_fisher_test_histone_tfbs_two_side.txt",row.names = F, col.names = T,quote =F,sep="\t")

#---------------alternative = "greater"
# all_fisher <-data.frame()
# for(state in marks){
fisher_greater <-function(state){
    aa <-filter(rs,Marker==state)%>%filter(jaacard_index >0)%>%group_by(Class)%>%summarise(table(Class)%>%as.data.frame())%>%as.data.frame()
    a <-filter(aa,Class=="hotspot")$Freq #jaacard_index >0
    r_c<-filter(aa,Class=="random")$Freq #jaacard_index >0
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
    # all_fisher <- bind_rows(all_fisher,tmp)
    print(state)
}
res  <-lapply(marks,fisher_greater)
result_g <- do.call(rbind,res)
write.table(result_g ,"./figure/0_0.176/compare_0_0.176_jaacard_index_fisher_test_histone_tfbs_greater.txt",row.names = F, col.names = T,quote =F,sep="\t")

