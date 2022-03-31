library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(Seurat)

setwd("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/enrichment/")
hotspot<-read.table("hotspot_cutoff_0.176_marker_jaccard_index_filter_3103833.txt.gz",header = T,sep = "\t") %>% as.data.frame()
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
all_res <-data.frame()
for(state in marks){
    aa <-filter(rs,Marker==state)%>%filter(jaacard_index >0)%>%group_by(Class)%>%summarise(table(Class)%>%as.data.frame())%>%as.data.frame()
    a <-filter(aa,Class=="hotspot")$Freq
    r_c<-filter(aa,Class=="random")$Freq
    pre_data <-matrix(c(a,r_c, nrow(hotspot)/10-a, nrow(hotspot)/10*1000-r_c), nrow = 2)

    # fisher_test_result <-fisher.test(pre_data,alternative = "greater")
    # fisher_test_result <-fisher.test(pre_data)
    tab_Xsqtest <- chisq.test(pre_data)
    x=tab_Xsqtest$statistic
    names(x)<-NULL
    p_value=round(tab_Xsqtest$p.value,4)
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
    tmp <-data.frame(marker =state,p_value=p_value,Chi=round(x,2),significant=annotation)
    all_res <- bind_rows(all_res,tmp)
    print(state)
}

write.table(all_res,"./figure/0_0.176/compare_0_0.176_jaacard_index_histone_tfbs_chisq_test.txt",row.names = F, col.names = T,quote =F,sep="\t")

