library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(Seurat)

# cancers <-c("BRCA","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","STAD","THCA")
cancers <-c("BRCA","KIRC","KIRP","LAML","LIHC","LUAD","OV","PAAD","PRAD","STAD","THCA")

c_p <-function(cancer =NULL){
# for(cancer in cancers){
    print(cancer)
    path = paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/",cancer,"/Cis_eQTL/anno/")
    setwd(path)
    print(path)
    cancer_file <- paste0(cancer,"_segment_hotspot_cutoff_0.176_marker_jaccard_index_label.txt.gz")
    org <- read.table(cancer_file,header = T,sep = "\t") %>% as.data.frame()
    tri_markers <-c("CHROMATIN_Accessibility","TFBS","H3K27ac","H3K9ac","H3K36me3","H3K4me1","H3K4me3","H3K27me3","H3K9me3")
    unique_marker <-unique(org$marker)
    overlap_marker <- unique_marker[unique_marker%in%tri_markers]
    # table_length = nrow(org)/length(unique_marker)
    fisher_t <-function(state=NULL){
        bb <-filter(org,marker==state)
        class_number <-table(bb$Class)
        all_specific_n <-as.numeric(class_number[2])
        all_other_n <-as.numeric(class_number[1])
        # aa <-bb%>%filter(jaccard_index >0)%>%group_by(Class)%>%summarise(table(Class)%>%as.data.frame())%>%as.data.frame()
        cccccc <- bb%>%filter(jaccard_index >0)%>%group_by(Class)%>%as.data.frame()
        bbbb <-table(cccccc$Class)%>%as.data.frame()
        colnames(bbbb)[1] <-"Class"
        a <-filter(bbbb,Class=="Specific")$Freq
        r_c<-filter(bbbb,Class=="Others")$Freq 
        pre_data <-matrix(c(a,r_c, all_specific_n-a, all_other_n-r_c), nrow = 2)
        #------------------------------------
        fisher_test_result <-fisher.test(pre_data)

        p_value =fisher_test_result$p.value
        conf_int1 <-fisher_test_result$conf.int[1]
        conf_int2 <-fisher_test_result$conf.int[2]
        od <-fisher_test_result$estimate
        names(od)=NULL

        # if(p_value<2.2e-16){
        #     p_value ="< 2.2e-16"
        # }else{
        #     p_value=p_value
        # }
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
    }
    #----------------------
    rs <-lapply(overlap_marker,fisher_t)
    fisher_test_result <-do.call(rbind,rs)
    write.table(fisher_test_result,"29_share_specific_fisher_test_two_side.txt",row.names = F, col.names = T,quote =F,sep="\t")
}

lapply(cancers,c_p)
