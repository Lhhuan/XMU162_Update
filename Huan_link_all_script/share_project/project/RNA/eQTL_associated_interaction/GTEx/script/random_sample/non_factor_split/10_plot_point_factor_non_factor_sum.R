library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(parallel)


setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/")
org<-read.table("081_sum_random_result.txt",header = T,sep = "\t") %>% as.data.frame()
library(Hmisc)
factors <- c("promoter","enhancer","TFBS","CHROMATIN_Accessibility","HISTONE_modification")
tissues = c("Lung","Whole_Blood")
for(tissue in tissues){
    figure_list<-list()
    i=1
    for(factor in factors){
        factor_org <-dplyr::filter(org,Tissue==tissue & Factor == factor )
        factor = tolower(factor)
        factor <-capitalize(factor) 
        # pdf(paste0(factor,"_sum_factor.pdf"))
        
        factor <-str_replace(factor,"_"," ")
        print(factor)
        figure_list[[i]] <-ggplot(factor_org,mapping = aes(x = FPR, y = TPR,colour = Cutoff)) + geom_point(size = 2)+ggtitle(factor)+theme(plot.title = element_text(hjust = 0.5))
        # figure_list[[i]] <-ggplot(factor_org,mapping = aes(x = FPR, y = TPR)) + geom_point(size = 3)+ggtitle(factor)+theme(plot.title = element_text(hjust = 0.5))
        # print(p)
        # dev.off()
        i=i+1
    }
    #------------------
    pdf(paste0("ROC_cis_eQTL_in_",tissue,".pdf"),height = 7,width = 12)
    tissue <-str_replace(tissue,"_"," ")
    top_title <-paste0("Random select sum Cis eQTL in ",tissue)
    p1<-marrangeGrob(figure_list,nrow=2,ncol=3,top = top_title)  
    print(p1)
    dev.off()
}

#--------------------------




# enhancer <-dplyr::filter(org,Factor == "enhancer")
# promoter <-dplyr::filter(org,Factor == "promoter")
# pdf("enhancer_factor.pdf")
# p <-ggplot(enhancer,mapping = aes(x = FPR, y = TPR,colour = Cutoff)) + geom_point(size = 3)+ggtitle("Enhancer")+theme(plot.title = element_text(hjust = 0.5))
# print(p)
# dev.off()

# pdf("promoter_factor.pdf")
# p <-ggplot(promoter,mapping = aes(x = FPR, y = TPR,colour = Cutoff)) + geom_point(size = 3)+ggtitle("Promoter")+theme(plot.title = element_text(hjust = 0.5))
# print(p)
# dev.off()