library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(parallel)


setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_select_exclude_eQTL/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/")
org<-read.table("08_prepare_number_ROC_factor_count.txt",header = T,sep = "\t") %>% as.data.frame()
library(Hmisc)
nums <-c(1:10)
tissues = c("Lung","Whole_Blood")
factors <- c("promoter","enhancer","TFBS","CHROMATIN_Accessibility","HISTONE_modification")
for (tissue in tissues){
    setwd(paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_select_exclude_eQTL/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/figure/",tissue))
    for (num in nums){
        figure_list<-list()
        i=1
        for(factor in factors){
            factor_org <-dplyr::filter(org,Tissue==tissue & Random_number==num & Factor==factor)
            factor <-capitalize(factor) 
            # pdf(paste0(factor,"_factor.pdf"))
            figure_list[[i]] <-ggplot(factor_org,mapping = aes(x = FPR, y = TPR,colour = Cutoff)) + geom_point(size = 3)+ggtitle(factor)+theme(plot.title = element_text(hjust = 0.5))
        }
        top_title <-paste0("Random select ",num," Cis eQTL in ",tissue)
        pdf(paste0("ROC_cis_eQTL_in_",tissue,".pdf"))
        p1<-marrangeGrob(figure_list,nrow=2,ncol=3,top = top_title)  
        print(p1)
        dev.off()
    }
}



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