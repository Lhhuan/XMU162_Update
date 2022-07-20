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
library(parallel)

# setwd("/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/hotspot/")
org<-read.csv("/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/hotspot/6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()
colnames(org)[1] <-"hotspot"
org2 <-melt(org,"hotspot")
colnames(org2)[2] <-"seq"

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/")
all_random<-read.csv("6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()
colnames(all_random)[1] <-"hotspot"
all_random2 <-melt(all_random,"hotspot")
colnames(all_random2)[2] <-"seq"

unique_hotspot_kmer_need_test<-as.character(unique(org2$seq))
rs <-data.frame()
for(kmer in unique_hotspot_kmer_need_test){
    # ProcessBedGz_test <-function(hotspot_kmer_need_test_value =NULL, all_random_kmer_need_test_value =NULL,kmer= NULL){
    hotspot <-filter(org2,seq == kmer)
    hotspot_value <- hotspot$value
    random <-filter(all_random2 ,seq == kmer)
    random_value <- random$value
    re <- wilcox.test(hotspot_value, random_value,alternative = "greater")
    p_value = re$p.value
    #---------------------
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
    #-------------
    W <-re$statistic
    names(W)=NULL
    w_re <-data.frame(seq = kmer,W=W,p_value = re$p.value,annotation=annotation)
    rs <-bind_rows(rs,w_re)
    print(kmer)
    # gc()
    # return(w_re)
}
#------------------------

write.table(rs,"./figure/hotspot_10_neighbour_random_kmer_wilcox.txt",row.names = F, col.names = T,quote =F,sep="\t")


# for(kmer in unique_hotspot_kmer_need_test){
ProcessBedGz_test <-function(kmer= NULL){
    hotspot <-filter(org2,seq == kmer)
    hotspot_value <- hotspot$value
    random <-filter(all_random2 ,seq == kmer)
    random_value <- random$value
    re <- wilcox.test(hotspot_value, random_value,alternative = "greater")
    p_value = re$p.value
    #---------------------
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
    #-------------
    W <-re$statistic
    names(W)=NULL
    w_re <-data.frame(seq = kmer,W=W,p_value = re$p.value,annotation=annotation)
    rs <-bind_rows(rs,w_re)
    print(kmer)
    # gc()
    # return(w_re)
    return(rs)
}

AA <-mclapply(unique_hotspot_kmer_need_test,ProcessBedGz_test,mc.cores = 20)