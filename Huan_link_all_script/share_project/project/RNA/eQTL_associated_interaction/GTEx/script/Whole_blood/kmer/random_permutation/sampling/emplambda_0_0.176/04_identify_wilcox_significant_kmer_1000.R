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

load("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/figure/hotspot_kmer_need_test_value.Rdata")
unique_hotspot_kmer_need_test <-as.character(unique(hotspot_kmer_need_test_value$seq))
#---------
setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/random/background_original_random/0_0.176/")


#--------------
ProcessBedGz<-function(i=NULL,kmer=NULL){
    setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/random/background_original_random/0_0.176/")

    print(i)
    file_name <-paste0(i,"_6mers_uc_us_no_log.csv.gz")
    org<-read.csv(file_name,header = T,sep = ",") %>% as.data.frame()
    colnames(org)[1] <-"hotspot"
    org2 <-melt(org,"hotspot")
    colnames(org2)[2] <-"seq"
    org3 <-dplyr::select(org2,-hotspot)
    rm(org2)
    random_kmer_need_test_value <-filter(org3,seq==kmer)
    gc()
    return(random_kmer_need_test_value)
}



# unique_hotspot_kmer_need_test111<-unique_hotspot_kmer_need_test[1:2]
# for(kmer in unique_hotspot_kmer_need_test){
# w_result <- mclapply(unique_hotspot_kmer_need_test111,function(kmer){
w_result <- mclapply(unique_hotspot_kmer_need_test,function(kmer){
    hotspot <-filter(hotspot_kmer_need_test_value,seq == kmer)
    hotspot_value <- hotspot$value
    #--------------random
    a <- mclapply(c(1:1000), function(i){
        ProcessBedGz(i=i,kmer= kmer)}, mc.cores = 15)
    random <-do.call(rbind,a)
    random_value <- random$value
    #-------------
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
    gc()
    return(w_re)},
    mc.cores=2
)



all_wilcox_test <-do.call(rbind,w_result)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/random_permutation/sampling/emplambda_0_0.176/")

write.table(all_wilcox_test,"./figure/hotspot_emplambda_0_0.176_random_kmer_wilcox_1000.txt",row.names = F, col.names = T,quote =F,sep="\t")
print("Have a nice day!")

# gzip("./figure/hit_hospot_ratio_kmer_value.txt")
# gzfile("df1.txt.gz", "w")

# write.table(hotspot_kmer_need_test_value, gzfile("./figure/hit_hospot_ratio_kmer_value.txt.gz", "w"),row.names = F, col.names = T,quote =F,sep="\t")
# nrow(filter(seq_count,ratio >0.1))