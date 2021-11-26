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
setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/random/background_original_random/0/")


#--------------
ProcessBedGz<-function(i=NULL,kmer=NULL){
    setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/random/background_original_random/0/")

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

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/random_permutation/sampling/emplambda_0/figure/wilcox_test/")
#--------------------------------------
# w_result <- mclapply(unique_hotspot_kmer_need_test[1:30],function(kmer){
rs <-data.frame()
# for(kmer in unique_hotspot_kmer_need_test[1:30]){
for(kmer in unique_hotspot_kmer_need_test[31:100]){
    # for(kmer in unique_hotspot_kmer_need_test[1]){
    # w_result <- mclapply(unique_hotspot_kmer_need_test[1:2],function(kmer){
    hotspot <-filter(hotspot_kmer_need_test_value,seq == kmer)
    hotspot_value <- hotspot$value
    #--------------random
    a <- mclapply(c(1:1000), function(i){
        ProcessBedGz(i=i,kmer= kmer)}, mc.cores = 25)
    random <-do.call(rbind,a)
    random_value <- random$value
    file_random_value <- paste0(kmer,"_random_value_1000.Rdata")

    setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/random_permutation/sampling/emplambda_0/figure/wilcox_test/")
    save(random,file = file_random_value )
    #-------------
    re <- wilcox.test(hotspot_value, random_value,alternative = "greater")
    file_wilcox = paste0(kmer,"_file_wilcox_1000.Rdata")
    save(re,file = file_wilcox )
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
    w_file = paste0(kmer,"_hotspot_emplambda_0_random_wilcox_1000.txt")
    write.table(w_re,w_file,row.names = F, col.names = T,quote =F,sep="\t")
    gc()
    # return(w_re)},
    # mc.cores=2
}

# write.table(rs,"hotspot_emplambda_0_random_kmer_wilcox_1000_1_30.txt",row.names = F, col.names = T,quote =F,sep="\t")
write.table(rs,"hotspot_emplambda_0_random_kmer_wilcox_1000_31_100.txt",row.names = F, col.names = T,quote =F,sep="\t")




# unique_hotspot_kmer_need_test111<-unique_hotspot_kmer_need_test[1:2]
# for(kmer in unique_hotspot_kmer_need_test){
# w_result <- mclapply(unique_hotspot_kmer_need_test111,function(kmer){

w_result <- mclapply(unique_hotspot_kmer_need_test[200:202],function(kmer){
# w_result <- mclapply(unique_hotspot_kmer_need_test[1:2],function(kmer){
    hotspot <-filter(hotspot_kmer_need_test_value,seq == kmer)
    hotspot_value <- hotspot$value
    #--------------random
    a <- mclapply(c(1:200), function(i){
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

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/random_permutation/sampling/emplambda_0/")

write.table(all_wilcox_test,"./figure/hotspot_emplambda_0_random_kmer_wilcox_1000.txt",row.names = F, col.names = T,quote =F,sep="\t")
print("Have a nice day!")

# gzip("./figure/hit_hospot_ratio_kmer_value.txt")
# gzfile("df1.txt.gz", "w")

# write.table(hotspot_kmer_need_test_value, gzfile("./figure/hit_hospot_ratio_kmer_value.txt.gz", "w"),row.names = F, col.names = T,quote =F,sep="\t")
# nrow(filter(seq_count,ratio >0.1))