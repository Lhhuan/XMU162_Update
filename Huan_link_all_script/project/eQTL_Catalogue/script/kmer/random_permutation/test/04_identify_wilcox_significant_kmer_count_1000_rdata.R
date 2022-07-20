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

unique_hotspot_kmer_need_test <-as.data.frame(unique(hotspot_kmer_need_test_value$seq))
colnames(unique_hotspot_kmer_need_test)[1] <-"seq"
#---------
setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/random/original_random/")

#-------------
ProcessBedGz<-function(i=NULL){
    setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/random/original_random/")
    print(i)
    file_name <-paste0(i,"_6mers_uc_us_no_log.csv.gz")
    org<-read.csv(file_name,header = T,sep = ",") %>% as.data.frame()
    colnames(org)[1] <-"hotspot"
    org2 <-melt(org,"hotspot")
    colnames(org2)[2] <-"seq"
    org3 <-dplyr::select(org2,-hotspot)
    rm(org2)
    random_kmer_need_test_value <-inner_join(org3,unique_hotspot_kmer_need_test,by="seq")
    gc()
    return(random_kmer_need_test_value)
}

# a <-mclapply(c(1:100), ProcessBedGz, mc.cores = 15)
# all_random_kmer_need_test_value <-do.call(rbind,a)
# setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/random_permutation/figure/")
# save(all_random_kmer_need_test_value,file ="all_random_kmer_need_test_value_1_100.Rdata") 
# #---------------------
# unique_hotspot_kmer_need_test1 <-as.character(unique(unique_hotspot_kmer_need_test$seq))
# rs <-data.frame()
# for(kmer in unique_hotspot_kmer_need_test1){
#     # ProcessBedGz_test <-function(hotspot_kmer_need_test_value =NULL, all_random_kmer_need_test_value =NULL,kmer= NULL){
#     hotspot <-filter(hotspot_kmer_need_test_value,seq == kmer)
#     hotspot_value <- hotspot$value
#     random <-filter(all_random_kmer_need_test_value ,seq == kmer)
#     random_value <- random$value
#     re <- wilcox.test(hotspot_value, random_value,alternative = "greater")
#     p_value = re$p.value
#     #---------------------
#     if(p_value <= 0.0001){
#         annotation = "****"
#     }else if(p_value <= 0.001){
#         annotation = "***"
#     }else if(p_value <= 0.01){
#         annotation = "**"
#     }else if(p_value <= 0.05){
#         annotation = "*"
#     }else{
#         annotation = "ns"
#     }
#     #-------------
#     W <-re$statistic
#     names(W)=NULL
#     w_re <-data.frame(seq = kmer,W=W,p_value = re$p.value,annotation=annotation)
#     rs <-bind_rows(rs,w_re)
#     print(kmer)
#     # gc()
#     # return(w_re)
# }
# #-------------------

# setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/random_permutation/")

# write.table(rs,"./figure/hotspot_original_random_kmer_wilcox.txt",row.names = F, col.names = T,quote =F,sep="\t")



# --------------------
# save(all_random_kmer_need_test_value,file ="./figure/all_random_kmer_need_test_value_1_100.Rdata") 
# --------------
number_list <-list()
j=0
for (i in seq(0,1000,99)){
    if(tail(seq(0,1000,99),1)==i){
            ind1 = i+1
            ind2 = 1000
         #a = mclapply((i+1):1000, ProcessBedGz, mc.cores = 20 )
    }else{
         ind1 = i+1
         ind2 = i+99
    }
    a = mclapply(ind1:ind2, ProcessBedGz, mc.cores = 15 )
    b <-do.call(rbind,a)
    setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/random_permutation/figure/")
    save(b,file =paste0("all_random_kmer_need_test_value","_",ind1,"_",ind2,".Rdata") )
    j=j+1
    number = paste0(ind1,"_",ind2)
    number_list[[j]] <- number
}
#-------------

# # #---------

# # setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/random_permutation/figure/")
# # all_hotspot_kmer_need_test_value <-data.frame()
# # # kmer_random_List <-list()
# # for(i in c(1:length(number_list))){
# #     print(i)
# #     file_names =paste0("all_random_kmer_need_test_value","_",number_list[[i]],".Rdata")
# #     load(file_names)
# #     all_hotspot_kmer_need_test_value <-bind_rows(all_hotspot_kmer_need_test_value,b)
# #     # aaa <-b
# #     # kmer_random_List[[i]] <-aaa
# #     gc()
# # }



# # # load("all_random_kmer_need_test_value_991_1000.Rdata")

# # # load("all_random_kmer_need_test_value_1_99.Rdata")
# # # aaaaaaaa <-load("all_random_kmer_need_test_value_991_1000.Rdata")
# # # kmer_random_List[[1]] <-load("all_random_kmer_need_test_value_991_1000.Rdata")


# # all_hotspot_kmer_need_test_value <-do.call(rbind,kmer_random_List)

# unique_hotspot_kmer_need_test1 <-as.character(unique(unique_hotspot_kmer_need_test$seq))
# rs <-data.frame()
# for(kmer in unique_hotspot_kmer_need_test1){
#     # ProcessBedGz_test <-function(hotspot_kmer_need_test_value =NULL, all_random_kmer_need_test_value =NULL,kmer= NULL){
#     hotspot <-filter(hotspot_kmer_need_test_value,seq == kmer)
#     hotspot_value <- hotspot$value
#     random <-filter(all_random_kmer_need_test_value ,seq == kmer)
#     random_value <- random$value
#     re <- wilcox.test(hotspot_value, random_value,alternative = "greater")
#     p_value = re$p.value
#     #---------------------
#     if(p_value <= 0.0001){
#         annotation = "****"
#     }else if(p_value <= 0.001){
#         annotation = "***"
#     }else if(p_value <= 0.01){
#         annotation = "**"
#     }else if(p_value <= 0.05){
#         annotation = "*"
#     }else{
#         annotation = "ns"
#     }
#     #-------------
#     W <-re$statistic
#     names(W)=NULL
#     w_re <-data.frame(seq = kmer,W=W,p_value = re$p.value,annotation=annotation)
#     rs <-bind_rows(rs,w_re)
#     print(kmer)
#     # gc()
#     # return(w_re)
# }

# w_result <- mclapply(unique_hotspot_kmer_need_test1, function(kmer){
#         ProcessBedGz_test(hotspot_kmer_need_test_value =hotspot_kmer_need_test_value, all_random_kmer_need_test_value =all_random_kmer_need_test_value,kmer= kmer)}, mc.cores = 10)

# all_wilcox_test <-do.call(rbind,w_result)





# setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/random_permutation/")

# write.table(all_wilcox_test,"./figure/hotspot_original_random_kmer_wilcox.txt",row.names = F, col.names = T,quote =F,sep="\t")

# # gzip("./figure/hit_hospot_ratio_kmer_value.txt")
# # gzfile("df1.txt.gz", "w")

# # write.table(hotspot_kmer_need_test_value, gzfile("./figure/hit_hospot_ratio_kmer_value.txt.gz", "w"),row.names = F, col.names = T,quote =F,sep="\t")
# # nrow(filter(seq_count,ratio >0.1))