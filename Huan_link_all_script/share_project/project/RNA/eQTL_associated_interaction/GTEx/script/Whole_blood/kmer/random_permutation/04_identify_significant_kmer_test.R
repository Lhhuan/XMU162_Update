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


# load("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/figure/hotspot_kmer_need_test_value.Rdata")

setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/")
org<-read.csv("6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()
colnames(org)[1] <-"hotspot"
org2 <-melt(org,"hotspot")
colnames(org2)[2] <-"seq"

sig_hotspot <-filter(org2, value>0)
seq_count <- sig_hotspot%>%group_by(seq)%>%summarise(count=n())%>%as.data.frame()
seq_count$ratio <-seq_count$count/nrow(org)

rs <-data.frame()

cutoffs<-seq(0,0.25,0.01)

for(cutoff in cutoffs){
    aa <-filter(seq_count,ratio >cutoff)
    kmer_ratio =nrow(aa)/(4^6)
    kmer_count = nrow(aa)
    tmp <-data.frame(hotspot_ratio =cutoff,kmer_ratio = kmer_ratio,kmer_count=kmer_count)
    rs<-bind_rows(rs,tmp)
    print(cutoff)
}
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/")
write.table(rs,"./figure/hit_hospot_ratio_kmer_count.txt",row.names = F, col.names = T,quote =F,sep="\t")
write.table(seq_count,"./figure/hit_hospot_ratio_kmer.txt",row.names = F, col.names = T,quote =F,sep="\t")
#------
hotspot_kmer_need_test <-filter(seq_count,ratio >0.1)
# hotspot_kmer_need_test <-hotspot_kmer_need_test

unique_hotspot_kmer_need_test <-as.data.frame(unique(hotspot_kmer_need_test$seq))
colnames(unique_hotspot_kmer_need_test)[1] <-"seq"
org3 <-dplyr::select(org2,-hotspot)
hotspot_kmer_need_test_value <-inner_join(org3,unique_hotspot_kmer_need_test,by="seq")

save(hotspot_kmer_need_test_value,file="./figure/hotspot_kmer_need_test_value.Rdata")
write.table(hotspot_kmer_need_test_value,"./figure/hit_hospot_ratio_kmer_value.txt",row.names = F, col.names = T,quote =F,sep="\t")

gzip("./figure/hit_hospot_ratio_kmer_value.txt")
# gzfile("df1.txt.gz", "w")

# write.table(hotspot_kmer_need_test_value, gzfile("./figure/hit_hospot_ratio_kmer_value.txt.gz", "w"),row.names = F, col.names = T,quote =F,sep="\t")
# nrow(filter(seq_count,ratio >0.1))