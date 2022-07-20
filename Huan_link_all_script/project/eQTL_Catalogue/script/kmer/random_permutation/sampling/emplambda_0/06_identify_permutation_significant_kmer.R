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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/emplambda_0/")
org1 <- read.table("../../../figure/hit_hospot_ratio_kmer.txt",header = T,sep = "\t") %>% as.data.frame()
kmer_need <-filter(org1,ratio>0.5)


ProcessBedGz<-function(i=NULL){
    setwd("/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/random/emp0/")
    print(i)
    file_name <-paste0(i,"_6mers_uc_us_no_log.csv.gz")
    org<-read.csv(file_name,header = T,sep = ",") %>% as.data.frame()
    colnames(org)[1] <-"hotspot"
    org2 <-melt(org,"hotspot")
    colnames(org2)[2] <-"seq"
    sig_hotspot <-filter(org2, value>0)
    seq_count <- sig_hotspot%>%group_by(seq)%>%dplyr::summarise(count=n())%>%as.data.frame()
    seq_count_in <-filter(seq_count,seq %in% kmer_need$seq) 
    gc()
    # rs<-bind_rows(rs,sig_hotspot)
    # print(i)
    seq_count_in$random_ID=i
    return(seq_count_in)
}

# all_s_h <-mclapply(c(1:10), ProcessBedGz, mc.cores = 20)
# all_s_h <-mclapply(c(1:1000), ProcessBedGz, mc.cores = 20)
# all_6kmers <-do.call(rbind,all_s_h)
a <-lapply(c(1:100), ProcessBedGz)



all_random_kmer_need_test_value <-do.call(rbind,a)
all_random_kmer_need_test_ratio <-all_random_kmer_need_test_value
save(all_random_kmer_need_test_ratio,file ="all_random_kmer_need_test_ratio_1_100.Rdata") 
hotf<-read.table("/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t")%>% as.data.frame()

random_more0_count <- all_random_kmer_need_test_ratio
random_more0_count$ratio <-random_more0_count$count/nrow(hotf)
random <-filter(random_more0_count,seq %in%kmer_need$seq)
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/emplambda_0/")
save(random,file="all_random_kmer_need_test_ratio_1_100.Rdata")

colnames(kmer_need)[2:3] <-c("hotspot_count","hotspot_ratio")
dat <-inner_join(random,kmer_need,by="seq")
dat$diff <- dat$ratio -dat$hotspot_ratio 
dat$diff0 <-NA
dat[which(dat$diff>=0),"diff0"]=1
dat[which(dat$diff<0),"diff0"]=0
fdat <-dat%>%group_by(seq)%>%summarise(n=sum(diff0))%>%data.frame()
fdat$p_value <-fdat$n/100

fdat0 <-filter(fdat,p_value<0.05)
save(fdat0,file="06_permutation_test_100_sig_kmer.Rdata")
