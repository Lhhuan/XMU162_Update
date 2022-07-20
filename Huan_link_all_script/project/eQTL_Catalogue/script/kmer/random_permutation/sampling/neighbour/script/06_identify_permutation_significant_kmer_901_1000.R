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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/")
wilcox <-read.table("./output/figure/03_hotspot_10_neighbour_random_kmer_wilcox_twoside_fdr_sig_foldchange.txt",header= T,sep = "\t") %>% as.data.frame()
sigK <- filter(wilcox,Fold_change>=1.15 | Fold_change <=0.9)
kmer_need <-sigK

ProcessBedGz<-function(i=NULL){
    setwd("/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/random/10_fold_neighbour/")
    print(i)
    file_name <-paste0(i,"_6mers_uc_us_no_log.csv.gz")
    org<-read.csv(file_name,header = T,sep = ",") %>% as.data.frame()
    rownames(org) <- org[,1]
    Sorg <-org[,which(colnames(org) %in% sigK$seq)]
    Sorg$hotspot <- rownames(Sorg)
    org2 <-melt(Sorg,"hotspot")
    colnames(org2)[2] <-"seq"
    # sig_hotspot <-filter(org2, value>0 & seq %in% sigK$seq)
    sig_hotspot <-filter(org2, value>0)
    seq_count_in <- sig_hotspot%>%group_by(seq)%>%dplyr::summarise(count=n())%>%as.data.frame()
    # seq_count_in <-filter(seq_count,seq %in% org1$seq) 
    gc()
    # rs<-bind_rows(rs,sig_hotspot)
    print(paste(i,"end",sep=" "))
    seq_count_in$random_ID=i
    return(seq_count_in)
}

# all_s_h <-mclapply(c(1:100), ProcessBedGz, mc.cores = 20)
# all_s_h <-mclapply(c(1:1000), ProcessBedGz, mc.cores = 20)
# all_6kmers <-do.call(rbind,all_s_h)
# a <-lapply(c(1:1000), ProcessBedGz)
a <-lapply(c(901:1000), ProcessBedGz)
random_901_1000<-do.call(rbind,a)
save(random_901_1000,file="./output/06_permutation_random_901_1000.Rdata")
