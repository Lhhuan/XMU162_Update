library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(parallel)
library(conflicted)






setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/output/ALL_eQTL/cis_trans/hotspot")

intervals=c(12,15,18)

types= c("cis_1MB","cis_10MB")
for(interval in intervals){
  for(type in types){
    org<-read.table(paste0("interval_",interval,"_cutoff_7.3_",type,"_eQTL_segment_hotspot_length_more_than6.bed.gz"),header = F,sep = "\t") %>% as.data.frame()
    rs<-data.frame()
    for(i in c(1:22)){
      chr <-paste0("chr",i)
      org_chr<-dplyr::filter(org,V1==chr)
      set.seed(100)
      a<-round(nrow(org_chr)*0.1)
      random_number <-sample(x=c(1:nrow(org_chr)), a,replace = F)
      org_a<-org_chr[random_number,]
      rs<-bind_rows(rs,org_a)
    }
      # file <-paste("/random_select/interval_18_cutoff_7.3_cis_10MB_eQTL_segment_hotspot_length_more_than6_random.bed.gz")
      # gz1 <- gzfile(paste0("interval_",interval,"_cutoff_7.3_",type,"_eQTL_segment_hotspot_length_more_than6_random_select_fraction_0.1.bed.gz"), "w")
      # write.table(rs,gz1,row.names = F, col.names = F,quote =F,sep="\t")
      # close(gz1)
      write.table(rs, gzfile(paste0("random_select/interval_",interval,"_cutoff_7.3_",type,"_eQTL_segment_hotspot_length_more_than6_random_select_fraction_0.1.bed.gz")),
        r = F, c = F, quote =F,sep="\t")
  }
}

# write.table(rs, gzfile(paste0("interval_",interval,"_cutoff_7.3_",type,"_eQTL_segment_hotspot_length_more_than6_random_select_fraction_0.1.bed.gz")),
        # r = F, c = F, quote =F,sep="\t")