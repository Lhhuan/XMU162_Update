library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(VennDiagram)
setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/chr1_6/6kmer/5_community/specific")

c5 <-read.table("/home/huanhuan/project/GTEx/script/Tissue_merge/chr1_6/6kmer/5_community/knowResult_specific.txt",header = T,sep = "\t") %>% as.data.frame()
# c5$keys = paste(c5$V1,c5$V2,c5$V3,sep=":")
# load("/home/huanhuan/project/GTEx/script/Tissue_merge/chr1_6/6kmer/5_community/bind_tf5.Rdata")
# tf$key=paste(tf$Motif_Name,tf$Consensus,tf$community)
for(i in 1:5){
    cc= filter(c5,community==i)
    file_n = paste0(i,"com_specific.txt")
    write.table(cc,file_n,col.names=T,row.names=F,quote=F,sep="\t")
}



# read_html("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_6/kmer/6/5_community/homer/1/knownResults.html")
# path = "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_6/kmer/6/5_community/homer/1/knownResults.html"

# path %>% 
#   read_html() %>% 
#   html_element(css = css_selector_paragraph) %>% 
#   html_text()

