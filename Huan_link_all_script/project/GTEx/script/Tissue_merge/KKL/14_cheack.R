library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(VennDiagram)
setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/KKL/homer/specific/")

c5 <-read.table("/home/huanhuan/project/GTEx/script/Tissue_merge/KKL/homer/knowResult_specific.txt",header = T,sep = "\t") %>% as.data.frame()

for(i in 1:6){
    cc= filter(c5,community==i)
    file_n = paste0(i,"com_specific.txt")
    write.table(cc,file_n,col.names=T,row.names=F,quote=F,sep="\t")
}

