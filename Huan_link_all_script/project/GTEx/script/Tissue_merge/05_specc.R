library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tibble)
library(Seurat)
library(R.utils)
library(reshape2)
library(parallel)


p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))

#---------------------------------------------------------
library(Hmisc)
setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/figure/")

load("04_hotspot_mark_fraction.Rdata")

# save(dat,file="04_hotspot_mark_fraction.Rdata")

data  <- dat[,-1]

#------
#----------------------------------specc
library("kernlab")
dat1 <-as.matrix(data)
dat2 <-head(dat1,1000)
sc <-specc(dat1,centers=10)




load("04_10_markers_euclidean_dist.Rdata")

dis <- specc(eu_dist,centers=10)

specc_n <-function(i){
    print(i)
    sc <-specc(data,centers=i)
}
re <-lapply(c(2:10),specc_n)
re <-mclapply(c(2:10),specc_n,mc.cores=9)

#---------------
