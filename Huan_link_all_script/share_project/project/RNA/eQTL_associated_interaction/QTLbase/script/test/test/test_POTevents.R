library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

# R >= 3.1.1 is needed

if (!require("NHPoisson")) {
    install.packages("NHPoisson")
    library("NHPoisson")
}

# setting a random seed for reproducibility
set.seed(1138)


org<-read.table("test.txt",header = T,sep = "\t") %>% as.data.frame()

BarEv <- POTevents.fun(T = -log10(org$min_p), thres = 7, date = org$SNP_pos)

