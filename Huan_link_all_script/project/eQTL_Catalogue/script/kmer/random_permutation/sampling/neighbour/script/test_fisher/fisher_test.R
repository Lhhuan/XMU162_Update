library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tibble)
library(R.utils)

system("bedtools fisher -a b.bed -b a.bed -g t.genome >ab_result.txt")
org <- read.table("ab_result.txt",header = T,sep = "\t") %>% as.data.frame()

