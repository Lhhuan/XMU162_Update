library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(pheatmap)
library(reshape2)
library(tidyverse)

setwd("/home/huanhuan/project/link_database/COXPRESdb/")
# All <-read.table("Hsa-r.c5-0.expression.combat.txt",header = F,sep = "\t") %>% as.data.frame()
All <-read_tsv("Hsa-r.c5-0.expression.combat.txt",col_names = F) %>% as.data.frame()


