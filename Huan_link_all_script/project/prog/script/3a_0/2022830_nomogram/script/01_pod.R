library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)
library(cowplot)
library(ggpubr)
library(Hmisc)
library(reshape2)
setwd("/home/huanhuan/project/prog/script/3a_0/2022830_nomogram/output/")
load("/home/huanhuan/project/prog/data/04_add_age_raw_pfs_os_filter_grade_refine_ldh_b2m.Rdata")
dat <-filter(dat1,grade!="3b")
library(regplot)


