library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)
library(cowplot)
library(ggpubr)
library(Hmisc)
library("survival")
library("survminer")

setwd("/home/huanhuan/project/prog/output/")
load("final_3a_data20220512.Rdata")

dat <-dat1 
library(magrittr)
pdf("10_4_3a_all_calibration_plot.pdf")
p <- calibration_plot(data = dat, obs = "new_pod_total", pred = "Ours", title = "Calibration plot for development data")
p
dev.off()
test <-dat[test_set_number,]
