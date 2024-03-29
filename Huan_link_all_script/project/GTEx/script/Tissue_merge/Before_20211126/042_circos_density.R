
# source("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/script/OmicCircos_circos_refine_revise_color.R")
library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(tidyverse)
library(reshape2)
library(Seurat)
library(circlize)

setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/")


tissue_merge <- read.table("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_sorted.bed.gz",header = F,sep = "\t") %>% as.data.frame()
#-----


#--------------------------
circos.clear()
pdf("./figure/042_circos_hotspot_density_all_tissue_0.2.pdf")
par(mar = c(2, 2, 2, 2))
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species= "hg19",chromosome.index = paste0("chr", 1:22))
# circos.genomicDensity(tissue_merge, col = c("#00EAD3"), track.height = 0.1,ylim.force=T,window.size = 1e4)
circos.genomicDensity(tissue_merge, col = c("#F5A962"), track.height = 0.2,ylim.force=T,window.size = 1e4)

dev.off()
print("finish")

#----------------------------
circos.clear()
pdf("./figure/042_circos_hotspot_density_all_tissue_0.3.pdf")
par(mar = c(2, 2, 2, 2))
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species= "hg19",chromosome.index = paste0("chr", 1:22))
# circos.genomicDensity(tissue_merge, col = c("#00EAD3"), track.height = 0.1,ylim.force=T,window.size = 1e4)
circos.genomicDensity(tissue_merge, col = c("#F5A962"), track.height = 0.3,ylim.force=T,window.size = 1e4)

dev.off()
print("finish")

circos.clear()
pdf("./figure/042_circos_hotspot_density_all_tissue_0.4.pdf")
par(mar = c(2, 2, 2, 2))
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species= "hg19",chromosome.index = paste0("chr", 1:22))
# circos.genomicDensity(tissue_merge, col = c("#00EAD3"), track.height = 0.1,ylim.force=T,window.size = 1e4)
circos.genomicDensity(tissue_merge, col = c("#F5A962"), track.height = 0.4,ylim.force=T,window.size = 1e4)

dev.off()
print("finish")
#---------------------------
# Breast <-read.table("../../output/Breast_Mammary_Tissue/Cis_eQTL/hotspot_cis_eQTL/interval_18/Breast_Mammary_Tissue_segment_hotspot_cutoff_0.176.bed.gz",header = F,sep = "\t") %>% as.data.frame()
# Liver <-read.table("../../output/Liver/Cis_eQTL/hotspot_cis_eQTL/interval_18/Liver_segment_hotspot_cutoff_0.176.bed.gz",header = F,sep = "\t") %>% as.data.frame()
# Lung <-read.table("../../output/Lung/Cis_eQTL/hotspot_cis_eQTL/interval_18/Lung_segment_hotspot_cutoff_0.176.bed.gz",header = F,sep = "\t") %>% as.data.frame()
# Ovary <-read.table("../../output/Ovary/Cis_eQTL/hotspot_cis_eQTL/interval_18/Ovary_segment_hotspot_cutoff_0.176.bed.gz",header = F,sep = "\t") %>% as.data.frame()
# Pancreas <-read.table("../../output/Pancreas/Cis_eQTL/hotspot_cis_eQTL/interval_18/Pancreas_segment_hotspot_cutoff_0.176.bed.gz",header = F,sep = "\t") %>% as.data.frame()
# Prostate <-read.table("../../output/Prostate/Cis_eQTL/hotspot_cis_eQTL/interval_18/Prostate_segment_hotspot_cutoff_0.176.bed.gz",header = F,sep = "\t") %>% as.data.frame()
# Thyroid <-read.table("../../output/Thyroid/Cis_eQTL/hotspot_cis_eQTL/interval_18/Thyroid_segment_hotspot_cutoff_0.176.bed.gz",header = F,sep = "\t") %>% as.data.frame()
# Whole_Blood <-read.table("../../output/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18/whole_blood_segment_hotspot_cutoff_0.176.bed.gz",header = F,sep = "\t") %>% as.data.frame()


#----------------------
# circos.genomicDensity(Breast, col = c("#B6C867"), track.height = 0.1,ylim.force=T,window.size = 1e4)
# # circos.genomicDensity(Liver, col = c("#FFC074"), track.height = 0.1)
# circos.genomicDensity(Lung, col = c("#01937C"), track.height = 0.1,ylim.force=T,window.size = 1e4)
# # circos.genomicDensity(Ovary, col = c("#125D98"), track.height = 0.1)
# circos.genomicDensity(Pancreas, col = c("#3C8DAD"), track.height = 0.1,ylim.force=T,window.size = 1e4)
# # circos.genomicDensity(Prostate, col = c("#F5A962"), track.height = 0.1)
# circos.genomicDensity(Thyroid, col = c("#2541B2"), track.height = 0.1,ylim.force=T,window.size = 1e4)
# circos.genomicDensity(Whole_Blood, col = c("#C84B31"), track.height = 0.1,ylim.force=T,window.size = 1e4)