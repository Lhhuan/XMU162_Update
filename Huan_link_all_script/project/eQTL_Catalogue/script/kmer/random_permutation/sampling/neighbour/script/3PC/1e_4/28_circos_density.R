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
library(circlize)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/")
org<-read.table("./GWAS/17_whole_genome_leiden_pca3_k50_resolution1e-04_sorted.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(org)<-c("CHR","Start","End","Cluster")
tmp <-lapply(1:6,function(i){
    cl <- filter(org,Cluster==i)
    return(cl)
})

#--------------------------
circos.clear()
pdf("28_circos_hotspot_density_6cluster_10kb.pdf")
par(mar = c(2, 2, 2, 2))
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species= "hg38",chromosome.index = paste0("chr", 1:22))
circos.genomicDensity(tmp[[1]], col = c("#1E77B4"), track.height = 0.1,ylim.force=T,window.size = 1e4)
circos.genomicDensity(tmp[[2]], col = c("#FF7F0E"), track.height = 0.1,ylim.force=T,window.size = 1e4)
circos.genomicDensity(tmp[[3]], col = c("#2CA02C"), track.height = 0.1,ylim.force=T,window.size = 1e4)
circos.genomicDensity(tmp[[4]], col = c("#C22324"), track.height = 0.1,ylim.force=T,window.size = 1e4)
circos.genomicDensity(tmp[[5]], col = c("#9567BD"), track.height = 0.1,ylim.force=T,window.size = 1e4)
circos.genomicDensity(tmp[[6]], col = c("#8C554B"), track.height = 0.1,ylim.force=T,window.size = 1e4)

dev.off()
print("finish")
#========================
circos.clear()
pdf("28_circos_hotspot_density_6cluster_1e7bp.pdf")
par(mar = c(2, 2, 2, 2))
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species= "hg38",chromosome.index = paste0("chr", 1:22))
circos.genomicDensity(tmp[[1]], col = c("#1E77B4"), track.height = 0.1,ylim.force=T,window.size = 1e7)
circos.genomicDensity(tmp[[2]], col = c("#FF7F0E"), track.height = 0.1,ylim.force=T,window.size = 1e7)
circos.genomicDensity(tmp[[3]], col = c("#2CA02C"), track.height = 0.1,ylim.force=T,window.size = 1e7)
circos.genomicDensity(tmp[[4]], col = c("#C22324"), track.height = 0.1,ylim.force=T,window.size = 1e7)
circos.genomicDensity(tmp[[5]], col = c("#9567BD"), track.height = 0.1,ylim.force=T,window.size = 1e7)
circos.genomicDensity(tmp[[6]], col = c("#8C554B"), track.height = 0.1,ylim.force=T,window.size = 1e7)

dev.off()
print("finish")
#==========================

#========================
circos.clear()
pdf("28_circos_hotspot_density_6cluster_1e6bp.pdf")
par(mar = c(2, 2, 2, 2))
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species= "hg38",chromosome.index = paste0("chr", 1:22))
circos.genomicDensity(tmp[[1]], col = c("#1E77B4"), track.height = 0.1,window.size = 1e6)
circos.genomicDensity(tmp[[2]], col = c("#FF7F0E"), track.height = 0.1,window.size = 1e6)
circos.genomicDensity(tmp[[3]], col = c("#2CA02C"), track.height = 0.1,window.size = 1e6)
circos.genomicDensity(tmp[[4]], col = c("#C22324"), track.height = 0.1,window.size = 1e6)
circos.genomicDensity(tmp[[5]], col = c("#9567BD"), track.height = 0.1,window.size = 1e6)
circos.genomicDensity(tmp[[6]], col = c("#8C554B"), track.height = 0.1,window.size = 1e6)

dev.off()
print("finish")
#========================
circos.clear()
pdf("28_circos_hotspot_density_6cluster_1e5bp.pdf")
par(mar = c(2, 2, 2, 2))
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species= "hg38",chromosome.index = paste0("chr", 1:22))
circos.genomicDensity(tmp[[1]], col = c("#1E77B4"), track.height = 0.1,window.size = 1e5)
circos.genomicDensity(tmp[[2]], col = c("#FF7F0E"), track.height = 0.1,window.size = 1e5)
circos.genomicDensity(tmp[[3]], col = c("#2CA02C"), track.height = 0.1,window.size = 1e5)
circos.genomicDensity(tmp[[4]], col = c("#C22324"), track.height = 0.1,window.size = 1e5)
circos.genomicDensity(tmp[[5]], col = c("#9567BD"), track.height = 0.1,window.size = 1e5)
circos.genomicDensity(tmp[[6]], col = c("#8C554B"), track.height = 0.1,window.size = 1e5)

dev.off()
print("finish")
#==========================
#========================
circos.clear()
pdf("28_circos_hotspot_density_6cluster_1e4bp.pdf")
par(mar = c(2, 2, 2, 2))
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species= "hg38",chromosome.index = paste0("chr", 1:22))
circos.genomicDensity(tmp[[1]], col = c("#1E77B4"), track.height = 0.1,window.size = 1e4)
circos.genomicDensity(tmp[[2]], col = c("#FF7F0E"), track.height = 0.1,window.size = 1e4)
circos.genomicDensity(tmp[[3]], col = c("#2CA02C"), track.height = 0.1,window.size = 1e4)
circos.genomicDensity(tmp[[4]], col = c("#C22324"), track.height = 0.1,window.size = 1e4)
circos.genomicDensity(tmp[[5]], col = c("#9567BD"), track.height = 0.1,window.size = 1e4)
circos.genomicDensity(tmp[[6]], col = c("#8C554B"), track.height = 0.1,window.size = 1e4)

dev.off()
print("finish")
#==========================

#=============================

       --colors "#1E77B4" "#FF7F0E" "#2CA02C" "#C22324" "#9567BD" "#8C554B"\
#----------------------------
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