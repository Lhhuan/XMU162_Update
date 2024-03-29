# library(circlize)

# load(paste(system.file(package = "circlize"), "/extdata/DMR.RData", sep=""))
# circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22))
# circos.genomicDensity(DMR_hyper, col = c("#FF000080"), track.height = 0.1)
# circos.genomicDensity(DMR_hyper, col = c("#FF000080"), count_by = "number", track.height = 0.1)



#---------------------
library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(Seurat)
library(circlize)
cancers <-c("BRCA","KIRP")

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/cancer_total/trans/")
    #--------BRCA
# cancer = "BRCA"
for (cancer in cancers){
    trans_file <- paste0("../../",cancer,"/Trans_eQTL/hotspot_trans_eQTL/interval_18/",cancer,"_segment_hotspot_cutoff_0.176.bed.gz")
    cis_file <- paste0("../../",cancer,"/Cis_eQTL/hotspot_cis_eQTL/interval_18/",cancer,"_segment_hotspot_cutoff_0.176.bed.gz")
    tran_gene_file <- paste0("../../",cancer,"/Trans_eQTL/hotspot_trans_eQTL/interval_18/hotspot_gene.bed.gz")

    cis <- read.table(cis_file,header = F,sep = "\t") %>% as.data.frame()
    trans <- read.table(trans_file,header = F,sep = "\t") %>% as.data.frame()
    trans_gene <-read.table(tran_gene_file,header = F,sep = "\t") %>% as.data.frame()
    gene_pos <-read.table("/home/huanhuan/ref_data/Ensembl/ensemble_v104_hg19_gene_pos.txt",header=T,sep="\t")%>%as.data.frame()
    gene_pos <-gene_pos%>%dplyr::select(-c(1:2))
    colnames(gene_pos) <-c("gene","gene_chr","gene_start","gene_end")
    colnames(cis) <-c("chr","start","end")
    colnames(trans) <-c("chr","start","end")
    colnames(trans_gene) <-c("chr","start","end","gene","cancer")

    trans_gene_Pos1 <- inner_join(trans_gene,gene_pos,by="gene")
    genome_chr <-c(1:22)
    trans_gene_Pos <-filter(trans_gene_Pos1,gene_chr%in%genome_chr)
    trans_gene_Pos$gene_chr <-paste0("chr",trans_gene_Pos$gene_chr)
    # trans_gene_Pos[is.na(trans_gene_Pos)] <- 0
    # trans_gene_Pos$
    trans_qtl <-trans_gene_Pos%>%dplyr::select(chr,start,end)
    # trans_qtl$value <-1
    trans_egene <-trans_gene_Pos%>%dplyr::select(gene_chr,gene_start,gene_end)
    # trans_egene$value <-1
    circos.clear()
    set.seed(80)
    figure_name <-paste0("09_",cancer,"_cis_trans_circos_assemble.pdf")
    pdf(figure_name)
    par(mar = c(1, 1, 1, 1))
    circos.par(start.degree = 90)
    circos.initializeWithIdeogram(species= "hg19",chromosome.index = paste0("chr", 1:22))
    circos.genomicDensity(cis, col = c("#3C8DAD"), track.height = 0.15)
    circos.genomicDensity(trans, col = c("#F5A962"), track.height = 0.15)
    # set.seed(16)
    circos.genomicLink(trans_qtl, trans_egene, col = sample(nrow(trans_qtl), nrow(trans_qtl), replace = TRUE))
    # circos.genomicLink(trans_qtl, trans_egene)
    dev.off()
    circos.clear()
}
    # dev.off()

#-----cancer normal
cancer = "BRCA"
normal <-read.table("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Breast_Mammary_Tissue/Cis_eQTL/hotspot_cis_eQTL/interval_18/Breast_Mammary_Tissue_segment_hotspot_cutoff_0.176.bed.gz",header = F,sep = "\t") %>% as.data.frame()
trans_file <- paste0("../../",cancer,"/Trans_eQTL/hotspot_trans_eQTL/interval_18/",cancer,"_segment_hotspot_cutoff_0.176.bed.gz")
cis_file <- paste0("../../",cancer,"/Cis_eQTL/hotspot_cis_eQTL/interval_18/",cancer,"_segment_hotspot_cutoff_0.176.bed.gz")
tran_gene_file <- paste0("../../",cancer,"/Trans_eQTL/hotspot_trans_eQTL/interval_18/hotspot_gene.bed.gz")

cis <- read.table(cis_file,header = F,sep = "\t") %>% as.data.frame()
trans <- read.table(trans_file,header = F,sep = "\t") %>% as.data.frame()
trans_gene <-read.table(tran_gene_file,header = F,sep = "\t") %>% as.data.frame()
gene_pos <-read.table("/home/huanhuan/ref_data/Ensembl/ensemble_v104_hg19_gene_pos.txt",header=T,sep="\t")%>%as.data.frame()
gene_pos <-gene_pos%>%dplyr::select(-c(1:2))
colnames(gene_pos) <-c("gene","gene_chr","gene_start","gene_end")
colnames(cis) <-c("chr","start","end")
colnames(trans) <-c("chr","start","end")
colnames(trans_gene) <-c("chr","start","end","gene","cancer")

trans_gene_Pos1 <- inner_join(trans_gene,gene_pos,by="gene")
genome_chr <-c(1:22)
trans_gene_Pos <-filter(trans_gene_Pos1,gene_chr%in%genome_chr)
trans_gene_Pos$gene_chr <-paste0("chr",trans_gene_Pos$gene_chr)
trans_qtl <-trans_gene_Pos%>%dplyr::select(chr,start,end)
trans_egene <-trans_gene_Pos%>%dplyr::select(gene_chr,gene_start,gene_end)
circos.clear()
set.seed(80)
figure_name <-paste0("09_",cancer,"_normal_cis_trans_circos_assemble.pdf")
pdf(figure_name)
par(mar = c(1, 1, 1, 1))
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species= "hg19",chromosome.index = paste0("chr", 1:22))
circos.genomicDensity(normal, col = c("#01937C"), track.height = 0.12,ylim.force=T,window.size = 1e3)
circos.genomicDensity(cis, col = c("#3C8DAD"), track.height = 0.12,ylim.force=T,window.size = 1e3)
circos.genomicDensity(trans, col = c("#F5A962"), track.height = 0.12,ylim.force=T,window.size = 1e3)
# set.seed(16)
circos.genomicLink(trans_qtl, trans_egene, col = sample(nrow(trans_qtl), nrow(trans_qtl), replace = TRUE))
# circos.genomicLink(trans_qtl, trans_egene)
dev.off()
# circos.clear()

#--------------------
circos.clear()
set.seed(80)
figure_name <-paste0("09_",cancer,"_cis_trans_circos_assemble_1kb.pdf")
pdf(figure_name)
par(mar = c(1, 1, 1, 1))
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species= "hg19",chromosome.index = paste0("chr", 1:22))
# circos.genomicDensity(normal, col = c("#01937C"), track.height = 0.12,ylim.force=T,window.size = 1e3)
circos.genomicDensity(cis, col = c("#3C8DAD"), track.height = 0.12,ylim.force=T,window.size = 1e3)
circos.genomicDensity(trans, col = c("#F5A962"), track.height = 0.12,ylim.force=T,window.size = 1e3)
# set.seed(16)
circos.genomicLink(trans_qtl, trans_egene, col = sample(nrow(trans_qtl), nrow(trans_qtl), replace = TRUE))
# circos.genomicLink(trans_qtl, trans_egene)
dev.off()
print("aaa")

circos.clear()
set.seed(80)
figure_name <-paste0("09_",cancer,"_normal_cis_trans_circos_assemble_10kb.pdf")
pdf(figure_name)
par(mar = c(1, 1, 1, 1))
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species= "hg19",chromosome.index = paste0("chr", 1:22))
circos.genomicDensity(normal, col = c("#01937C"), track.height = 0.12,ylim.force=T,window.size = 1e4)
circos.genomicDensity(cis, col = c("#3C8DAD"), track.height = 0.12,ylim.force=T,window.size = 1e4)
circos.genomicDensity(trans, col = c("#F5A962"), track.height = 0.12,ylim.force=T,window.size = 1e4)
# set.seed(16)
circos.genomicLink(trans_qtl, trans_egene, col = sample(nrow(trans_qtl), nrow(trans_qtl), replace = TRUE))
# circos.genomicLink(trans_qtl, trans_egene)
dev.off()
print("aaa")