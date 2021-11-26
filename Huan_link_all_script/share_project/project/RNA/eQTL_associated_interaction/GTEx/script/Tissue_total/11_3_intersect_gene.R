library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tidyverse)
library(pheatmap)
library(reshape2)
library(mygene)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_total/gene")
read_f<-function(tissue=NULL){
    file_name <-paste0("49_share_hotspot_",tissue,"_gene.txt.gz")
    org <- read.table(file_name,header = F,sep = "\t") %>% as.data.frame()
    colnames(org)[4] <-"Gene_ENSG"
    unique_gene <-unique(org$Gene_ENSG)
    return(unique_gene)
}
tissues <-c("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Cells_EBV-transformed_lymphocytes","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Kidney_Cortex","Muscle_Skeletal","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Uterus","Prostate","Brain_Cerebellar_Hemisphere","Testis","Brain_Nucleus_accumbens_basal_ganglia","Minor_Salivary_Gland","Cells_Cultured_fibroblasts","Pituitary","Vagina","Thyroid","Artery_Tibial","Artery_Coronary","Brain_Hypothalamus","Nerve_Tibial","Brain_Putamen_basal_ganglia","Brain_Amygdala","Breast_Mammary_Tissue","Liver","Lung","Ovary","Pancreas","Whole_Blood")


tmp <-lapply(tissues,read_f)
all_overlap_gene <-Reduce(intersect,tmp)%>%data.frame()
colnames(all_overlap_gene) <-"ENSG_ID"
en<-all_overlap_gene$ENSG_ID%>%unique()
df.0 <- queryMany(en, scopes="ensembl.gene", fields=c("symbol"), species="human")%>%as.data.frame()
final_gene <-df.0 %>%dplyr::select(query,symbol)
final_symbol <-unique(df.0$symbol)%>%unique()%>%as.data.frame()
# id<-df.0 %>% dplyr::select(query,symbol)
write.table(final_gene ,"11_3_49_overlap_gene.txt",row.names = F, col.names = T,quote =F,sep="\t")
write.table(final_symbol ,"11_3_49_overlap_gene_symbol.txt",row.names = F, col.names = F,quote =F,sep="\t")


all_union_gene <-Reduce(union,tmp)%>%data.frame()


df.0 <- queryMany(all_union_gene, scopes="ensembl.gene", fields=c("symbol"), species="human")%>%as.data.frame()
final_gene <-df.0 %>%dplyr::select(query,symbol)
final_symbol <-unique(df.0$symbol)%>%unique()%>%as.data.frame()

write.table(final_gene ,"11_3_49_union_gene.txt",row.names = F, col.names = T,quote =F,sep="\t")
write.table(final_symbol ,"11_3_49_union_gene_symbol.txt",row.names = F, col.names = F,quote =F,sep="\t")

# a <-Reduce(intersect,list(f,f1,f2,f3))