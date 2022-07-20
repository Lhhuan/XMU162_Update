library(tidyverse)
library(stringr)
library(GenomicFeatures)
library(readr)
library(mygene)
setwd("~/project/link_database/reactome_FI/output")



# tfs <- read_tsv(file_path('data/kinase.txt', col_names=F)
org<-read.table("../FIsInGene_122921_with_annotations.txt",header = T,sep = "\t") %>% as.data.frame()
eg <-read.table("021_query_gene_ensembl_pos.txt",header = T,sep = "\t") %>% as.data.frame()
eg <-filter(eg,ensembl !="NULL")
colnames(eg)[1] <-"Gene1"
org1 <-inner_join(org,eg,by="Gene1")
colnames(org1)[6:9] <-c("Gene1_ensembl","Gene1_chr","Gene1_start","Gene1_end")
colnames(eg)[1] <-"Gene2"
org2 <-inner_join(org1,eg,by="Gene2")
colnames(org2)[10:13] <-c("Gene2_ensembl","Gene2_chr","Gene2_start","Gene2_end")
org2 <-filter(org2,Gene2_chr ==Gene1_chr)
# org2 <-filter(org2,Gene1_ensembl!="NULL" & Gene2_ensembl!="NULL" )
unpred <- filter(org2,Annotation != "predicted")

write.table(org2,"03_FIsInGene_122921_with_annotations_ENSG.txt",quote = F,sep = "\t",row.names = F)
write.table(unpred,"03_FIsInGene_122921_with_annotations_ENSG_filter_predicted.txt",quote = F,sep = "\t",row.names = F)