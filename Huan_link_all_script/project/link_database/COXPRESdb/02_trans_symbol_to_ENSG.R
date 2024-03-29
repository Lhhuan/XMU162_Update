library(tidyverse)
library(stringr)
library(GenomicFeatures)
library(readr)
library(mygene)
setwd("~/project/link_database/COXPRESdb/")
org<-read.table("01_G16808_S85825_get_gene_file_list.txt.gz",header = F,sep = "\t") %>% as.data.frame()

g <-c(org$V1)%>%unique()
df.0 <- queryMany(g, scopes="entrezgene", fields=c("symbol","ensembl.gene"), species="human")%>% as.data.frame()
df.1 <- as_tibble(df.0) %>% filter(is.na(notfound)) %>% group_by(query) %>% top_n(1,X_score)%>%data.frame()
write.table(df.1,"02_G16808_S85825_entrezgene_symbol_ensembl.gene.txt",quote = F,sep = "\t",row.names = F)
