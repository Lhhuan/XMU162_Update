library(tidyverse)
library(stringr)
library(GenomicFeatures)
library(readr)
library(mygene)
setwd("~/project/kinase/")
di<-"~/project/kinase"


# tfs <- read_tsv(file_path('data/kinase.txt', col_names=F)
tfs<-read.table(file.path(di, "kinase_ppi/data/kinase.txt"),header = T,sep = "\t") %>% as.data.frame()
g<-tfs$gene %>% unique()
df.0 <- queryMany(g, scopes="symbol", fields=c("symbol","entrezgene","ensembl.gene"), species="human")%>% as.data.frame()
# df.1 <- as_tibble(df.0) %>% filter(is.na(notfound)) %>% group_by(query) %>% top_n(1,X_score)
# df.tmp <- as_tibble(df.0) %>% filter(!is.na(notfound))#查找不成功的
id<-df.0 %>% dplyr::select(query,entrezgene,symbol,ensembl)

# df.00 <- queryMany(df.tmp$query, scopes="alias", fields=c("symbol"), species="human") #查找不成功的用别名重新查找。
# df.2 <- as_tibble(df.00) %>%filter(is.na(notfound)) %>% group_by(query) %>% top_n(1,X_score) #注意这里面有分数相同的，在perl里读的时候只读第一个
# df.err <- as_tibble(df.00) %>% filter(!is.na(notfound))


# id.1 <- df.1 %>% dplyr::select(query,entrezgene)
# id.2 <- df.2 %>% dplyr::select(query,entrezgene)
# id.e <- df.err %>% dplyr::select(query,entrezgene)
# id <- bind_rows(id.1,id.2) 

write.table(id,file.path(di,"2_strong_gene-impact/output/03_transfrom_kinsae_gene_to_enstrez_ENSG.txt"),quote = F,sep = "\t",row.names = F)

