library(tidyverse)
library(stringr)
library(GenomicFeatures)
library(readr)
library(mygene)
setwd("~/project/link_database/reactome_FI/")



# tfs <- read_tsv(file_path('data/kinase.txt', col_names=F)
org<-read.table("FIsInGene_122921_with_annotations.txt",header = T,sep = "\t") %>% as.data.frame()

g <-c(org$Gene1,org$Gene2)%>%unique()
df.0 <- queryMany(g, scopes="symbol", fields=c("symbol","entrezgene","ensembl.gene"), species="human")%>% as.data.frame()
df.tmp <- as_tibble(df.0) %>% filter(!is.na(notfound))#查找不成功的
df.00 <- queryMany(df.tmp$query, scopes="alias", fields=c("symbol","entrezgene","ensembl.gene"), species="human")%>% as.data.frame() #查找不成功的用别名重新查找
df.err <- as_tibble(df.00) %>% filter(!is.na(notfound))
df.true <- as_tibble(df.00) %>% filter(is.na(notfound))
re <-bind_rows(df.0,df.true)
# re <- re %>%d
df.1 <- as_tibble(re) %>% filter(is.na(notfound)) %>% group_by(query) %>% top_n(1,X_score)%>%data.frame()

# a <-lapply(1:nrow(df.1),function(x) df.1$ensembl[[x]][[1,1]])
# df.1$ensg <-unlist(a)
write.table(df.1,"./output/01_symbol_ENSG.txt",quote = F,sep = "\t",row.names = F)

# id<-df.0 %>% dplyr::select(query,entrezgene,symbol,ensembl)
# df.00 <- queryMany(df.tmp$query, scopes="alias", fields=c("symbol"), species="human") #查找不成功的用别名重新查找。
# df.2 <- as_tibble(df.00) %>%filter(is.na(notfound)) %>% group_by(query) %>% top_n(1,X_score) #注意这里面有分数相同的，在perl里读的时候只读第一个
# df.err <- as_tibble(df.00) %>% filter(!is.na(notfound
