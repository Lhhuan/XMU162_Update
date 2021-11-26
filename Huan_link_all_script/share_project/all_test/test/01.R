library(tidyverse)
library(stringr)
library(mygene)
setwd("~/test")
tfs1 <- read_tsv("07_somatic_snv_indel_mutationID_gene_geneLevel.txt",col_names = T)

en <- tfs1$Gene%>%unique()
df.0 <- queryMany(en, scopes="ensembl.gene", fields=c("entrezgene","symbol"), species="human")
write.table(df.0,"08_ensg_to_entrezid.txt",row.names = F, col.names = T,quote = F, sep = "\t")
