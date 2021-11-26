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

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/cancer_total/share/total/gene")
read_f<-function(tissue=NULL){
    file_name <-paste0("19_share_hotspot_",tissue,"_gene.txt.gz")
    org <- read.table(file_name,header = F,sep = "\t") %>% as.data.frame()
    colnames(org)[4] <-"Gene_symbol"
    unique_gene <-unique(org$Gene_symbol)
    return(unique_gene)
}
tissues <-c("ACC","BRCA","COAD","ESCA","KICH","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","THCA","UCEC","UCS")


tmp <-lapply(tissues,read_f)
all_union_gene <-Reduce(union,tmp)%>%data.frame()

write.table(all_union_gene ,"35_3_19_share_hotspot_union_gene.txt",row.names = F, col.names = F,quote =F,sep="\t")





# all_overlap_gene <-Reduce(intersect,tmp)%>%data.frame()
# colnames(all_overlap_gene) <-"ENSG_ID"
# en<-all_overlap_gene$ENSG_ID%>%unique()
# df.0 <- queryMany(en, scopes="ensembl.gene", fields=c("symbol"), species="human")%>%as.data.frame()
# final_gene <-df.0 %>%dplyr::select(query,symbol)
# final_symbol <-unique(df.0$symbol)%>%unique()%>%as.data.frame()
# # id<-df.0 %>% dplyr::select(query,symbol)
# write.table(final_gene ,"11_3_49_overlap_gene.txt",row.names = F, col.names = T,quote =F,sep="\t")
# write.table(final_symbol ,"11_3_49_overlap_gene_symbol.txt",row.names = F, col.names = F,quote =F,sep="\t")


# all_union_gene <-Reduce(union,tmp)%>%data.frame()


# df.0 <- queryMany(all_union_gene, scopes="ensembl.gene", fields=c("symbol"), species="human")%>%as.data.frame()
# final_gene <-df.0 %>%dplyr::select(query,symbol)
# final_symbol <-unique(df.0$symbol)%>%unique()%>%as.data.frame()

# write.table(final_gene ,"11_3_49_union_gene.txt",row.names = F, col.names = T,quote =F,sep="\t")
# write.table(final_symbol ,"11_3_49_union_gene_symbol.txt",row.names = F, col.names = F,quote =F,sep="\t")

# # a <-Reduce(intersect,list(f,f1,f2,f3))