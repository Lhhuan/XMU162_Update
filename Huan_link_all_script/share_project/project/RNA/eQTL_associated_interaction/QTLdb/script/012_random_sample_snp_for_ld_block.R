library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)



setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/script/")
# org<-read.table("../output/merge_QTL_all_QTLtype_pop.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org<-read.table("../output/011_eur_1kg_phase3_v5_hg19_snp.txt.gz",header = T,sep = "\t") %>% as.data.frame()
colnames(org)[1] <-"SNP_chr"
colnames(org)[2] <-"SNP_pos"
unique_snp <-org%>%dplyr::select(SNP_chr,SNP_pos)%>%unique()

rs<-data.frame()
for(i in c(1:10)){
    j=i+40
set.seed(j)
new_df <- unique_snp %>% group_by(SNP_chr)%>% summarise(SNP_pos = paste(base::sample(SNP_pos, 2), collapse=","))%>% as.data.frame()
rs<-bind_rows(rs,new_df)
print(i)
}

rs_chr <- filter(rs,SNP_chr!=23,SNP_chr!=24)
tmp <- str_split_fixed(rs_chr$SNP_pos,pattern=',',n=2)
rs_chr$pos1 <- tmp[,1]
rs_chr$pos2 <- tmp[,2]
rs_chr$pos1 <-paste(rs_chr$SNP_chr,rs_chr$pos1,sep="_")
rs_chr$pos2 <-paste(rs_chr$SNP_chr,rs_chr$pos2,sep="_")
A <-rs_chr$pos1%>%as.data.frame()
B <-rs_chr$pos2%>%as.data.frame()
final_snp <-bind_rows(A,B)

gz1 <- gzfile("../output/012_random_sample_snp_for_ld_block.txt.gz", "w")
#write.csv(df1, gz1)
write.table(final_snp,gz1,row.names = F, col.names = F,quote =F,sep="\t")
close(gz1)
print("aaa")
