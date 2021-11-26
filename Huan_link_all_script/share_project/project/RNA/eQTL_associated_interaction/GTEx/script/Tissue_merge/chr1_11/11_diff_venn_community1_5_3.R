library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(VennDiagram)
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/chr1_11/6kmer/")


c31 <-read.table("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_11/kmer/6/3_community/homer/communities_1.bed",header = F,sep = "\t") %>% as.data.frame()
c51 <-read.table("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_11/kmer/6/5_community/homer/communities_1.bed",header = F,sep = "\t") %>% as.data.frame()

A =paste0(c31$V1,c31$V2,c31$V3)
B =paste0(c51$V1,c51$V2,c51$V3)


# tf$key =paste(tf$Motif_Name,tf$Consensus,sep="_")
# ff <-function(i){
#     i_c  <- filter(tf,community==i)%>%select(Consensus)%>%as.matrix()%>%as.vector()    #%>%as.character()
#     return(i_c)
# }
# AAA <- list(Community1=ff(1),Community2=ff(2),Community3=ff(3))
# AAA <- list(C1=ff(1),C2=ff(2),C3=ff(3),C4=ff(4),C5=ff(5))
# my_c2=c(colors()[616], colors()[38], colors()[468],"#beca5c")
# my_c2=c(colors()[616], colors()[38], colors()[468],"#beca5c","#80ED99")
# my_c2 =c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")
AAA <-list(C3_1=A,C5_1=B)
my_c2 =c("dodgerblue", "goldenrod1")

pdf("10_Venn_community1_5_3.pdf",height=5,width=5)
 P1 <- venn.diagram(x=AAA,filename =NULL,lwd=1,lty=2,col=my_c2 ,fill=my_c2,reverse=TRUE,cex=0.8,fontface = "bold",cat.cex=0.9,cat.fontfamily = "serif")
#  print(P1)
grid.draw(P1)
dev.off()


#-----------------------------------

c31 <-read.table("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_11/kmer/6/3_community/homer/communities_0.bed",header = F,sep = "\t") %>% as.data.frame()
c51 <-read.table("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_11/kmer/6/5_community/homer/communities_0.bed",header = F,sep = "\t") %>% as.data.frame()

A =paste0(c31$V1,c31$V2,c31$V3)
B =paste0(c51$V1,c51$V2,c51$V3)
AAA <-list(C3_0=A,C5_0=B)

pdf("10_Venn_community0_5_3.pdf",height=5,width=5)
 P1 <- venn.diagram(x=AAA,filename =NULL,lwd=1,lty=2,col=my_c2 ,fill=my_c2,reverse=TRUE,cex=0.8,fontface = "bold",cat.cex=0.9,cat.fontfamily = "serif")
#  print(P1)
grid.draw(P1)
dev.off()



c31 <-read.table("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_11/kmer/6/6_community/homer/communities_0.bed",header = F,sep = "\t") %>% as.data.frame()
c51 <-read.table("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_11/kmer/6/5_community/homer/communities_0.bed",header = F,sep = "\t") %>% as.data.frame()

A =paste0(c31$V1,c31$V2,c31$V3)
B =paste0(c51$V1,c51$V2,c51$V3)
AAA <-list(C6_0=A,C5_0=B)

pdf("10_Venn_community0_5_6.pdf",height=5,width=5)
 P1 <- venn.diagram(x=AAA,filename =NULL,lwd=1,lty=2,col=my_c2 ,fill=my_c2,reverse=TRUE,cex=0.8,fontface = "bold",cat.cex=0.9,cat.fontfamily = "serif")
#  print(P1)
grid.draw(P1)
dev.off()