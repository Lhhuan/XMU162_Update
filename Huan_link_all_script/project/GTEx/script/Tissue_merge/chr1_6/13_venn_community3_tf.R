library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(VennDiagram)
setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/chr1_6/6kmer/3_community/")


tf <-read.table("knowResult_merge.txt",header = T,sep = "\t") %>% as.data.frame()
tf$key =paste(tf$Motif_Name,tf$Consensus,sep="_")
ff <-function(i){
    i_c  <- filter(tf,community==i)%>%select(key)%>%as.matrix()%>%as.vector()    #%>%as.character()
    return(i_c)
}
AAA <- list(Community1=ff(1),Community2=ff(2),Community3=ff(3))
my_c2 =c("dodgerblue", "goldenrod1", "darkorange1")

pdf("13_Venn_community3_motif.pdf",height=5,width=5)
 P1 <- venn.diagram(x=AAA,filename =NULL,lwd=1,lty=2,col=my_c2 ,fill=my_c2,reverse=TRUE,cex=0.8,fontface = "bold",cat.cex=0.9,cat.fontfamily = "serif")
#  print(P1)
grid.draw(P1)
dev.off()
