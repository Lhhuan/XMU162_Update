
library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(VennDiagram)
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/homer_200/")
rm(list=ls())
q<-0.05

tf <-read.table(paste0("homer200_result_qvalue",q,".txt"),header = T,sep = "\t") %>% as.data.frame()
# tf$Cluster <-paste0("Cluster",tf$Cluster)
# tf$key =paste(tf$Motif_Name,tf$Consensus,sep="_")
ff <-function(i){
    i_c  <- filter(tf,Cluster==i)%>%select(Motif_name)%>%as.matrix()%>%as.vector()    #%>%as.character()
    return(i_c)
}
AAA <- list(C1=ff(1),C2=ff(2),C3=ff(3),C4=ff(4),C5=ff(5),C6=ff(6))
# my_c2 =c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")
library(venn)

pdf(paste0("12_02_Venn_Cluster_qvalue",q,"_motif.pdf"),height=3.1,width=3.1)
 venn(AAA, ilab = TRUE, zcolor = "style")
dev.off()

# pdf(paste0("12_02_Venn_Cluster7_qvalue",q,"_motif.pdf"),height=5,width=5)
#  P1 <- venn.diagram(x=AAA,filename =NULL,lwd=1,lty=2,col=my_c2 ,fill=my_c2,reverse=TRUE,cex=0.8,fontface = "bold",cat.cex=0.9,cat.fontfamily = "serif")
# #  print(P1)
# grid.draw(P1)
# dev.off()




# set.seed(12345)
# x <- list(First = 1:20, Second = 10:30, Third = sample(25:50, 15))
# venn(x, ilab = T, zcolor = "style")