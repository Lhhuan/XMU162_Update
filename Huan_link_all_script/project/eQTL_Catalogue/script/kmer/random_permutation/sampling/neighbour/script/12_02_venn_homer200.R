
library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(VennDiagram)
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/homer_200")
rm(list=ls())
q<-0.1

tf <-read.table(paste0("homer200_result_qvalue",q,".txt"),header = T,sep = "\t") %>% as.data.frame()
# tf$Cluster <-paste0("Cluster",tf$Cluster)
# tf$key =paste(tf$Motif_Name,tf$Consensus,sep="_")
ff <-function(i){
    i_c  <- filter(tf,Cluster==i)%>%select(Motif_name)%>%as.matrix()%>%as.vector()    #%>%as.character()
    return(i_c)
}
AAA <- list(C1=ff(1),C2=ff(2),C3=ff(3),C6=ff(6),C7=ff(7))
# AAA <- list(C1=ff(1),C2=ff(2),C3=ff(3),C4=ff(4),C5=ff(5))
# my_c2=c(colors()[616], colors()[38], colors()[468],"#beca5c")
# my_c2=c(colors()[616], colors()[38], colors()[468],"#beca5c","#80ED99")
my_c2 =c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")
# my_c2 =c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3")

pdf(paste0("12_02_Venn_Cluster7_qvalue",q,"_motif.pdf"),height=5,width=5)
 P1 <- venn.diagram(x=AAA,filename =NULL,lwd=1,lty=2,col=my_c2 ,fill=my_c2,reverse=TRUE,cex=0.8,fontface = "bold",cat.cex=0.9,cat.fontfamily = "serif")
#  print(P1)
grid.draw(P1)
dev.off()





