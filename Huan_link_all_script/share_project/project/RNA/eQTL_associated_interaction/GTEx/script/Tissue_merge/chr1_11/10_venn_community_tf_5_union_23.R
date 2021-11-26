library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(VennDiagram)
library(tidyverse)
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/chr1_11/6kmer/5_community")


tf1 <-read.table("knowResult_merge.txt",header = T,sep = "\t") %>% as.data.frame()

tf2 <-read_tsv("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_11/kmer/6/5_community/homer/2_3/knownResults.txt",col_names=TRUE) %>% as.data.frame()
colnames(tf2)[c(1,2,5)] <-c("Motif_Name","Consensus","q_value")
tf2 <-filter(tf2,q_value <0.05)%>%select(Motif_Name,Consensus,q_value)
tf2$community <-"2_3"
tf1 <-filter(tf1,community ==1 |community>3)
tf1$community <-as.character(tf1$community)
tf <-bind_rows(tf1,tf2)
tf$key =paste(tf$Motif_Name,tf$Consensus,sep="_")
ff <-function(i){
    i_c  <- filter(tf,community==i)%>%select(Consensus)%>%as.matrix()%>%as.vector()    #%>%as.character()
    return(i_c)
}
AAA <- list(Community1=ff("1"),Community2_3=ff("2_3"),Community4=ff("4"),Community5=ff("5"))
# AAA <- list(C1=ff(1),C2=ff(2),C3=ff(3),C4=ff(4),C5=ff(5))
# my_c2=c(colors()[616], colors()[38], colors()[468],"#beca5c")
# my_c2=c(colors()[616], colors()[38], colors()[468],"#beca5c","#80ED99")
# my_c2 =c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")
my_c2 =c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3")

pdf("10_Venn_community_motif_2_3.pdf",height=5,width=5)
 P1 <- venn.diagram(x=AAA,filename =NULL,lwd=1,lty=2,col=my_c2 ,fill=my_c2,reverse=TRUE,cex=0.8,fontface = "bold",cat.cex=0.9,cat.fontfamily = "serif")
#  print(P1)
grid.draw(P1)
dev.off()

# pdf("10_Venn_community_motif1.pdf",height=5,width=5)
#  P1 <- venn.diagram(x=AAA,filename =NULL,lwd=1,lty=2,col=my_c2 ,fill=my_c2,reverse=TRUE,cex=0.8,fontface = "bold",cat.cex=0.9)
# #  print(P1)
# grid.draw(P1)
# dev.off()

# P1 <- venn.diagram(x=AAA,"bbbb.png",height = 1500, width = 1600, resolution =300, imagetype="png", lwd=1,lty=2,col=my_c2 ,fill=my_c2,reverse=TRUE,fontface = "bold",cat.fontfamily = "serif")

# pdf("aaaa.pdf",height=5,width=5)
#  P1 <- venn.diagram(x=AAA,filename ="NULL")
# #  print(P1)
# grid.draw(P1)
# dev.off()






# pdf("aaaa.pdf",height=5,width=5)
#  P1 <- venn.diagram(x=AAA,filename =NULL,lwd=1,lty=2,col=my_c2 ,fill=my_c2,reverse=TRUE,cex=0.8,fontface = "bold",cat.cex=0.8,cat.fontfamily = "serif")
# #  print(P1)
# grid.draw(P1)
# dev.off()

# P1 <- venn.diagram(x=AAA,"bbbb.png",height = 1500, width = 1500, resolution =300, imagetype="png", lwd=1,lty=2,col=my_c2 ,fill=my_c2,reverse=TRUE,fontface = "bold",cat.fontfamily = "serif")



# my_c2=c(colors()[616], colors()[38], colors()[468],"#beca5c","#80ED99")