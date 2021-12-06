library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(VennDiagram)
library(venn)
setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/KKL/homer/")


tf <-read.table("knowResult_merge.txt",header = T,sep = "\t") %>% as.data.frame()
tf$key =paste(tf$Motif_Name,tf$Consensus,sep="_")
ff <-function(i){
    i_c  <- filter(tf,community==i)%>%select(key)%>%as.matrix()%>%as.vector()    #%>%as.character()
    return(i_c)
}

AAA <- list(Community1=ff(1),Community2=ff(2),Community3=ff(3),Community4=ff(4),Community5=ff(5),Community6=ff(6))
pdf("05_Venn_community_motif.pdf",height=5,width=5)
venn(AAA,ilab=TRUE, zcolor = "style")
dev.off()
# my_c2 =c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3","blue")
# pdf("05_Venn_community_motif.pdf",height=5,width=5)
#  P1 <- venn.diagram(x=AAA,color=my_c2)
# #  print(P1)
# grid.draw(P1)
# dev.off()
# ff <-function(i){
#     i_c  <- filter(tf,community==i)%>%select(Consensus)    #%>%as.character()
#     write.table(i_c,paste0(i,"consensus.txt"),col.names=F,row.names=F,quote=F,sep="\t")
# }

# lapply(1:6,ff)







# AAA <- list(C1=ff(1),C2=ff(2),C3=ff(3),C4=ff(4),C5=ff(5))
# my_c2=c(colors()[616], colors()[38], colors()[468],"#beca5c")
# my_c2=c(colors()[616], colors()[38], colors()[468],"#beca5c","#80ED99")
my_c2 =c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3","blue")
pdf("05_Venn_community_motif.pdf",height=5,width=5)
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