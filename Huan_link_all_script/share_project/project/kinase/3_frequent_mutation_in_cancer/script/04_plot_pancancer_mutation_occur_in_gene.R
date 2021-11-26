library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
setwd("/home/huanhuan/project/kinase/3_frequent_mutation_in_cancer/figure/")
dit <-"/home/huanhuan/project/kinase/3_frequent_mutation_in_cancer/output/"

org<-read.table(file.path(dit, "02_pancancer_count_mutation_occur_in_gene.txt"),header = T,sep = "\t") %>% as.data.frame()

pdf("04_pancancer_mutation_occur_in_gene.pdf",height = 5,width = 4)
# p1<-ggplot(org, aes(x=SYMBOL, y=number, group=1)) + geom_line(linetype="solid", color ="#4a9ff5") +geom_point(size=0.6, color ="#4a9ff5" )
p1<-ggplot(org, aes(x = reorder(SYMBOL,number), y=number, group=1)) +geom_point(size=0.6, color ="black" ) 
p1<-p1 +coord_flip() #x和y交换坐标轴
p2<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.title.y = element_text(size = 13),axis.title.x = element_text(size = 13),axis.line = element_line(colour = "black"),
    axis.text.y = element_text(size=5)) 
p2<-p2+xlab("SYMBOL") + ylab("Mutation number")
p2
dev.off()