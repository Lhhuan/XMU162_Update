library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
#-------------------------------------------------------
tissue = "Lung"

org<-read.table(paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/",tissue,"/",tissue,"_cis_eQTL_in_hotspot_ratio.gz"),header = T,sep = "\t") %>% as.data.frame()

# org1 <-filter(org,Cutoff >0.25)
org$Cutoff<-as.factor(org$Cutoff)
setwd(paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/figure/",tissue))
pdf("23_eQTL_ratio_in_hotspot.pdf",height = 7,width =10)
p1 <- ggplot(org, aes(x =Cutoff,y=eQTL_ratio_in_hotspot)) +geom_boxplot(outlier.colour = NA,fill = "#5698c3")
p1 <-p1 +p_theme + ylab("eQTL ratio in hotspot")
print(p1)
dev.off()

