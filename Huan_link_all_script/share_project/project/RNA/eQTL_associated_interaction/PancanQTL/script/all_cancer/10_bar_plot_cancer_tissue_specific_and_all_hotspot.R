library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(patchwork)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
#-------------------------------------------------------
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/cancer_total/")
#-------------normal
normal <- read.table("09_count_cancer_related_tissue_specific_hotspot.txt.gz",header = T,sep = "\t") %>% as.data.frame()

normal$cancer_tissue <-paste(normal$cancer,normal$tissue,sep =":")
normal$diff <- normal$number_of_cancer_relate_tissue_all_hotspot - normal$number_of_cancer_relate_tissue_specific_hotspot

normal_sp <-normal%>%dplyr::select(cancer_tissue,number_of_cancer_relate_tissue_specific_hotspot)
colnames(normal_sp)[2] <-"Number"
normal_sp$Class <- "Specific"
normal_sp$cancer_tissue<-factor(normal_sp$cancer_tissue,levels= sort(normal_sp$cancer_tissue,decreasing = T))

normal_all <-normal%>%dplyr::select(cancer_tissue,diff)
colnames(normal_all)[2] <-"Number"
normal_all$Class <- "All"
normal_all$cancer_tissue<-factor(normal_all$cancer_tissue,levels= sort(normal_all$cancer_tissue,decreasing = T))

f_normal <-bind_rows(normal_sp,normal_all)

p1 <-ggplot(f_normal,aes(x=cancer_tissue,y=Number,fill=Class))+
    geom_bar(stat='identity',width=0.5 )+
    scale_fill_manual(values=c("#f29191","#1eae98"))+# fill  color
    p_theme+
    labs(x="",y="",title="Number of Normal hotspot")+
    coord_flip()
print(p1)  
dev.off()
p11 <- p1 + theme(
        axis.text.x = element_text(color="black",angle=90,vjust=1,hjust=1),
        axis.ticks.y = element_blank(), 
        axis.text.y=element_blank(),
        axis.line.y=element_blank(),
        legend.title = element_blank(), 
        plot.title=element_text(hjust = 0.5))+
        scale_y_continuous (limits = c(0,66000),breaks = c(0,20000,40000,60000),labels = c(0,"20,000","40,000","60,000"),expand=c(0,0))

print(p11)  
dev.off()
#----------------------------------------------cancer

cancer_all <- read.table("08_number_of_cancer_all_hotspot.txt.gz",header = T,sep = "\t") %>% as.data.frame()
cancer_sp <- read.table("08_number_of_cancer_specific_hotspot.txt.gz",header = T,sep = "\t") %>% as.data.frame()
cancer_sp$cancer_tissue <- paste(cancer_sp$cancer,cancer_sp$tissue,sep =":")

cancer <-inner_join(cancer_all,cancer_sp,by="cancer")
cancer$diff <-cancer$Number.x - cancer$Number.y
cancer_diff  <- cancer%>%dplyr::select(cancer_tissue,diff)
colnames(cancer_diff)[2] <- "Number"
cancer_diff$Class <- "All"

cancer_sp <- cancer_sp %>% dplyr::select(cancer_tissue,Number)
cancer_sp$Class <- "Specific"
cancer_diff$cancer_tissue<-factor(cancer_diff$cancer_tissue,levels= sort(cancer_diff$cancer_tissue,decreasing = T))
cancer_sp$cancer_tissue<-factor(cancer_sp$cancer_tissue,levels= sort(cancer_sp$cancer_tissue,decreasing = T))

cancer_f <- bind_rows(cancer_diff,cancer_sp)
cancer_f$negative <--cancer_f$Number
p2 <-ggplot(cancer_f,aes(x=cancer_tissue,y=negative,fill=Class))+
    geom_bar(stat='identity',width=0.5 )+
    scale_fill_manual(values=c("#f29191","#1eae98"))+# fill  color
    p_theme+
    labs(x="",y="",title="Number of Cancer hotspot")+
    coord_flip()+
    theme(
        axis.text.x = element_text(color="black",angle=90,vjust=0,hjust=1),
        axis.ticks.y = element_blank(), 
        axis.text.y=element_text(color="black",size = 7.5),
        axis.line.y=element_blank(),
        legend.position = "none", 
        plot.title=element_text(hjust = 0.5))+
        scale_y_continuous (limits = c(-66000,0),breaks = c(-60000,-40000,-20000,0),labels = c("60,000","40,000","20,000",0),expand=c(0,0))
print(p2)
dev.off()

p3 <-p2+p11
pdf("10_barplot_of_cancer_normal_hotspot.pdf",height = 6,width = 8.5)
print(p3)
dev.off()

