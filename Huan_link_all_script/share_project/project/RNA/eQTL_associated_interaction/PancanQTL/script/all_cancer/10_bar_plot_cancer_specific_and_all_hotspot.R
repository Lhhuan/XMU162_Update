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

cancer_all <- read.table("08_number_of_cancer_all_hotspot.txt.gz",header = T,sep = "\t") %>% as.data.frame()
cancer_sp <- read.table("08_number_of_cancer_specific_hotspot.txt.gz",header = T,sep = "\t") %>% as.data.frame()
cancer_sp$cancer_tissue <- paste(cancer_sp$cancer,cancer_sp$tissue,sep =":")

cancer <-inner_join(cancer_all,cancer_sp,by="cancer")
cancer$diff <-cancer$Number.x - cancer$Number.y
cancer_diff  <- cancer%>%dplyr::select(cancer_tissue,diff)
colnames(cancer_diff)[2] <- "Number"
cancer_diff$Class <- "All hotspot"

cancer_sp <- cancer_sp %>% dplyr::select(cancer_tissue,Number)
cancer_sp$Class <- "Cancer specific hotspot"
cancer_diff$cancer_tissue<-factor(cancer_diff$cancer_tissue,levels= sort(cancer_diff$cancer_tissue,decreasing = T))
cancer_sp$cancer_tissue<-factor(cancer_sp$cancer_tissue,levels= sort(cancer_sp$cancer_tissue,decreasing = T))

cancer_f <- bind_rows(cancer_diff,cancer_sp)
pdf("10_bar_plot_cancer_specific_and_all_tissue.pdf")

# pdf("13_tissue_specific_and_all_tissue.pdf",height = 7,width = 7.5)
p1 <-ggplot(cancer_f,aes(x=cancer_tissue,y=Number,fill=Class))+
    geom_bar(stat='identity',width=0.5 )+
    scale_fill_manual(values=c("#f29191","#1eae98"))+
    p_theme+xlab("Cancer") +ylab("Number of hotspots") +
    coord_flip()+
    scale_y_continuous(expand=c(0,0))+
    theme(axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black"),
    legend.title = element_blank(),
    legend.text = element_text(size=7))

print(p1)   
dev.off()

cancer <-inner_join(cancer_all,cancer_sp,by="cancer")
colnames(cancer)[2] <-"number_of_all_hotspot"
colnames(cancer)[4] <-"number_of_cancer_specific_hotspot"
cancer$cancer_specific_ratio <-cancer$number_of_cancer_specific_hotspot/cancer$number_of_all_hotspot *100

write.table(cancer ,"10_cancer_specific_ratio.txt",row.names = F, col.names = T,quote =F,sep="\t")






# p2 <-ggplot(cancer_f,aes(x=cancer_tissue,y=negative,fill=Class))+
#     geom_bar(stat='identity',width=0.5 )+
#     scale_fill_manual(values=c("#f29191","#1eae98"))+# fill  color
#     p_theme+
#     labs(x="",y="",title="Number of Cancer hotspot")+
#     coord_flip()+
#     theme(
#         axis.text.x = element_text(color="black",angle=90,vjust=0,hjust=1),
#         axis.ticks.y = element_blank(), 
#         axis.text.y=element_text(color="black",size = 7.5),
#         axis.line.y=element_blank(),
#         legend.position = "none", 
#         plot.title=element_text(hjust = 0.5))+
#         scale_y_continuous (limits = c(-66000,0),breaks = c(-60000,-40000,-20000,0),labels = c("60,000","40,000","20,000",0),expand=c(0,0))
# print(p2)
# dev.off()

# p3 <-p2+p11
# pdf("10_barplot_of_cancer_normal_hotspot.pdf",height = 6,width = 8.5)
# print(p3)
# dev.off()

