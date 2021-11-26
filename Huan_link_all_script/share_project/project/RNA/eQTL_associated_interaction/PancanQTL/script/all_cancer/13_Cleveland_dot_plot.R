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
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/cancer_total/")

org<-read.table("11_count_cancer_tissue_share_hotspot.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org$cancer_tissue <-paste(org$cancer,org$tissue,sep=":")

cancer_specific <- read.table("08_number_of_cancer_specific_hotspot.txt.gz",header = T,sep = "\t") %>% as.data.frame()
cancer_all <- read.table("08_number_of_cancer_all_hotspot.txt.gz",header = T,sep = "\t") %>% as.data.frame()

cancer <- inner_join(cancer_specific,cancer_all, by="cancer")
colnames(cancer)[1] <-"number_of_cancer_specfic_hotspot"
colnames(cancer)[4] <-"number_of_all_cancer_hotspot"
# cancer$cancer_diff <-cancer$number_of_all_cancer_hotspot - cancer$number_of_cancer_specfic_hotspot
cancer$cancer_tissue <-paste(cancer$cancer,cancer$tissue,sep =":")

all_cancer_info <-inner_join(cancer,org,by = "cancer_tissue")
all_cancer_info <- all_cancer_info %>% dplyr::select(-c(2,3,7,8))
all_cancer_info$relatively_share <- all_cancer_info$number_of_all_cancer_hotspot - all_cancer_info$number_of_cancer_specfic_hotspot -all_cancer_info$Number #Number: cancer tissue share 

# all_cancer <-all_cancer_info%>%dplyr::select(cancer_tissue,number_of_all_cancer_hotspot)
# colnames(all_cancer)[2] <-"Number"
# all_cancer$Class <- "All"
# all_cancer$cancer_tissue<-factor(all_cancer$cancer_tissue,levels= sort(all_cancer$cancer_tissue,decreasing = T))

cancer_specific <- all_cancer_info%>%dplyr::select(cancer_tissue,number_of_cancer_specfic_hotspot)
colnames(cancer_specific)[2] <-"Number"
cancer_specific$Class <- "Specific"
cancer_specific$cancer_tissue<-factor(cancer_specific$cancer_tissue,levels= sort(cancer_specific$cancer_tissue,decreasing = T))

cancer_share <- all_cancer_info%>%dplyr::select(cancer_tissue,Number)
cancer_share$Class <- "Share"
cancer_share$cancer_tissue<-factor(cancer_share$cancer_tissue,levels= sort(cancer_share$cancer_tissue,decreasing = T))

#---------
cancer_other <- all_cancer_info%>%dplyr::select(cancer_tissue,relatively_share)
colnames(cancer_other)[2] <-"Number"
cancer_other$Class <- "Other"
cancer_other$cancer_tissue<-factor(cancer_other$cancer_tissue,levels= sort(cancer_other$cancer_tissue,decreasing = T))

aaaaa <-bind_rows(cancer_specific,cancer_share,cancer_other)

color = "#5f939a"

p2 <-ggplot(aaaaa,mapping=aes(x=cancer_tissue,y=Number,fill=Class)) +
    p_theme+geom_segment(aes(xend =cancer_tissue),yend = 0,stat = "identity",position=position_dodge(width=0.5)) #position=position_dodge(0.7)

print(p2)
dev.off()


# p2 <-ggplot(aaaaa,mapping=aes(x=cancer_tissue,y=Number,fill=Class)) +
#     p_theme+geom_segment(aes(xend =cancer_tissue),yend = 0,stat = "identity" ,position_dodge(width=0.5) )  + # ,position="dodge" #(
#     geom_point(size = 2, color= "#008c9e")+
#     xlab("")+
#     theme(axis.text.y = element_text(size =9,color="black"),
#     axis.text.x = element_text(size =5.4,color="black",angle = 30,hjust=1))
#     # scale_y_continuous(expand = c(0,0))
# pdf("12_Cleveland_dot_plot_number_of_cancer_share_hotspot.pdf",height = 3,width = 5.5)
# print(p2)
# dev.off()


# org$cancer_tissue<-factor(org$cancer_tissue,levels= sort(org$cancer_tissue,decreasing = T))
co =c("#eba83a","#bb371a","#005792")
p2 <-ggplot(aaaaa,mapping=aes(x=cancer_tissue,y=Number,colour=Class)) +
    p_theme+geom_segment(aes(xend =cancer_tissue),yend = 0,colour = 'grey35')+
    geom_point(size = 0.8 )+#, color= "#008c9e"
    scale_colour_manual(values = co)+
    xlab("")+
    theme(axis.text.y = element_text(size =6,color="black"),
    axis.title.y = element_text(size =6,color="black"),
    axis.text.x = element_text(size =4.2,color="black",angle = 30,hjust=1),
    legend.title=element_blank(),
    legend.text=element_text(size=6),
    legend.key.size=unit(0.3, "lines"))
    # scale_y_continuous(expand = c(0,0))
pdf("13_Cleveland_dot_plot_number_of_cancer_hotspot.pdf",height = 3,width = 5.5)
print(p2)
dev.off()





# org$cancer_tissue<-factor(org$cancer_tissue,levels= sort(org$cancer_tissue,decreasing = T))
# p2 <-ggplot(org,mapping=aes(x=cancer_tissue,y=Number)) +
#     p_theme+geom_segment(aes(xend =cancer_tissue),yend = 0,colour = 'grey35')+
#     geom_point(size = 2, color= "#008c9e")+
#     xlab("")+
#     theme(axis.text.y = element_text(size =9,color="black"),
#     axis.text.x = element_text(size =5.4,color="black",angle = 30,hjust=1))
#     # scale_y_continuous(expand = c(0,0))
# pdf("12_Cleveland_dot_plot_number_of_cancer_share_hotspot___test.pdf",height = 3,width = 5.5)
# print(p2)
# dev.off()