library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(devtools)
library(cowplot)
library(gridExtra)
library(ggpubr)

setwd("/home/huanhuan/project/kinase/3_frequent_mutation_in_cancer/figure/")
dit <-"/home/huanhuan/project/kinase/3_frequent_mutation_in_cancer/output/"

org<-read.table(file.path(dit, "03_cancer_specific_count_mutation_occur_in_gene.txt"),header = T,sep = "\t")

type <-org %>% dplyr::select(cancer_type)%>%unique()
n1<-nrow(type)
figure_list <-list ()
# par(mfrow=c(6,6))
for(i in c(1:n1)){
    Cancer = type[i,1]
    cancer_specific<-filter(org,cancer_type == Cancer)
    # p<-Cancer
    figure_name<-paste("04_",Cancer,"_mutation_occur_in_gene.pdf",sep = "")
    title_name<-paste(Cancer)
    # p<-paste("p",i,sep = "")
    # pdf(figure_name,height = 5,width = 9)
    pdf(figure_name,height = 5,width = 3.5)
    # p<-ggplot(org, aes(x=SYMBOL, y=number, group=1)) + geom_line(linetype="solid", color ="#4a9ff5") +geom_point(size=0.6, color ="#4a9ff5" )
    p<-ggplot(cancer_specific, aes(x = reorder(SYMBOL,number), y=number,fill = Cancer)) +geom_point(size=0.6,color ="#4a9ff5") 
    p<-p +coord_flip() #x和y交换坐标轴
    p<-p+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.title.y = element_text(size = 9),axis.title.x = element_text(size = 9),axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=6),axis.text.x = element_text(size=9)) 
    p<-p+xlab("SYMBOL") + ylab("Mutation number")
    p<-p+guides(fill=FALSE) #删除图例
    # p<-p + labs("Cancer") #修改图例
    # p<-p+theme(legend.title = element_text(size=8),legend.text = element_text(size = 8)) #,legend.position = 'top'
    # p2<-p2+guides(fill=)
    p<-p+ggtitle(title_name)+theme(plot.title=element_text(size=9))
    print(p)
    dev.off()
    figure_list[[i]] <- p #把图片存进list
    # append(figure_list,figure_name)
}
#---------------------------------------------------把多张图组在一起
pdf("04_1_6_data.pdf",height = 12,width = 10) 
marrangeGrob(figure_list[1:6],nrow=2,ncol=3)  #指定list的具体位置，也可以直接用整个list组装，比如：figure_list
dev.off()
pdf("04_7_12_data.pdf",height = 12,width = 10)
marrangeGrob(figure_list[7:12],nrow=2,ncol=3)
dev.off()
pdf("04_13_18_data.pdf",height = 12,width = 10)
marrangeGrob(figure_list[13:18],nrow=2,ncol=3)
dev.off()
pdf("04_19_24_data.pdf",height = 12,width = 10)
marrangeGrob(figure_list[19:24],nrow=2,ncol=3)
dev.off()
pdf("04_25_30_data.pdf",height = 12,width = 10)
marrangeGrob(figure_list[25:30],nrow=2,ncol=3)
dev.off()
pdf("04_31_33_data.pdf",height = 6,width = 10)
marrangeGrob(figure_list[31:33],nrow=1,ncol=3)
dev.off()
# figure_list[[1]]
# pdf("04_all_data.pdf",height = 10,width = 10)
# pp<-multiplot(p1, p2, p3, p4, cols=2)
# print(pp)
# dev.off()
# plotlist p[1:4]