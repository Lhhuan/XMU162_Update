library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/")

org<-read.table("../output/count_of_different_distance_of_region.txt",header = T,sep = "\t") %>% as.data.frame()


pdf("bar_plot_distance_of_locus_region.pdf",height = 4,width = 5) #把图片存下来
p<-ggplot(org,aes(x = distance, y = number)) + geom_bar(stat = 'identity',fill = "black")
p1<-p+scale_x_discrete(limits= c("0-20000","20000-40000","40000-60000","60000-80000","80000-100000","0-100000","100000-200000","200000-300000","300000-400000","400000-500000")) #调整坐标轴上的顺序
p1<-p1+xlab("Distance of locus region")+ylab("Locus number") #修改坐标轴标签的文本
p1<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),axis.text.x=element_text(angle = 90,size =7),axis.text.y=element_text(size =7),
                                                axis.line = element_line(colour = "black")) #去背景
# p1<-p1+scale_fill_discrete(limits= c("SNV/Indel","CNV", "Deletion","Duplication","Inversion","Translocation")) #修改图例顺序
p1
dev.off()

org1<-read.table("../output/count_of_different_distance_of_region_QTL_specific.txt",header = T,sep = "\t") %>% as.data.frame()
#-------------------------all QTL type
pdf("bar_plot_distance_of_locus_region_pop_specific.pdf",height = 3.5,width = 5) #把图片存下来
p<-ggplot(org1,aes(x = distance, y = number,fill = QTL_type)) + geom_bar(stat = 'identity', position = 'dodge')
p1<-p+scale_x_discrete(limits= c("0-20000","20000-40000","40000-60000","60000-80000","80000-100000","0-100000","100000-200000","200000-300000","300000-400000","400000-500000")) #调整坐标轴上的顺序
p1<-p1+xlab("Distance of locus region")+ylab("Locus number") #修改坐标轴标签的文本
p1<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),axis.text.x=element_text(angle = 90,size =6),axis.text.y=element_text(size =6),
                                                axis.line = element_line(colour = "black")) #去背景
# p1<-p1+scale_fill_discrete(limits= c("SNV/Indel","CNV", "Deletion","Duplication","Inversion","Translocation")) #修改图例顺序
p1

dev.off()


#--------------------------------------------
#-----------------------------------------------------------------------------------qtl specific
type <-org1 %>% dplyr::select(QTL_type)%>%unique()
n1<-nrow(type)
figure_list <-list ()
#----------------------------------------------------------
for(i in c(1:n1)){
    qtl = type[i,1]
    QTL_specific<-filter(org1,QTL_type == qtl)
    # p<-Cancer
    figure_name<-paste(qtl,"_bar_plot_distance_of_locus_region.pdf",sep = "")
    title_name<-paste(qtl)
    pdf(figure_name,height = 4,width = 5) #把图片存下来
    p<-ggplot(QTL_specific,aes(x = distance, y = number)) + geom_bar(stat = 'identity', fill = "black")
    p1<-p+scale_x_discrete(limits= c("0-20000","20000-40000","40000-60000","60000-80000","80000-100000","0-100000","100000-200000","200000-300000","300000-400000","400000-500000")) #调整坐标轴上的顺序
    p1<-p1+xlab("Distance of locus region")+ylab("Locus number") #修改坐标轴标签的文本
    p1<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                    panel.background = element_blank(), axis.title.y = element_text(size = 9),
                                                    axis.title.x = element_text(size = 9),axis.text.x=element_text(angle = 90,size =6),axis.text.y=element_text(size =6),
                                                    axis.line = element_line(colour = "black")) #去背景
    # p1<-p1+scale_fill_discrete(limits= c("SNV/Indel","CNV", "Deletion","Duplication","Inversion","Translocation")) #修改图例顺序
    p1<-p1+ggtitle(title_name)+theme(plot.title=element_text(size=9,hjust = 0.5))
    print(p1)
    dev.off()
    figure_list[[i]] <- p #把图片存进list
}

#--------------------------------------------

p.theme <- theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                    panel.background = element_blank(), axis.title.y = element_text(size = 8),
                                                    axis.title.x = element_text(size = 8),axis.text.x=element_text(angle = 90,size =5),axis.text.y=element_text(size =8),
                                             axis.line = element_line(colour = "black"))


#-------------------------------存到一起

for(i in c(1:n1)){
    qtl = type[i,1]
    cancer_specific<-filter(org1,QTL_type == qtl)
    # p<-Cancer
    figure_name<-paste(qtl,"_bar_plot_distance_of_locus_region.pdf",sep = "")
    title_name<-paste(qtl)
    #pdf(figure_name,height = 4,width = 5) #把图片存下来
    figure_list[[i]] <- print(ggplot(cancer_specific,aes(x = distance, y = number)) + 
        geom_bar(stat = 'identity', fill = "black")+
        scale_x_discrete(limits= c("0-20000","20000-40000","40000-60000","60000-80000","80000-100000","0-100000","100000-200000","200000-300000","300000-400000","400000-500000")) + 
        xlab("Distance of locus region") +
        ylab("Locus number") +# theme_bw() + 
        p.theme+
        ggtitle(title_name)+theme(plot.title=element_text(size=9,hjust = 0.5))  )
    #print(p)
   # dev.off()
    #figure_list[[i]] <- p #把图片存进list
}

pdf("QTL_specific_bar_plot_distance_of_locus_region.pdf",height = 16,width = 16) 
marrangeGrob(figure_list,nrow=4,ncol=4)  #指定list的具体位置，也可以直接用整个list组装，比如：figure_list
# marrangeGrob(figure_list[1:13],nrow=4,ncol=4)  #指定list的具体位置，也可以直接用整个list组装，比如：figure_list
dev.off()
#--------------------------------------------------------

a <- lapply(c(1:n1),function(i){
    qtl = type[i,1]
    cancer_specific<-filter(org1,QTL_type == qtl)
    # p<-Cancer
    figure_name<-paste(qtl,"_bar_plot_distance_of_locus_region.pdf",sep = "")
    title_name<-paste(qtl)
    #pdf(figure_name,height = 4,width = 5) #把图片存下来
    print(ggplot(cancer_specific,aes(x = distance, y = number)) + 
        geom_bar(stat = 'identity', fill = "black")+
        scale_x_discrete(limits= c("0-20000","20000-40000","40000-60000","60000-80000","80000-100000","0-100000","100000-200000","200000-300000","300000-400000","400000-500000")) + 
        xlab("Distance of locus region") +
        ylab("Locus number") +# theme_bw() + 
        p.theme+
        ggtitle(title_name)+theme(plot.title=element_text(size=9,hjust = 0.5))  )

})

pdf("QTL_specific_bar_plot_distance_of_locus_region.pdf",height = 16,width = 16) 
marrangeGrob(a,nrow=4,ncol=4)
dev.off()


#------------------------------------------------------------

#---------------------------------------------------把多张图组在一起
pdf("QTL_specific_bar_plot_distance_of_locus_region.pdf",height = 16,width = 16) 
marrangeGrob(figure_list,nrow=4,ncol=4)  #指定list的具体位置，也可以直接用整个list组装，比如：figure_list
# marrangeGrob(figure_list[1:13],nrow=4,ncol=4)  #指定list的具体位置，也可以直接用整个list组装，比如：figure_list
dev.off()
pdf("QTL_specific_7_13_bar_plot_distance_of_locus_region.pdf",height = 12,width = 10)
marrangeGrob(figure_list[7:12],nrow=3,ncol=3)








#--------------------------------------------------------





