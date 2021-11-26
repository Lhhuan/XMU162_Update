library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(CMplot)
library(data.table)


QTLs = c("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL","QTL")
n<-length(QTLs)
# for (xQTL in QTLs)
xQTL = "mQTL"

p.theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
# fun=function(i){
#     xQTL =QTLs[i]
current_directory <- paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/ALL_",xQTL,"/", sep = "")
setwd(current_directory)

# cutoff =24.25
cutoff =7.3
intervals =1000
j=intervals 
input_file <-paste("../../output/ALL_", xQTL, "/NHPoisson_emplambda_interval_",j,"cutoff_",cutoff,"_all_",xQTL,".txt", sep = "")
org<-fread(input_file,header = T,sep = "\t") %>% as.data.frame()
org1<-filter(org, !is.na(emplambda))
pdf("cutoff7_density_plot_in_different_interval_all_eQTL.pdf",height = 8,width = 10) 
    p <-ggplot (org1,aes(x=emplambda))+ 
        geom_density(alpha=.2) + 
        xlab("Emplambda")+
        p.theme+
        ylab("Density")+
        ggtitle(title_name)+
        theme(plot.title = element_text(hjust = 0.5))
    

# }


# library(parallel)
# a=mclapply(1:n,function(x){fun(x)},mc.cores=9)
#-------------------------------------------------------------


p.theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
interval1<-seq(from=2000, to=8000, by=2000)
interval2<-seq(from=10000, to=60000, by=10000)
interval=c(interval1,interval2) #两个interval合并
#-----------------cutoff=7
figure_list <-list ()
i=1
for (j in interval){
    file_name <-paste("../output/ALL_eQTL/NHPoisson_emplambda_interval_",j,"_all_eQTLbase.txt", sep = "")
    org1<-read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
    title_name<-paste("Density of emplambda int",j, "bp",sep = " ")
    figure_list[[i]] <-ggplot (org1,aes(x=emplambda))+ 
        geom_density(alpha=.2) + 
        xlab("Emplambda")+
        p.theme+
        ylab("Density")+
        ggtitle(title_name)+
        theme(plot.title = element_text(hjust = 0.5))
        i =i+1
}
pdf("cutoff7_density_plot_in_different_interval_all_eQTL.pdf",height = 40,width = 56) 
marrangeGrob(figure_list,nrow=4,ncol=3,top = "cutoff = 7" )  #指定list的具体位置，也可以直接用整个list组装，比如：figure_list
dev.off()

# #-------cutoff=9
# interval<-seq(from=10000, to=60000, by=10000)
# figure_list <-list ()
# i=1
# for (j in interval){
#     # file_name <-paste("../output/NHPoisson_emplambda_interval",j,".txt", sep = "")
#     file_name <-paste("../output/NHPoisson_emplambda_interval_",j,"_cutoff_9.txt", sep = "")
#     org1<-read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
#     title_name<-paste("Density of emplambda int",j, "bp",sep = " ")
#     figure_list[[i]]  <-ggplot (org1,aes(x=emplambda))+ 
#         geom_density(alpha=.2) + 
#         xlab("Emplambda")+
#         p.theme+
#         ylab("Density")+
#         ggtitle(title_name)+
#         theme(plot.title = element_text(hjust = 0.5)) 
#         i=i+1   
#     # print(p1)
# }

# pdf("cutoff9_density_plot_in_different_interval.pdf",height = 20,width = 28) 
# marrangeGrob(figure_list,nrow=3,ncol=2,top = "cutoff = 9" )  #指定list的具体位置，也可以直接用整个list组装，比如：figure_list
# # marrangeGrob(figure_list[1:13],nrow=4,ncol=4)  #指定list的具体位置，也可以直接用整个list组装，比如：figure_list
# dev.off()