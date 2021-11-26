library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/")


p.theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
interval<-seq(from=10000, to=60000, by=10000)
#-----------------cutoff=7
for (j in interval){
    file_name <-paste("../output/NHPoisson_emplambda_interval",j,".txt", sep = "")
    org1<-read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
    # org<-filter(org1, emplambda != "NA")


    figure_name <-paste("Density_of_emplambda_int",j,".pdf", sep = "")
    title_name<-paste("Density of emplambda int",j, "bp",sep = " ")
    pdf(figure_name,height = 4,width = 5) 
    p1 <-ggplot (org1,aes(x=emplambda))+ #不把chr转成factor就不会分染色体，因为chr那列是纯数字  #,colour = as.factor(chr)
        geom_density(alpha=.2) + 
        xlab("Emplambda")+
        p.theme+
        ylab("Density")+
        ggtitle(title_name)+
        theme(plot.title = element_text(hjust = 0.5))
        # labs(colour = "Chr") ##修改图例名字
    print(p1)
    dev.off()
}


#-------cutoff=9

for (j in interval){
    # file_name <-paste("../output/NHPoisson_emplambda_interval",j,".txt", sep = "")
    file_name <-paste("../output/NHPoisson_emplambda_interval_cutoff_9_",j,".txt", sep = "")
    org1<-read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
    # org<-filter(org1, emplambda != "NA")


    figure_name <-paste("Density_of_emplambda_int",j,"_cutoff_9.pdf", sep = "")
    title_name<-paste("Density of emplambda int",j, "bp cutoff 9",sep = " ")
    pdf(figure_name,height = 4,width = 5) 
    p1 <-ggplot (org1,aes(x=emplambda))+ #不把chr转成factor就不会分染色体，因为chr那列是纯数字  #,colour = as.factor(chr)
        geom_density(alpha=.2) + 
        xlab("Emplambda")+
        p.theme+
        ylab("Density")+
        ggtitle(title_name)+
        theme(plot.title = element_text(hjust = 0.5))
        # labs(colour = "Chr") ##修改图例名字
    print(p1)
    dev.off()
}