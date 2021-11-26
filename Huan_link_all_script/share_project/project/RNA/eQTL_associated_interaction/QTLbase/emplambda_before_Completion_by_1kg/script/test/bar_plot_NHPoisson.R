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
for (j in interval){
    file_name <-paste("../output/NHPoisson_emplambda_interval",j,".txt", sep = "")
    org1<-read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
    # org<-filter(org1, emplambda != "NA")
    pos<-paste(org1$chr,org1$t,sep="_")
    date<-cbind(emplambda=org1$emplambda, pos)%>%as.data.frame()

    figure_name <-paste("Bar_plot_of_emplambda_int",j,".png", sep = "")
    title_name<-paste("Bar plot of emplambda int",j, "bp",sep = " ")
    png(figure_name,height = 800,width = 1000) 
    p1<-ggplot(date, aes(x = pos, y =emplambda ))+ 
        geom_bar(stat = "identity")+
        xlab("Emplambda")+
        p.theme+
        ylab("pos")+
        ggtitle(title_name)+
        theme(plot.title = element_text(hjust = 0.5))
    print(p1)
    dev.off()
}


#------------------------------------------------------