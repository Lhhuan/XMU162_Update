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
    for(i in c(1:22)){
    file_name <-paste("../output/NHPoisson_emplambda_interval",j,".txt", sep = "")
    org<-read.table(file_name,header = T,sep = "\t") %>% as.data.frame()

    # aa<-head(org1,n=10000)
    figure_name <-paste("Bar_plot_of_emplambda_int",j,"_chr",i,".png", sep = "")
    # figure_name <-paste("Bar_plot_of_emplambda_int",j,"_chr",i,".pdf", sep = "")
    title_name<-paste("Bar plot of emplambda int",j, "bp","chr",i,sep = " ")
    # pdf(figure_name,height = 7,width = 7) 
    png(figure_name,height = 800,width = 1000) 
    # p1<-ggplot(aa, aes(x = t, y =emplambda ))+ 
    p1<-ggplot(org1, aes(x = t, y =emplambda ))+
        geom_bar(stat = "identity")+
        xlab("Pos")+
        p.theme+
        ylab("Emplambda")+
        ggtitle(title_name)+
        theme(plot.title = element_text(hjust = 0.5))
    print(p1)
    dev.off()
}


#------------------------------------------------------