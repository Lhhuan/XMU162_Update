library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/script/")


p.theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))

population<-c("AFR","AMR","EAS","EUR","SAS","MIX")
for (pop in population){
    interval<-seq(from=10000, to=60000, by=10000)
    #-----------------cutoff=7
    figure_list <-list ()
    i=1
    for (j in interval){
        file_name <-paste("../output/population_specific/NHPoisson_emplambda_interval_",j,"_all_QTLbase_",pop,"_specific.txt", sep = "")
        org1<-read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
        title_name<-paste("Density of emplambda int",j, "bp_",pop,sep = " ")
        figure_list[[i]] <-ggplot (org1,aes(x=emplambda))+ 
            geom_density(alpha=.2) + 
            xlab("Emplambda")+
            p.theme+
            ylab("Density")+
            ggtitle(title_name)+
            theme(plot.title = element_text(hjust = 0.5))
            i =i+1
    }
    figure_name <-paste("../figure/population_specific/cutoff7_density_plot_in_different_interval_all_QTL_",pop,"_.pdf", sep = "")
    pdf(figure_name,height = 20,width = 28) 
    marrangeGrob(figure_list,nrow=3,ncol=2,top = "cutoff = 7" )  #指定list的具体位置，也可以直接用整个list组装，比如：figure_list
    dev.off()
}

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