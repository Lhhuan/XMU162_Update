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

# cutoffs= c(7, 24.25)

xQTL = "eQTL"
cutoffs= 7.3
#----------------------------------
current_directory <- paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/figure/ALL_",xQTL,"/", sep = "")
setwd(current_directory)

#---------------------------------
# interval<-seq(from=6, to=32, by=1)
# figure_list <-list ()
# Density_figure_list<-list()
# i=1
# for (j in interval){
#     input_file <-paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/output/ALL_eQTL/NHPoisson_emplambda_interval_",j,"cutoff_7.3_all_eQTL.txt.gz", sep = "")
#     org<-read.table(input_file,header = T,sep = "\t") %>% as.data.frame()
#     # org<-read.table("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/output/NHPoisson_emplambda_interval_6cutoff_7.3_all_eQTL.txt.gz",header = T,sep = "\t") %>% as.data.frame()

#     title_name<-paste("interval",j,sep = " ")
#     figure_list[[i]] <-ggplot(org, aes(x =emplambda)) +geom_histogram(position="identity")+
#         xlab("Emplambda") + ylab("Count")+p_theme+ggtitle(title_name)+xlim(0,1)+
#         theme(plot.title = element_text(hjust = 0.5))
        
#         # print(p1)
#     Density_figure_list[[i]] <- ggplot(org, aes(x =emplambda)) + geom_density()+
#         xlab("Emplambda") +ylab("Density")+p_theme+
#         ggtitle(title_name)+theme(plot.title = element_text(hjust = 0.5))
#     i=i+1
# }

interval2<-c(6,7,8,9,12,15,18)

figure_list2 <-list ()
Density_figure_list2<-list()

ii=1
for (j in interval2){
    input_file <-paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/output/ALL_eQTL/NHPoisson_emplambda_interval_",j,"cutoff_7.3_all_eQTL.txt.gz", sep = "")
    org<-read.table(input_file,header = T,sep = "\t") %>% as.data.frame()
    # org<-read.table("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/output/NHPoisson_emplambda_interval_6cutoff_7.3_all_eQTL.txt.gz",header = T,sep = "\t") %>% as.data.frame()

    title_name<-paste("interval",j,sep = " ")
    figure_list2[[ii]] <-ggplot(org, aes(x =emplambda)) +geom_histogram(position="identity")+
        xlab("Emplambda") + ylab("Count")+p_theme+ggtitle(title_name)+
        theme(plot.title = element_text(hjust = 0.5))
        
        # print(p1)
    Density_figure_list2[[ii]] <- ggplot(org, aes(x =emplambda)) + geom_density()+
        xlab("Emplambda") +ylab("Density")+p_theme+
        ggtitle(title_name)+theme(plot.title = element_text(hjust = 0.5))
    ii=ii+1
}

print("qqq")

# hist_figure_name1 <-paste("histogram_plot_cutoff_7.3_interval_6..32_all_",xQTL,".pdf", sep = "")
# pdf(hist_figure_name1,height = 7.5,width = 15) 
# p1<-marrangeGrob(figure_list,nrow=2,ncol=3)  
# print(p1)
# dev.off()

#-------------

# density_figure_name1 <-paste("density_plot_cutoff_7.3_interval_6..32_all_",xQTL,".pdf", sep = "")
# # pdf(density_figure_name1,height = 15,width = 30) 
# pdf(density_figure_name1,height = 7.5,width = 15) 
# p1<-marrangeGrob(Density_figure_list,nrow=2,ncol=3)  
# print(p1)
# dev.off()
#---------------

hist_figure_name2 <-paste("histogram_plot_cutoff_7.3_interval_ld_all_",xQTL,".pdf", sep = "")
pdf(hist_figure_name2,height = 7.5,width = 15) 
p1<-marrangeGrob(figure_list2,nrow=2,ncol=3)  
print(p1)
dev.off()

#-----------------------------
density_figure_name2 <-paste("density_plot_cutoff_7.3_interval_ld_all_",xQTL,".pdf", sep = "")
pdf(density_figure_name2,height = 7.5,width = 15) 
p1<-marrangeGrob(Density_figure_list2,nrow=2,ncol=3)  
print(p1)
dev.off()
print("qqq")
