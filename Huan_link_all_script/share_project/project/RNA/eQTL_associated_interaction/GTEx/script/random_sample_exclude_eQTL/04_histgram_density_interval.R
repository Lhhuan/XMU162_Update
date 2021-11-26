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

# xQTL = "eQTL"
cutoff= 7.3

tissues = c("Lung","Whole_Blood")
# intervals=c(6,7,8,9,12,15,18)
j=18
nums <-c(1:10)

library(Hmisc)



for(tissue in tissues){
    setwd(paste0("/share/data0/QTLbase/huan/GTEx/random_select_exclude_eQTL/",tissue,"/figure/"))
    figure_list <-list ()
    i=1
    print(i)
    for(number in nums){
        org<-read.table(paste0("../NHP/NHPoisson_emplambda_interval_",j,"_cutoff_",cutoff,"_",tissue,"_",number,".txt.gz"),header = T,sep = "\t") %>% as.data.frame()
        title_name <-paste0("Random ",number)
        figure_list[[i]] <-ggplot(org, aes(x =emplambda)) +geom_histogram(position="identity")+
        xlab("Emplambda") + ylab("Count")+p_theme+ggtitle(title_name)+coord_cartesian(xlim = c(0, 1))+
        theme(plot.title = element_text(hjust = 0.5))
        i = i+1
    }
    top_title <-paste0("Random select Cis eQTL in ",tissue)
    # pdf(paste0("histogram_plot_cutoff_7.3_cis_eQTL_in_",tissue,".pdf"),compress = TRUE) 
    pdf(paste0("histogram_plot_cutoff_7.3_cis_eQTL_in_",tissue,".pdf"))
    p1<-marrangeGrob(figure_list,nrow=2,ncol=3,top = top_title)  
    print(p1)
    dev.off()
}



