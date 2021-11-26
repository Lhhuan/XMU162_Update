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

intervals=c(6,7,8,9,12,15,18)

types= c("cis_1MB","cis_10MB","trans_1MB","trans_10MB")

library(Hmisc)

for(j in intervals){
    i=1
    figure_list <-list ()
    Density_figure_list<-list()
    setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/output/ALL_eQTL/cis_trans/NHP/")
    for (type in types){
        org<-read.table(paste0(type,"/NHPoisson_emplambda_interval_",j,"_cutoff_7.3_",type,"_eQTL.txt.gz"),header = T,sep = "\t") %>% as.data.frame()

            # org<-read.table(input_file,header = T,sep = "\t") %>% as.data.frame()
        # org<-read.table("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/NHPoisson_emplambda_interval_6cutoff_7.3_all_eQTL.txt.gz",header = T,sep = "\t") %>% as.data.frame()
        new_type <-str_replace(type,"_"," ")
        new_type <-capitalize(new_type)
        # title_name<-paste(new_type,": interval ",j,sep = "")
        title_name<- new_type
        figure_list[[i]] <-ggplot(org, aes(x =emplambda)) +geom_histogram(position="identity")+
            xlab("Emplambda") + ylab("Count")+p_theme+ggtitle(title_name)+coord_cartesian(xlim = c(0, 1))+
            theme(plot.title = element_text(hjust = 0.5))
            
            # print(p1)
        Density_figure_list[[i]] <- ggplot(org, aes(x =emplambda)) + geom_density()+
            xlab("Emplambda") +ylab("Density")+p_theme+
            ggtitle(title_name)+theme(plot.title = element_text(hjust = 0.5))
        i=i+1
    }


    setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/figure/ALL_eQTL/cis_trans/")
    top_title <-paste("Interval",j,sep=" ")
    pdf(paste0("histogram_plot_cutoff_7.3_interval_",j,"_cis_trans.pdf"),height = 7.5,width = 7) 
    p1<-marrangeGrob(figure_list,nrow=2,ncol=2,top = top_title)  
    print(p1)
    dev.off()

    pdf(paste0("density_plot_cutoff_7.3_interval_",j,"_cis_trans.pdf"),height = 7.5,width = 7) 
    p1<-marrangeGrob(Density_figure_list,nrow=2,ncol=2,top = top_title)  
    print(p1)
    dev.off()
    print(j)
}

