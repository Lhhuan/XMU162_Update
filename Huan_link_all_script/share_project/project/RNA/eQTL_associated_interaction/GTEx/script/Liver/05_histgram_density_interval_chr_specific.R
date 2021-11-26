library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(parallel)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
#-------------------------------------------------------

# cutoffs= c(7, 24.25)

xQTL = "eQTL"
cutoff= 7.3
tissue = "Lung"
setwd(paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/figure/",tissue,"/chr/"))
# setwd(paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/",tissue,"/Cis_eQTL/NHP/"))

# for (j in interval){
ProcessBedGz <- function(j = NULL){
    figure_list<-list()
    Density_figure_list <-list()
    i=1
    input_file <-paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/",tissue,"/Cis_eQTL/NHP/NHPoisson_emplambda_interval_",j,"_cutoff_",cutoff,"_",tissue,".txt.gz", sep = "")
    org<-read.table(input_file,header = T,sep = "\t") %>% as.data.frame()
    # org<-read.table("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/NHPoisson_emplambda_interval_6cutoff_7.3_all_eQTL.txt.gz",header = T,sep = "\t") %>% as.data.frame()
    # for(number in c(1:22)){
    for(number in c(1:22)){
        org_chr <-filter(org,chr==number)
        title_name<-paste("chr",i,sep = " ")
        figure_list[[i]] <-ggplot(org_chr, aes(x =emplambda)) +geom_histogram(position="identity")+
            xlab("Emplambda") + ylab("Count")+p_theme+ggtitle(title_name)+coord_cartesian(xlim = c(0, 1))+
            theme(plot.title = element_text(hjust = 0.5))
            
            # print(p1)
        Density_figure_list[[i]] <- ggplot(org_chr, aes(x =emplambda)) + geom_density()+
            xlab("Emplambda") +ylab("Density")+p_theme+coord_cartesian(xlim = c(0, 1))+
            ggtitle(title_name)+theme(plot.title = element_text(hjust = 0.5))
        i=i+1
    }
    hist_figure_name1 <-paste0("histogram_plot_cutoff_7.3_interval_",j,"_cis_eQTL_in_",tissue,"_chr_specific.pdf")
    pdf(hist_figure_name1) 
    p1<-marrangeGrob(figure_list,nrow=3,ncol=3)  
    print(p1)
    dev.off()
    hist_figure_name2 <-paste0("density_plot_cutoff_7.3_interval_",j,"_cis_eQTL_in_",tissue,"_chr_specific.pdf")
    pdf(hist_figure_name2) 
    p1<-marrangeGrob(Density_figure_list,nrow=3,ncol=4)  
    print(p1)
    dev.off()
    print("aaa")
}

# interval <-c(6,7,8,9,12,15,19)
# mclapply(interval, ProcessBedGz, mc.cores = 4)
interval <-c(18)
mclapply(interval, ProcessBedGz, mc.cores = 4)