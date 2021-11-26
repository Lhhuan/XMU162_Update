library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(Seurat)

# cancers <-c("BRCA","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","STAD","THCA")
cancers <-c("BRCA","KIRC","KIRP","LAML","LIHC","LUAD","OV","PAAD","PRAD","STAD","THCA")


p_theme<-theme(panel.grid =element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))



credplot.gg_main <- function(cancer=NULL){
    path = paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/",cancer,"/Cis_eQTL/anno/")
    setwd(path)
    org <- read.table("29_share_specific_fisher_test_two_side.txt",header = T,sep = "\t") %>% as.data.frame()
    org <-org%>%filter(p_value<=0.05)
    # org <-org%>%arrange(OR)
    if(nrow(org)>0){
        org <-org%>%arrange(OR)

    
        org$marker <- factor(org$marker, levels=org$marker) #marker ranking
        p1<-ggplot(data=org, aes(x=marker,y=log(OR), ymin=log(conf_int_bottom) , ymax=log(conf_int_up)))+
            geom_pointrange(size = 0.3)+
            coord_flip() +  # flip coordinates (puts labels on y axis)
            ylab("Odds ratio (log scale)")+xlab("Markers")+ 
            p_theme+theme(axis.text.x = element_text(color = "black"),axis.text.y = element_text(color = "black")) +
            geom_text(aes(x=marker,y=0.95,label = significant)) 
        
        figure_name <-paste0("30_",cancer,"_specific_share_fisher_test_forest_plot.pdf")
        setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/cancer_total/")
        pdf(figure_name,height = 3.5,width = 4) #
        print(p1)
        dev.off()
    }
}

lapply(cancers,credplot.gg_main)



cancer = "BRCA"

    path = paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/",cancer,"/Cis_eQTL/anno/")
    setwd(path)
    org <- read.table("29_share_specific_fisher_test_two_side.txt",header = T,sep = "\t") %>% as.data.frame()
    # org <-org%>%filter(p_value<=0.05)
    # org <-org%>%arrange(OR)
    if(nrow(org)>0){
        org <-org%>%arrange(OR)

        org$marker <-str_replace(org$marker,"CHROMATIN_Accessibility","CA")
        # org$marker <- factor(org$marker, levels=org$marker) #marker ranking
        org$marker <- factor(org$marker, levels=rev(c("CA","TFBS","H3K27ac","H3K4me1","H3K4me3","H3K9ac","H3K27me3")))
        p1<-ggplot(data=org, aes(x=marker,y=log(OR), ymin=log(conf_int_bottom) , ymax=log(conf_int_up)))+
            geom_pointrange(size = 0.3)+
            coord_flip() +  # flip coordinates (puts labels on y axis)
            ylab("Odds ratio (log scale)")+xlab("Markers")+ 
            p_theme+theme(axis.text.x = element_text(color = "black"),axis.text.y = element_text(color = "black"), plot.title = element_text(hjust = 0.5)) +
            geom_text(aes(x=marker,y=0.95,label = significant)) +
            ggtitle(cancer)
        
        figure_name <-paste0("30_",cancer,"_specific_share_fisher_test_forest_plot_A.pdf")
        setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/cancer_total/")
        pdf(figure_name,height = 3.5,width = 3.8) #
        print(p1)
        dev.off()
    }
#----------
cancer ="PRAD"
    path = paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/",cancer,"/Cis_eQTL/anno/")
    setwd(path)
    org <- read.table("29_share_specific_fisher_test_two_side.txt",header = T,sep = "\t") %>% as.data.frame()
    # org <-org%>%filter(p_value<=0.05)
    # org <-org%>%arrange(OR)
    if(nrow(org)>0){
        org <-org%>%arrange(OR)

        org$marker <-str_replace(org$marker,"CHROMATIN_Accessibility","CA")
        # org$marker <- factor(org$marker, levels=org$marker) #marker ranking
        org$marker <- factor(org$marker, levels=rev(c("CA","TFBS","H3K27ac","H3K4me1","H3K4me3","H3K36me3","H3K27me3","H3K9me3")))
        p1<-ggplot(data=org, aes(x=marker,y=log(OR), ymin=log(conf_int_bottom) , ymax=log(conf_int_up)))+
            geom_pointrange(size = 0.3)+
            coord_flip() +  # flip coordinates (puts labels on y axis)
            ylab("Odds ratio (log scale)")+xlab("Markers")+ 
            p_theme+theme(axis.text.x = element_text(color = "black"),axis.text.y = element_text(color = "black"), plot.title = element_text(hjust = 0.5)) +
            geom_text(aes(x=marker,y=2.95,label = significant)) +
            ggtitle(cancer)
        
        figure_name <-paste0("30_",cancer,"_specific_share_fisher_test_forest_plot_A.pdf")
        setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/cancer_total/")
        pdf(figure_name,height = 3.5,width = 3.8) #
        print(p1)
        dev.off()
    }