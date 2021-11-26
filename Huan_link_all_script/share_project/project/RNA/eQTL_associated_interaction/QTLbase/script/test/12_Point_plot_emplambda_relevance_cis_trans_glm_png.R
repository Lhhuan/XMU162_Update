library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(Hmisc)
library(gcookbook)
library(ggpubr)


files<- list.files("~/project/RNA/eQTL_associated_interaction/QTLbase/output/xQTL_merge/Cis_Trans_1MB","txt.gz")
# n<-length(files)
fitting="glm"
i=1
for (file in files){
    name <- str_replace(file,"NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_1MB_","")  #replace
    name <- str_replace(name,".txt.gz","")  #replace
    # tmp <- str_split_fixed(names(kmpre$CellLabel),pattern='-',n=2)
    tmp <- str_split_fixed(name,pattern='_',n=4)
    xQTL1 = paste(tmp[1],tmp[2],sep =" ")
    xQTL2 = paste(tmp[3],tmp[4],sep =" ")
    x_QTL <-paste(xQTL1,"emplambda", sep = " ")
    y_QTL <-paste(xQTL2,"emplambda", sep = " ") 
    x_QTL <- capitalize(x_QTL)
    y_QTL <- capitalize(y_QTL)
    setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/xQTL_merge/Cis_Trans_1MB/")
    org<-fread(file,header = T,sep = "\t") %>% as.data.frame()
    n_r<-nrow(org)
    if (n_r >1){
        colnames(org)[1] <- "x_QTL"
        colnames(org)[2] <- "y_QTL" #change colname name 
        # round(org$x_QTL,4)
        setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/xQTL_merge/Cis_Trans_1MB_fitting/Cis_Trans_glm_1MB/")
        figure_name <-name <- str_replace(file,".txt.gz",".png")
        png(figure_name)

        p1<-ggplot(org, aes(x =x_QTL, y=y_QTL)) +geom_point(size=0.1, color ="black" ) 
        p1<-p1+stat_smooth(method=fitting)
        p1<-p1+stat_cor(data=org, method = "pearson")
        p2<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
            axis.line = element_line(colour = "black") )
        p2<-p2+xlab(x_QTL) + ylab(y_QTL)
        p2<-p2+ggtitle("1MB glm") 
        print(p2)
        dev.off()
        print(i)
        i=i+1
    }
}

files<- list.files("~/project//RNA/eQTL_associated_interaction/QTLbase/output/xQTL_merge/Cis_Trans_10MB","txt.gz")
# n<-length(files)
i=1
for (file in files){
    name <- str_replace(file,"NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_10MB_","")  #replace
    name <- str_replace(name,".txt.gz","")  #replace
    # tmp <- str_split_fixed(names(kmpre$CellLabel),pattern='-',n=2)
    tmp <- str_split_fixed(name,pattern='_',n=4)
    xQTL1 = paste(tmp[1],tmp[2],sep =" ")
    xQTL2 = paste(tmp[3],tmp[4],sep =" ")
    x_QTL <-paste(xQTL1,"emplambda", sep = " ")
    y_QTL <-paste(xQTL2,"emplambda", sep = " ") 
    x_QTL <- capitalize(x_QTL)
    y_QTL <- capitalize(y_QTL)
    setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/xQTL_merge/Cis_Trans_10MB/")
    org<-fread(file,header = T,sep = "\t") %>% as.data.frame()
    n_r<-nrow(org)
    if (n_r >1){
        colnames(org)[1] <- "x_QTL"
        colnames(org)[2] <- "y_QTL" #change colname name 
        # round(org$x_QTL,4)
        setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/xQTL_merge/Cis_Trans_10MB_fitting/Cis_Trans_glm_10MB/")
        figure_name <-name <- str_replace(file,".txt.gz",".png")
        png(figure_name)

        p1<-ggplot(org, aes(x =x_QTL, y=y_QTL)) +geom_point(size=0.1, color ="black" ) 
        p1<-p1+stat_smooth(method=fitting)
        p1<-p1+stat_cor(data=org, method = "pearson")
        p2<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
            axis.line = element_line(colour = "black") )
        p2<-p2+xlab(x_QTL) + ylab(y_QTL)
        p2<-p2+ggtitle("10MB glm") 
        print(p2)
        dev.off()
        print(i)
        i=i+1
    }
}

