library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(ggpubr)


files<- list.files("~/project/RNA/eQTL_associated_interaction/QTLbase/output/xQTL_merge/Cis_1MB","txt.gz")
# n<-length(files)
i=1
j=1
for (file in files){
    name <- str_replace(file,"NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_cis_1MB_","")  #replace
    name <- str_replace(name,".txt.gz","")  #replace
    # tmp <- str_split_fixed(names(kmpre$CellLabel),pattern='-',n=2)
    tmp <- str_split_fixed(name,pattern='_',n=2)
    xQTL1 = tmp[1]
    xQTL2 = tmp[2]
    x_QTL <-paste(xQTL1," lambda", sep = "")
    y_QTL <-paste(xQTL2," lambda", sep = "")
    setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/xQTL_merge/Cis_1MB/")
    org<-fread(file,header = T,sep = "\t") %>% as.data.frame()
    n_r<-nrow(org)
    if (n_r >1){
        colnames(org)[1] <- "x_QTL"
        colnames(org)[2] <- "y_QTL" #change colname name 
        # round(org$x_QTL,4)
        figure_name <-name <- str_replace(file,".txt.gz",".png")
        fittings= c("lm","glm","loess")
        for(fitting in fittings){
            dir = paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/xQTL_merge/Cis_1MB_fitting/Cis_1MB_",fitting,sep="")
            title_name = paste("Cis 1MB", fitting,sep=" ")
            setwd(dir)
            png(figure_name)

            p1<-ggplot(org, aes(x =x_QTL, y=y_QTL)) +geom_point(size=0.1, color ="black" ) 
            p1<-p1+stat_smooth(method=fitting)
            p1<-p1+stat_cor(data=org, method = "pearson")
            p2<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                axis.line = element_line(colour = "black") )
            p2<-p2+xlab(x_QTL) + ylab(y_QTL)
            p2<-p2+ggtitle(title_name) 
            print(p2)
            dev.off()
            i=i+1
            print(i)
            print(fitting)
        }
        j=j+1
        print(j)
    }

}

files<- list.files("~/project//RNA/eQTL_associated_interaction/QTLbase/output/xQTL_merge/Cis_10MB","txt.gz")
# n<-length(files)
i=1
j=1
for (file in files){
    name <- str_replace(file,"NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_cis_10MB_","")  #replace
    name <- str_replace(name,".txt.gz","")  #replace
    # tmp <- str_split_fixed(names(kmpre$CellLabel),pattern='-',n=2)
    tmp <- str_split_fixed(name,pattern='_',n=2)
    xQTL1 = tmp[1]
    xQTL2 = tmp[2]
    x_QTL <-paste(xQTL1," lambda", sep = "")
    y_QTL <-paste(xQTL2," lambda", sep = "")
    setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/xQTL_merge/Cis_10MB/")
    org<-fread(file,header = T,sep = "\t") %>% as.data.frame()
    n_r<-nrow(org)
    if (n_r >1){
        colnames(org)[1] <- "x_QTL"
        colnames(org)[2] <- "y_QTL" #change colname name 
        # round(org$x_QTL,4)
        figure_name <-name <- str_replace(file,".txt.gz",".png")
        fittings= c("lm","glm","loess")
        
        for(fitting in fittings){
            dir = paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/xQTL_merge/Cis_10MB_fitting/Cis_10MB_",fitting,sep="")
            title_name = paste("Cis 10MB", fitting,sep=" ")
            setwd(dir)
            png(figure_name)

            p1<-ggplot(org, aes(x =x_QTL, y=y_QTL)) +geom_point(size=0.1, color ="black" ) 
            p1<-p1+stat_smooth(method=fitting)
            p1<-p1+stat_cor(data=org, method = "pearson")
            p2<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                axis.line = element_line(colour = "black") )
            p2<-p2+xlab(x_QTL) + ylab(y_QTL)
            p2<-p2+ggtitle(title_name) 
            print(p2)
            dev.off()
            i=i+1
            print(i)
        }
        j=j+1
        print(j)
    }
}
