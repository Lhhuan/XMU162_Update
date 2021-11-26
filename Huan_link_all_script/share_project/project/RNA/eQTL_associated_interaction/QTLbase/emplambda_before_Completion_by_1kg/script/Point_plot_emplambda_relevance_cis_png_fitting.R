library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)


files<- list.files("~/project/RNA/eQTL_associated_interaction/QTLbase/output/xQTL_merge/Trans_1MB","txt.gz")
# n<-length(files)
i=1
j=1
for (file in files){
    name <- str_replace(file,"NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_trans_1MB_","")  #replace
    name <- str_replace(name,".txt.gz","")  #replace
    # tmp <- str_split_fixed(names(kmpre$CellLabel),pattern='-',n=2)
    tmp <- str_split_fixed(name,pattern='_',n=2)
    xQTL1 = tmp[1]
    xQTL2 = tmp[2]
    x_QTL <-paste(xQTL1," emplambda", sep = "")
    y_QTL <-paste(xQTL2," emplambda", sep = "")
    setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/xQTL_merge/Trans_1MB/")
    org<-fread(file,header = T,sep = "\t") %>% as.data.frame()
    n_r<-nrow(org)
    if (n_r >1){
        colnames(org)[1] <- "x_QTL"
        colnames(org)[2] <- "y_QTL" #change colname name 
        # round(org$x_QTL,4)
        figure_name <-name <- str_replace(file,".txt.gz",".png")
        # fittings= c("lm","glm","loess")
        fittings= c("lm","glm")
        for(fitting in fittings){
            dir = paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/xQTL_merge/Trans_1MB_fitting/Trans_1MB_",fitting,sep="")
            title_name = paste("Trans 1MB", fitting,sep=" ")
            setwd(dir)
            png(figure_name,height = 1000,width = 1000)

            p1<-ggplot(org, aes(x =x_QTL, y=y_QTL)) +geom_point(size=1, color ="black" ) 
            p1<-p1+stat_smooth(method=loess)
            p1<-p1+stat_smooth(method=fitting)
            p2<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20),axis.line = element_line(colour = "black"),
                axis.text.y = element_text(size=20),axis.text.x = element_text(size=20)) 
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

files<- list.files("~/project//RNA/eQTL_associated_interaction/QTLbase/output/xQTL_merge/Trans_10MB","txt.gz")
# n<-length(files)
i=1
j=1
for (file in files){
    name <- str_replace(file,"NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_trans_10MB_","")  #replace
    name <- str_replace(name,".txt.gz","")  #replace
    # tmp <- str_split_fixed(names(kmpre$CellLabel),pattern='-',n=2)
    tmp <- str_split_fixed(name,pattern='_',n=2)
    xQTL1 = tmp[1]
    xQTL2 = tmp[2]
    x_QTL <-paste(xQTL1," emplambda", sep = "")
    y_QTL <-paste(xQTL2," emplambda", sep = "")
    setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/xQTL_merge/Trans_10MB/")
    org<-fread(file,header = T,sep = "\t") %>% as.data.frame()
    n_r<-nrow(org)
    if (n_r >1){
        colnames(org)[1] <- "x_QTL"
        colnames(org)[2] <- "y_QTL" #change colname name 
        # round(org$x_QTL,4)
        figure_name <-name <- str_replace(file,".txt.gz",".png")
        # fittings= c("lm","glm","loess")
        fittings= c("lm","glm")
        for(fitting in fittings){
            dir = paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/xQTL_merge/Trans_10MB_fitting/Trans_10MB_",fitting,sep="")
            title_name = paste("Trans 10MB", fitting,sep=" ")
            setwd(dir)
            png(figure_name,height = 1000,width = 1000)

            p1<-ggplot(org, aes(x =x_QTL, y=y_QTL)) +geom_point(size=1, color ="black" ) 
            p1<-p1+stat_smooth(method=loess)
            p1<-p1+stat_smooth(method=fitting)
            p2<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20),axis.line = element_line(colour = "black"),
                axis.text.y = element_text(size=20),axis.text.x = element_text(size=20)) 
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
