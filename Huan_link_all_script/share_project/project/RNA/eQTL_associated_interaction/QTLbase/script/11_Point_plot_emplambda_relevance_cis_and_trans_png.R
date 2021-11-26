library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(Hmisc) #capitalize
library(ggpubr)#pearson

types=c("Cis_1MB","Cis_10MB","Trans_1MB","Trans_10MB")
for(type in types){
    files_dir<-paste("~/project/RNA/eQTL_associated_interaction/QTLbase/output/xQTL_merge/",type,"/",sep="")
    figure_dir <- paste("~/project/RNA/eQTL_associated_interaction/QTLbase/figure/xQTL_merge/",type,"/",sep="")
    files<- list.files(files_dir,"txt.gz")
    # n<-length(files)
    i=0
    print(type)
    for (file in files){
        name <- str_replace(file,"NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_","")  #replace
        name <- str_replace(name,".txt.gz","")  #replace
        # tmp <- str_split_fixed(names(kmpre$CellLabel),pattern='-',n=2)
        tmp <- str_split_fixed(name,pattern='_',n=4)
        c_o_t = tmp[1]
        distance_1 = tmp[2]
        xQTL1 = tmp[3]
        xQTL2 = tmp[4]
        x_QTL <-paste(xQTL1,"lambda", sep = " ")
        y_QTL <-paste(xQTL2,"lambda", sep = " ")
        setwd(files_dir)
        org<-fread(file,header = T,sep = "\t") %>% as.data.frame()
        n_r<-nrow(org)
        if (n_r >1){
            #-----------------------------------
            f1_name<-paste("../../../output/ALL_",xQTL1,"/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_",xQTL1,".txt.gz",sep="")
            f1<-fread(f1_name,header = T,sep = "\t") %>% as.data.frame()
            f2<-na.omit(f1)
            max_x =max(f2$emplambda)
            #-------------------------------------
            colnames(org)[1] <- "x_QTL"
            colnames(org)[2] <- "y_QTL" #change colname name 
            colnames(org)[5]<-"type"
            tmp2<-str_split_fixed(org$type,pattern='_',n=2)
            org$distance<-tmp2[,2]
            org$final_type<-capitalize(org$type)
            org1<-org[,-5]
            org2<-filter(org1, distance==distance_1)
            setwd(figure_dir)
            figure_name <-name <- str_replace(file,".txt.gz",".png")
            title_name <- paste(xQTL1,"type:",c_o_t, distance_1,sep=" ")
            labs_name <- paste(xQTL2,"type",sep=" ")
            png(figure_name)
            p1<-ggplot(org2, aes(x =x_QTL, y=y_QTL)) +geom_point(size=0.5,shape=1,aes(colour = factor(final_type))) 
            p2<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                axis.line = element_line(colour = "black")) 
            #-------------------------------
            # p2<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
            #     axis.title.y = element_text(size = 15),axis.title.x = element_text(size = 15),axis.line = element_line(colour = "black"),
            #     axis.text.y = element_text(size=10),axis.text.x = element_text(size=10)) 
            p2<-p2+xlab(x_QTL) + ylab(y_QTL) 
            p2<-p2+ggtitle(title_name)+theme(plot.title = element_text(hjust = 0.5))
            # p2<-p2+guides(fill=guide_legend(title=NULL))
            # p2<-theme(legend.title=element_text())
            p2<-p2+labs(colour=labs_name) #legend name
            # p2<-p2+scale_fill_discrete(name=labs_name)
            # p2<-p2+guides(colour=guide_legend(title=labs_name))
            p2<-p2+theme(legend.position="top") #legend.position
            p2<-p2+scale_x_continuous(limits=c(0, max_x))
            p2<-p2+stat_cor(data=org, method = "pearson")
            print(p2)
            dev.off()
            i=i+1
        }
        print("finish")
        cat(xQTL1, xQTL2, "\n",sep=" ")
        # print(xQTL1)
        print(i)
    }
}