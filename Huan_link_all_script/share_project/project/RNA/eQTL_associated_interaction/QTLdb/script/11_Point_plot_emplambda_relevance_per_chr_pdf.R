library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)


files<- list.files("~/project//RNA/eQTL_associated_interaction/QTLbase/output/xQTL_merge/ALL/","txt.gz")
# n<-length(files)
i=1
for (file in files){
    name <- str_replace(file,"NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_","")  #replace
    name <- str_replace(name,".txt.gz","")  #replace
    # tmp <- str_split_fixed(names(kmpre$CellLabel),pattern='-',n=2)
    tmp <- str_split_fixed(name,pattern='_',n=2)
    xQTL1 = tmp[1]
    xQTL2 = tmp[2]
    x_QTL <-paste(xQTL1," emplambda", sep = "")
    y_QTL <-paste(xQTL2," emplambda", sep = "")
    setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/xQTL_merge/ALL/")
    org<-fread(file,header = T,sep = "\t") %>% as.data.frame()
    n_r<-nrow(org)
    if (n_r >1){
        colnames(org)[1] <- "x_QTL"
        colnames(org)[2] <- "y_QTL" #change colname name 
        # round(org$x_QTL,4)
        n_chr<-unique(org$chr)
        for (j in n_chr){
            pwd_dir =paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/xQTL_merge/Chr",j,"/",sep="")
            if(dir.exists(pwd_dir)){
                print("dir_name exist")
            } else{
                dir.create(pwd_dir)
            }
            setwd(pwd_dir)
            org_chr<-filter(org,chr==j)
            figure_name <-name <- str_replace(file,".txt.gz",".pdf")
            title<-paste("Chr",j,sep="")
            pdf(figure_name,height = 7,width = 7)
            p1<-ggplot(org_chr, aes(x =x_QTL, y=y_QTL)) +geom_point(size=0.01, color ="black" ) 

            p2<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                axis.title.y = element_text(size = 10),axis.title.x = element_text(size = 10),axis.line = element_line(colour = "black"),
                axis.text.y = element_text(size=10),axis.text.x = element_text(size=10)) 
            p2<-p2+xlab(x_QTL) + ylab(y_QTL) 
            p2<-p2+ggtitle(title)
            print(p2)
            dev.off()
            print("j=")
            print(j)
        }
        #--------------------------------------
        # figure_name <-name <- str_replace(file,".txt.gz",".png")
        # png(figure_name,height = 1000,width = 1000)

        # p1<-ggplot(org, aes(x =x_QTL, y=y_QTL)) +geom_point(size=0.1, color ="black" ) 

        # p2<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        #     axis.title.y = element_text(size = 10),axis.title.x = element_text(size = 10),axis.line = element_line(colour = "black"),
        #     axis.text.y = element_text(size=10),axis.text.x = element_text(size=10)) 
        # p2<-p2+xlab(x_QTL) + ylab(y_QTL) 
        # print(p2)
        # dev.off()
        i=i+1
    }
    print(i)
}
