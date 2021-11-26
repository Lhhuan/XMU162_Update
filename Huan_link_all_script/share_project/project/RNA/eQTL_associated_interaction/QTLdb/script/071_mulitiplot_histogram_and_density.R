library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.title.y = element_text(size = 15),axis.title.x = element_text(size = 15),axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10),axis.text.x = element_text(size=10))

qtls = c("eQTL","riboQTL","reQTL","miQTL","metaQTL","lncRNAQTL","edQTL","cerQTL","pQTL","caQTL","sQTL","hQTL","mQTL")

for(qtl in qtls){
        dir<-paste("~/project/RNA/eQTL_associated_interaction/QTLbase/figure/ALL_",qtl,sep="")
        setwd(dir)
        file_name<-paste("../../output/ALL_",qtl,"/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_",qtl,".txt.gz",sep="")
        org<-fread(file_name,header=T,sep="\t")%>%as.data.frame()
        top_title<-qtl
        histogram_name<-paste("histogram_interval_1000_cutoff_7.3_",qtl,".pdf",sep="")
        pdf(histogram_name,height = 6,width = 7)
        p1<-ggplot(org, aes(x =emplambda)) +geom_histogram(position="identity")
        p1<-p1+xlab("Emplambda") + ylab("Count")+p_theme
        p1<-p1+ggtitle(top_title)+theme(plot.title = element_text(hjust = 0.5))
        print(p1)
        dev.off()
        #----------------------
        density_name<-paste("Density_interval_1000_cutoff_7.3_",qtl,".pdf",sep="")
        pdf(density_name,height = 6,width = 7)
        p1 <- ggplot(org, aes(x =emplambda)) + geom_density()
        p1<-p1+xlab("Emplambda") +ylab("Density")+p_theme
        p1<-p1+ggtitle(top_title)+theme(plot.title = element_text(hjust = 0.5))
        print(p1)
        dev.off()
        print(qtl)
}

types =c("trans_10MB","trans_1MB","cis_10MB","cis_1MB")

for(qtl in qtls){
        for(type in types){
                dir<-paste("~/project/RNA/eQTL_associated_interaction/QTLbase/figure/ALL_",qtl,sep="")
                setwd(dir)
                file_name<-paste("../../output/ALL_",qtl,"/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_",type,"_",qtl,".txt.gz",sep="")
                org<-fread(file_name,header=T,sep="\t")%>%as.data.frame()
                x<-nrow(org)
                if (x>0){
                        # print(file_name)
                        top_title<-paste(qtl,type,sep=": ")
                        top_title<-str_replace(top_title,"_"," ")
                        histogram_name<-paste("histogram_interval_1000_cutoff_7.3_",type,"_",qtl,".pdf",sep="")
                        pdf(histogram_name,height = 6,width = 7)
                        p1<-ggplot(org, aes(x =emplambda)) +geom_histogram(position="identity")
                        p1<-p1+xlab("Emplambda") + ylab("Count")+p_theme
                        p1<-p1+ggtitle(top_title)+theme(plot.title = element_text(hjust = 0.5))
                        print(p1)
                        dev.off()
                        #----------------------
                        density_name<-paste("Density_interval_1000_cutoff_7.3_",type,"_",qtl,".pdf",sep="")
                        pdf(density_name,height = 6,width = 7)
                        p1 <- ggplot(org, aes(x =emplambda)) + geom_density()
                        p1<-p1+xlab("Emplambda") +ylab("Density")+p_theme
                        p1<-p1+ggtitle(top_title)+theme(plot.title = element_text(hjust = 0.5))
                        print(p1)
                        dev.off()
                        print(qtl)
                }
        }
}
