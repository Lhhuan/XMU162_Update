library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(parallel)
library(conflicted)
library(gridExtra)
library(Hmisc)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/script/ALL_eQTL/cis_trans/kmer/interval_18/")
group1s <-c("cis_1MB","cis_10MB")
group2s <-c("trans_1MB","trans_10MB")
figure_list<-list()
i=1
for(group in group1s){
    org<-read.table(paste0(group,"/communities.csv"),header =T,sep=",")%>%data.frame()
    a<-table(org$Group)
    # mycolor<-c("#f9c00c","#00b9f1","#7200da","#f9320c","#75D701","black")
    refine_group <-str_replace(group,"_"," ")
    refine_group <-capitalize(refine_group)
    # mycolor<-c("#efb08c","#d08752","#ffd369","#adce74","#5aa469","#16697a")
    mycolor<-c("#93abd3","#9ad3bc","#f3eac2","#f5b461","#ec524b","#fca3cc")
    # mycolor<-c( "#70A1D7","#C86B85","#FFD2A5", "#C79ECF", "#C1C0B9", "#A1DE93")
    #30E3CA #C86B85 #C79ECF #FFDED,E #A1DE93 #70A1D7
    b <-data.frame(a)
    colnames(b)[1] <-"community"
    pdf(paste(refine_group,".pdf",sep=""))
    pie(b$Freq, cex=2,col = mycolor,labels = b$Freq, radius = 1,main=refine_group,cex.main=2)
    legend("topright", c("0","1","2","3","4","5"), cex = 1, fill = mycolor)
    dev.off()
}

for(group in group2s){
    org<-read.table(paste0(group,"/communities.csv"),header =T,sep=",")%>%data.frame()
    a<-table(org$Group)
    # mycolor<-c("#f9c00c","#00b9f1","#7200da","#f9320c","#75D701","black")
    refine_group <-str_replace(group,"_"," ")
    refine_group <-capitalize(refine_group)
    # mycolor<-c("#efb08c","#d08752","#ffd369","#adce74","#5aa469","#16697a")
    # mycolor<-c("#93abd3","#9ad3bc","#f3eac2","#f5b461","#ec524b","#fca3cc")
    mycolor<-c( "#70A1D7","#C86B85","#FFD2A5", "#C79ECF", "#C1C0B9", "#A1DE93")
    #30E3CA #C86B85 #C79ECF #FFDED,E #A1DE93 #70A1D7
    b <-data.frame(a)
    colnames(b)[1] <-"community"
    pdf(paste(refine_group,".pdf",sep=""))
    pie(b$Freq, cex=2,col = mycolor,labels = b$Freq, radius = 1,main=refine_group,cex.main=2)
    legend("topright", c("0","1","2","3","4","5"), cex = 1, fill = mycolor)
    dev.off()
}