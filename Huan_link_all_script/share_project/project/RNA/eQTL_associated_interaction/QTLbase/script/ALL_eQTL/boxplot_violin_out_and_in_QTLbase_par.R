library(ggplot2)
library(Rcpp)
# library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gcookbook)
library(gridExtra)
library(parallel)

#---------------------------------pdf
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/ALL_eQTL/")


xQTL = "eQTL"
# figure_list <-list ()
# for(i in c(6:32)){
ProcessBedGz <- function(i = NULL){  
  file1<-paste("../../output/ALL_",xQTL,"/QTLbase_NHPoisson_emplambda_interval_",i,"_cutoff_7.3_all_",xQTL,".txt.gz",sep ="")
  file2<-paste("../../output/ALL_",xQTL,"/1kg_out_QTLbase_NHPoisson_emplambda_interval_",i,"_cutoff_7.3_all_",xQTL,".txt.gz",sep ="")

  org1<-read.table(file1,header = T,sep = "\t") %>% as.data.frame()
  org2<-read.table(file2,header = T,sep = "\t") %>% as.data.frame()


  org1$Class <-"YES"
  org2$Class <- "NO"

  rs<- bind_rows(org1,org2)
  top_title<-paste(xQTL,"interval",i,sep = " ")

  figure_name <-paste("./diff/boxplot_violin_plot_cutoff_7.3_interval_",i,"_",xQTL,"_in_and_out_QTLbase.pdf", sep = "")

  pdf(figure_name ,height = 4,width = 5) 
  p <- ggplot(rs,aes(x=Class,y=emplambda,fill=Class))+geom_violin(width=0.8)+geom_boxplot(width=0.08)+
    scale_y_continuous(name = "lambda") + scale_x_discrete(name = "Class")+
    ggtitle(top_title)+theme(plot.title = element_text(hjust = 0.5))

  print(p)
  dev.off()
  print("aaa")
  reture(i)
}

number_list<-c(6:32)
aaaaa <- mclapply(number_list, ProcessBedGz, mc.cores = 40)

hist_figure_name1 <-paste("boxplot_violin_plot_cutoff_7.3_interval_6_32_all_",xQTL,"in_and_out_QTLbase.pdf", sep = "")
pdf(hist_figure_name1,height = 7.5,width = 15) 
p1<-marrangeGrob(figure_list,ncol=3,nrow=2)  
print(p1)
dev.off()

#-----------------------------------------
#--------------------png
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/ALL_eQTL/")
xQTL = "eQTL"
# figure_list <-list ()
# for(i in c(6:32)){
ProcessBedGz <- function(i = NULL){  
  file1<-paste("../../output/ALL_",xQTL,"/QTLbase_NHPoisson_emplambda_interval_",i,"_cutoff_7.3_all_",xQTL,".txt.gz",sep ="")
  file2<-paste("../../output/ALL_",xQTL,"/1kg_out_QTLbase_NHPoisson_emplambda_interval_",i,"_cutoff_7.3_all_",xQTL,".txt.gz",sep ="")

  org1<-read.table(file1,header = T,sep = "\t") %>% as.data.frame()
  org2<-read.table(file2,header = T,sep = "\t") %>% as.data.frame()


  org1$Class <-"YES"
  org2$Class <- "NO"

  rs<- bind_rows(org1,org2)
  top_title<-paste(xQTL,"interval",i,sep = " ")

  figure_name <-paste("./diff/boxplot_violin_plot_cutoff_7.3_interval_",i,"_",xQTL,"_in_and_out_QTLbase.png", sep = "")

  png(figure_name) 
  p <- ggplot(rs,aes(x=Class,y=emplambda,fill=Class))+geom_violin(width=0.8)+geom_boxplot(width=0.08)+
    scale_y_continuous(name = "lambda") + scale_x_discrete(name = "Class")+
    ggtitle(top_title)+theme(plot.title = element_text(hjust = 0.5))

  print(p)
  dev.off()
  print("aaa")
  reture(i)
}
number_list<-c(6:16)
aaaaa <- mclapply(number_list, ProcessBedGz, mc.cores = 20)