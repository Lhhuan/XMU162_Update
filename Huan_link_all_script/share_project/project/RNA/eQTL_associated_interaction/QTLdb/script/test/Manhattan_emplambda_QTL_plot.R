library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(CMplot)

# xQTL = "eQTL"
QTLs = c("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL","QTL")
for (xQTL in QTLs){
    current_directory <- paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/ALL_",xQTL,"/", sep = "")
    setwd(current_directory)



    cutoff =24.25
    intervals =1000
    j=intervals 
    input_file <-paste("../../output/ALL_eQTL/NHPoisson_emplambda_interval_",j,"cutoff_",cutoff,"_all_",xQTL,".txt", sep = "")
    org<-read.table(input_file,header = T,sep = "\t") %>% as.data.frame()
    org1<-filter(org, emplambda!="NA")
    colnames(org1)[2] <- 'Position'
    colnames(org1)[3] <- 'Chromosome' #change colname name 
    SNP <-paste(org1$Position,"_",org1$Chromosome)
    org1$SNP <-SNP
    figure_name <-paste("cutoff_",cutoff,"_snp_number_in_different_interval_all_",xQTL, sep = "")
    CMplot(org1,plot.type="m",LOG10=FALSE,threshold=NULL,chr.den.col=NULL,file="pdf",memo="",dpi=300, ylab = "Emplambda",file.output=TRUE,verbose=TRUE)
}
