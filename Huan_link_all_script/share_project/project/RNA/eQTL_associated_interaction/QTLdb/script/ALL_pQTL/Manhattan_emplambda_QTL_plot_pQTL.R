library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(CMplot)
library(data.table)



xQTL = "pQTL"
current_directory <- paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/ALL_",xQTL,"/", sep = "")
setwd(current_directory)

# cutoff =24.25
cutoff =7.3
intervals =1000
j=intervals 
input_file <-paste("../../output/ALL_", xQTL, "/NHPoisson_emplambda_interval_",j,"cutoff_",cutoff,"_all_",xQTL,".txt", sep = "")
org<-fread(input_file,header = T,sep = "\t") %>% as.data.frame()
org2<-filter(org, !is.na(emplambda))
org1 <- filter(org2,emplambda < 1) # discard outlier
colnames(org1)[2] <- 'Position'
colnames(org1)[3] <- 'Chromosome' #change colname name 
SNP <-paste(org1$Position,"_",org1$Chromosome)
org1$SNP <-SNP
org1=org1[,c(4,3,2,1)]
figure_name <-paste("cutoff_",cutoff,"_snp_number_in_different_interval_all_",xQTL,"_cex_0.3_single",sep = "")
CMplot(org1,plot.type="m",LOG10=FALSE,threshold=NULL,chr.den.col=NULL,file="jpg",memo=figure_name,dpi=300, cex=0.3,ylab = "Emplambda",file.output=TRUE,verbose=TRUE)
