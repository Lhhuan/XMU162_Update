library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(data.table)
library(parallel)
library(R.utils)
# R >= 3.1.1 is needed

if (!require("NHPoisson")) {
    install.packages("NHPoisson")
    library("NHPoisson")
}

# setting a random seed for reproducibility
set.seed(1138)

POTevents.fun <- function (T, thres, date = NULL)
{
    if (is.null(date))
        date <- c(1:length(T))
    date <- as.matrix(date)
    if ((is.null(date) == FALSE) & (dim(date)[1] != length(T)))
        stop("T and date must have the same number of observations")
    exc <- (T > thres)
    inrachtx <- (diff(c(0, exc)) == 1)
    Pi <- c(1:length(T))[inrachtx == 1]
    numerachtx <- cumsum(inrachtx)[exc == 1]
    intentx <- (T - thres)[exc == 1]
    Im <- tapply(intentx, INDEX = numerachtx, FUN = mean)
    Ix <- tapply(intentx, INDEX = numerachtx, FUN = max)
    L <- tapply(intentx, INDEX = numerachtx, FUN = length)
    Px <- Pi + tapply(intentx, INDEX = numerachtx, FUN = which.max) - 1
    inddat <- 1 - exc
    inddat[Px] <- 1
    datePi <- date[Pi, ]
    datePx <- date[Px, ]
    cat("Number of events: ", length(Im), fill = TRUE)
    cat("Number of excesses over threshold", thres, ":", sum(exc),
        fill = TRUE)
    return(list(Pi = Pi, datePi = datePi, Px = Px, datePx = datePx,
        Im = Im, Ix = Ix, L = L, inddat = inddat, T = T, thres = thres,
        date = date))
}

#----mkdir and set pwd
#----------------

get_nhp <- function(j = NULL, org = NULL, cutoff = NULL, name = NULL){
    rs1 <- data.frame() #create the table
    rs2 <- data.frame()
    rs3 <- data.frame()
    for(i in c(1:22)){
        # ProcessBedGz <- function(i = NULL){
        org2<-filter(org, SNP_chr==i)
        org2$SNP_pos<-as.numeric(as.character(org2$SNP_pos))
        org2$SNP_chr<-as.numeric(as.character(org2$SNP_chr))
        #-------------------------------------- snp has unique p，min_p
        org_p <- org2%>%dplyr::select(SNP_pos,Pvalue)%>%unique()
        rm(org2)
        org1<-group_by(org_p, SNP_pos)%>%summarise(min_p = min(Pvalue))%>%as.data.frame()
        #-----------------------------------
        org1 <- org1[order(org1$SNP_pos), ]
        BarEv <- POTevents.fun(T = -log10(org1$min_p), thres = cutoff, date = org1$SNP_pos) #POTevents.fun
        tB <-org1$SNP_pos
        title <-paste ("int",j,"bp_chr",i ,sep = "")
        # figure_name <-paste("emplambda_all_eQTL_int_",j,"bp_chr",i,".png", sep = "")
        # png(figure_name,height = 800,width = 1000)
        emplambdaB <- emplambda.fun(posE = BarEv$Px,inddat=BarEv$inddat, t = tB, lint = j, tit = title,plotEmp=FALSE) #emplambda.fun  #BarEv$Px 为校正后的矩阵，只保留了连续显著的最显著的点
        # dev.off()
        tmp1<-cbind(Px=BarEv$Px,datePx=BarEv$datePx,chr=i)%>% as.data.frame()
        tmp3<-cbind(inddat=BarEv$inddat,T=BarEv$T,chr=i)%>% as.data.frame()
        tmp2<-cbind(emplambda=emplambdaB$emplambda, t=emplambdaB$t,chr=i) %>% as.data.frame()
        rs1 <- bind_rows(rs1,tmp1)
        rs2 <- bind_rows(rs2,tmp2)
        rs3 <- bind_rows(rs3,tmp3)
        print(i)
    }
    # print(file)
    output_file1<-paste("NHPoisson_POTevents_Px_cutoff_7.3_",name,".txt", sep = "")
    output_file2<-paste("NHPoisson_emplambda_interval_",j,"_cutoff_7.3_",name,".txt", sep = "")
    output_file3<-paste("NHPoisson_POTevents_inddat_cutoff_7.3_",name,".txt", sep = "")
    write.table(rs1,output_file1,row.names = F, col.names = T,quote =F,sep="\t")
    write.table(rs2,output_file2,row.names = F, col.names = T,quote =F,sep="\t")
    write.table(rs3,output_file3,row.names = F, col.names = T,quote =F,sep="\t")
    #-----------------------------
    gzip(output_file1)
    gzip(output_file2)
    gzip(output_file3)  
    # print(j)
}

#----------------------------------------------

tissue= "Tissue_merge"
print(c(tissue,"start"))
pwd_dir <-paste0("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/",tissue,"/Cis_eQTL/NHP/")
dir.create(pwd_dir,recursive = T)
setwd(pwd_dir)
org <- read.table(paste0("../01_merge_all_tissue_cis_eQTL_1kg_Completion.txt.gz"),header = F,sep = "\t") %>% as.data.frame()
colnames(org) <- c("SNP_chr","SNP_pos","Pvalue")
cutoff =-log10(5e-8)
num = 18
get_nhp(j = num, org = org, cutoff = cutoff, name = tissue)
print(c(tissue,"end"))

# print $O1 "SNP_chr\tSNP_pos\tPvalue\n";
