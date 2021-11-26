library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

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

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/ALL_trans_eQTL/10MB")

# org<-read.table("../../../output/merge_QTL_all_QTLtype_pop.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org<-read.table("../../../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org_pop <- filter(org,QTL_type == "eQTL" && cis_trans_10MB == "trans") #10MB
# org_pop <- org
#--------------------500k
interval1<-seq(from=2000, to=8000, by=2000)
interval2<-seq(from=10000, to=100000, by=10000)
interval=c(interval1,interval2) #两个interval合并
cutoff = 7
for (j in interval){
    dir_name = paste("all_10MB_trans_eQTL_cutoff_",cutoff,"_int",j, sep = "")
    #----------------------------------mkdir_dir
    if(dir.exists(dir_name)){
        print("dir_name exist")
    }
    else{
        dir.create(dir_name)
    }
    rs1 <- data.frame() #create the table
    rs2 <- data.frame()
    rs3 <- data.frame()
    for(i in c(1:22)){
        org2<-filter(org_pop, SNP_chr==i)
        #-------------------------------------- snp has unique p
        org_p <- org2%>%dplyr::select(SNP_pos,Pvalue)%>%unique()
        org_pg <- group_by(org_p, SNP_pos)
        org1<-summarise(org_pg,min_p = min(Pvalue))%>%as.data.frame()
        #-----------------------------------
        org1 <- org1[order(org1$SNP_pos), ]
        BarEv <- POTevents.fun(T = -log10(org1$min_p), thres = cutoff, date = org1$SNP_pos) #POTevents.fun
        tB <-org1$SNP_pos
        title <-paste ("int",j,"bp_chr",i ,sep = "")
        figure_name <-paste(dir_name,"/emplambda_all_trans_10MB_eQTL_int_",j,"bp_chr",i,".png", sep = "")
        png(figure_name,height = 800,width = 1000)
        emplambdaB <- emplambda.fun(posE = BarEv$Px,inddat=BarEv$inddat, t = tB, lint = j, tit = title, plotEmp=FALSE) #emplambda.fun
        dev.off()
        tmp1<-cbind(Px=BarEv$Px,datePx=BarEv$datePx,chr=i)%>% as.data.frame()
        tmp3<-cbind(inddat=BarEv$inddat,T=BarEv$T,chr=i)%>% as.data.frame()
        tmp2<-cbind(emplambda=emplambdaB$emplambda, t=emplambdaB$t,chr=i) %>% as.data.frame()
        rs1 <- bind_rows(rs1,tmp1)
        rs2 <- bind_rows(rs2,tmp2)
        rs3 <- bind_rows(rs3,tmp1)
        print(i)
    }
    output_file1<-paste("../../../output/ALL_trans_eQTL/10MB/NHPoisson_POTevents_Px_cutoff_",cutoff,"_all_trans_eQTL.txt", sep = "")
    output_file2<-paste("../../../output/ALL_trans_eQTL/10MB/NHPoisson_emplambda_interval_",j,"_all_trans_eQTL.txt", sep = "")
    output_file3<-paste("../../../output/ALL_trans_eQTL/10MB/NHPoisson_POTevents_inddat_cutoff_",cutoff,"_all_trans_eQTL.txt", sep = "")
    write.table(rs1,output_file1,row.names = F, col.names = T,quote =F,sep="\t")
    write.table(rs2,output_file2,row.names = F, col.names = T,quote =F,sep="\t")
    write.table(rs3,output_file3,row.names = F, col.names = T,quote =F,sep="\t")
}
#----------------------------------------------
