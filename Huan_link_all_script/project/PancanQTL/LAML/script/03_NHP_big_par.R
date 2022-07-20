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
        org_pg <- group_by(org_p, SNP_pos)
        org1<-summarise(org_pg,min_p = min(Pvalue))%>%as.data.frame()
        rm(org_pg)
        #-----------------------------------
        org1 <- org1[order(org1$SNP_pos), ]
        BarEv <- POTevents.fun(T = -log10(org1$min_p), thres = cutoff, date = org1$SNP_pos) #POTevents.fun
        tB <-org1$SNP_pos
        title <-paste ("int",j,"bp_chr",i ,sep = "")
        # figure_name <-paste("emplambda_all_eQTL_int_",j,"bp_chr",i,".png", sep = "")
        # png(figure_name,height = 800,width = 1000)
        emplambdaB <- emplambda.fun(posE = BarEv$Px,inddat=BarEv$inddat, t = tB, lint = j, tit = title,plotEmp=FALSE) #emplambda.fun
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
    output_file1<-paste("NHPoisson_POTevents_Px_cutoff_",cutoff,"_",name,".txt", sep = "")
    output_file2<-paste("NHPoisson_emplambda_interval_",j,"_cutoff_",cutoff,"_",name,".txt", sep = "")
    output_file3<-paste("NHPoisson_POTevents_inddat_cutoff_",cutoff,"_",name,".txt", sep = "")
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


# for(tissue in tissues){
ProcessBedGz <- function(tissue = NULL){
    print(c(tissue,"start"))
    pwd_dir <-"/home/huanhuan/project/PancanQTL/LAML/output/NHP"
    dir.create(pwd_dir,recursive = T)
    setwd(pwd_dir)
    org <- read.table("../LAML_cis_eQTL_1kg_Completion.txt.gz",header = T,sep = "\t") %>% as.data.frame()
    cutoff =-log10(5e-8)
    num = 18
    get_nhp(j = num, org = org, cutoff = cutoff, name = tissue)
    print(c(tissue,"end"))
}

# tissues <- c("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Cells_EBV-transformed_lymphocytes","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Kidney_Cortex","Muscle_Skeletal","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Uterus","Prostate")
tissue <- "LAML"

lapply(tissue,ProcessBedGz)
