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


# xQTL = "eQTL"
# QTLs = c("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL","QTL")
QTLs = c("eQTL","hQTL","sQTL","QTL","mQTL","caQTL")
n<-length(QTLs)
# for (xQTL in QTLs)
fun=function(i){
    xQTL =QTLs[i]
    current_directory <- paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/ALL_",xQTL,"/Manhattan/", sep = "")
    if(dir.exists(current_directory)){
    print("dir_name exist")
    } else{
        dir.create(current_directory)
    }
    setwd(current_directory)

    # cutoff =24.25
    cutoff =7.3
    intervals =1000
    j=intervals 
    input_file <-paste("../../../output/ALL_", xQTL, "/NHPoisson_emplambda_interval_",j,"cutoff_",cutoff,"_all_",xQTL,".txt", sep = "")
    org<-fread(input_file,header = T,sep = "\t") %>% as.data.frame()
    org1<-filter(org, !is.na(emplambda))
    colnames(org1)[2] <- 'Position'
    colnames(org1)[3] <- 'Chromosome' #change colname name 
    SNP <-paste(org1$Position,"_",org1$Chromosome)
    org1$SNP <-SNP
    org1=org1[,c(4,3,2,1)]
    unique_chr<-unique(org1$Chromosome)
    n_chr<-length(unique_chr)
    for (j in c(1:n_chr)){
        chr = unique_chr[j]
        print(chr)
        chr_file<-filter(org1, Chromosome==chr)
        figure_name <-paste("cutoff_",cutoff,"_snp_number_in_different_interval_all_",xQTL,"_chr_",chr,"_cex_0.08",sep = "")
        CMplot(chr_file,plot.type="m",LOG10=FALSE,threshold=NULL,chr.den.col=NULL,file="jpg",memo=figure_name,dpi=300, cex=0.08,ylab = "Emplambda",file.output=TRUE,verbose=TRUE)
    }

    # return(dim(org))
}

#---------------------------------parallel
library(parallel)
a=mclapply(1:n,function(x){fun(x)},mc.cores=6)


# library(doParallel)
# library(foreach)
# cores=detectCores() #检查可用核数
# # cl <- makeCluster(cores-1)
# cl <- makeCluster(9)
# registerDoParallel(cl)
# #--------------------------------
# cl_library = function(){
#     library(ggplot2)
#     library(dplyr)
#     library(Rcpp)
#     library(readxl)
#     library(stringr)
#     library(gcookbook)
#     library(gridExtra)
#     library(ggpubr)
#     library(CMplot)
#     library(data.table)
# }
# # clusterCall(cl, function() library(ggplot2)) #特殊需要
# # clusterCall(cl, function() library(dplyr))
# # clusterCall(cl, function() library(Rcpp))
# # clusterCall(cl, function() library(readxl))
# # clusterCall(cl, function() library(stringr))
# # clusterCall(cl, function() library(gcookbook))
# # clusterCall(cl, function() library(gridExtra))
# # clusterCall(cl, function() library(ggpubr))
# # clusterCall(cl, function() library(CMplot))
# clusterCall(cl, cl_library)

# #----------------------------------------
# result1 = foreach(i=1:n) %dopar% fun(i)
# # result1=foreach(j=1:n,.combine='r') %dopar% function(j)
# stopImplicitCluster(cl)