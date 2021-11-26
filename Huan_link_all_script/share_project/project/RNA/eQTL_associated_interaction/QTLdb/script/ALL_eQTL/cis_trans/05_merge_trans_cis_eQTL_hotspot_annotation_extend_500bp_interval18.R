library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(parallel)
library(conflicted)

interval = 18
ProcessBedGz <- function(file = NULL){
# for (file in file2){
  # name <- str_replace(file,paste("_interval_18_cutoff_7.3_",type,"_eQTL_segment_",group,".bed.gz"),"")  #replace
  name <- str_replace(file,"_interval_18_cutoff_7.3_","")  #replace #-------
  name <- str_replace(name,"cis_1MB_eQTL_segment_hotspot.bed.gz","")  #replace
  name <- str_replace(name,"cis_10MB_eQTL_segment_hotspot.bed.gz","")  #replace
  name <- str_replace(name,"trans_10MB_eQTL_segment_hotspot.bed.gz","")  #replace
  name <- str_replace(name,"trans_1MB_eQTL_segment_hotspot.bed.gz","")  #replace
  name <- str_replace(name,"trans_1MB_eQTL_segment_non_hotspot.bed.gz","")  #replace
  name <- str_replace(name,"trans_10MB_eQTL_segment_non_hotspot.bed.gz","")  #replace
  name <- str_replace(name,"cis_10MB_eQTL_segment_non_hotspot.bed.gz","")  #replace
  name <- str_replace(name,"cis_1MB_eQTL_segment_non_hotspot.bed.gz","")  #replace
  # name <- str_replace(name,"_narrow_peak","")  #replace
  # name <- str_replace(name,"_loops_1_2","")  #replace
  # print(file)
  print(name)
  a<-file.info(file)$size
  if (a > 20){ #the null compressed files is 20
    if (name != "hic"){
      print("klk")
      org2<-read.table(file,header = F,sep = "\t") %>% as.data.frame()
      org2$name <-paste(org2$V4,org2$V5,org2$V6,sep=":")
      # colnames(org2)[9]<-name
      org2$key1 <-paste(org2$V1,org2$V2,org2$V3,sep=":")
      org2<-org2[!duplicated(org2$key1),]
      rownames(org2)<-org2$key1
      # org2<-org2[,-c(1:6)]
      org2<-org2[,-c(1:8)]
      colnames(org2)[1]<-name
      #--------------------
    } else { #hi-C
      org2<-read.table(file,header = F,sep = "\t") %>% as.data.frame()
      # org2<-org2[,-c(7:10)]
      org2<-org2[,-c(7:12)]
      org2$name <-paste(org2$V4,org2$V5,org2$V6,sep=":")
      # colnames(org2)[9]<-name
      org2$key1 <-paste(org2$V1,org2$V2,org2$V3,sep=":")
      org2<-org2[!duplicated(org2$key1),]
      rownames(org2)<-org2$key1
      org2<-org2[,-c(1:6)]
      colnames(org2)[1]<-name
    }
  }else{ #a<=20,file is 0 line
    print(paste(name,"NULL",sep="\t"))
    org2<-data.frame(a=c(1:10),key1=c(11:20)) #----------
    rownames(org2)<-org2$key1
    colnames(org2)[1]<-name
  }
  gc() #remove space
  print("finish")
  return(org2)
}


intervals = c(6,7,8,9,12,15,18)
groups = c("hotspot","non_hotspot")
types= c("cis_1MB","cis_10MB","trans_1MB","trans_10MB")

# # fractions = c(1,0.9,0.8,0.7,0.6,0.5)
# fractions = c(0.1,1)
# groups = c("hotspot","non_hotspot")
# types= c("trans_1MB","trans_10MB")

# fractions = c(1,0.9,0.8,0.7,0.6,0.5)
fractions = c(0.1,1)
tryCatch({
  for (fraction in fractions){
    for(group in groups){
      #--------------------
      mclapply(types, function(type){
        setwd(paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/output/ALL_eQTL/cis_trans",group,sep="/")) #-----------background
        org1<-read.table(paste("interval_",interval,"_cutoff_7.3_",type,"_eQTL_segment_",group,".bed.gz",sep=""),header = F,sep = "\t") %>% as.data.frame()
        colnames(org1)[1]<-"chr"
        colnames(org1)[2]<-"start"
        colnames(org1)[3]<-"end"
        org1$key <-paste(org1$chr,org1$start,org1$end,sep = ":")
        tmp<-org1 #9904240
        # tmp<-tmp[,-c(1:4)]
        rownames(tmp)<-tmp$key
        tmp<-tmp%>%dplyr::select(-key) 
        #-------------------------------------
        setwd(paste0("/state/partition1/huan/QTLdb/output/ALL_eQTL/cis_trans/fisher_exact_test_extend_500bp/annotation/interval_",interval,"/",group,"/",type,"/",fraction,"/"))


        files<- list.files(paste0("/state/partition1/huan/QTLdb/output/ALL_eQTL/cis_trans/fisher_exact_test_extend_500bp/annotation/interval_",interval,"/",group,"/",type,"/",fraction,"/"),"bed.gz")

        #--------------- #sub function
        org2 <- mclapply(files, function(file){ProcessBedGz(file=file)}, mc.cores = 10)
        #-------------------
        for(i in 1:length(org2)){
          tmp<-cbind(tmp,org2[[i]][rownames(tmp),])
          j=i+3
          tmp[,j][!is.na(tmp[,j])]=1
          tmp[,j][is.na(tmp[,j])]=0
          # tmp<-tmp[,-length(tmp)]
          tmp <-tmp%>%dplyr::select(-key1)
          print(i)
        }

        gz1 <- gzfile(paste("/state/partition1/huan/QTLdb/output/ALL_eQTL/cis_trans/fisher_exact_test_extend_500bp/annotation_out/",group,"/05_annotation_merge_",type,"_interval_",interval,"_overlap_fraction_",fraction,".txt.gz",sep = ""), "w")
        write.table(tmp,gz1,row.names = F, col.names = T,quote =F,sep="\t")
        close(gz1)
        print("finish_org2")
        #-----------------
        }, mc.cores = 2) 
        #----------------------
        print(group)
        print(fraction)
        out <-paste(group,fraction,sep=" ")
        print(out)
    }
  }
},warning=function(w){
  # a<-paste("warning",group,type,fraction,sep=" ")
  print("warning")
},error=function(e){
  # bbb<-paste("error",group,type,fraction,sep=" ")
  print("error")
})