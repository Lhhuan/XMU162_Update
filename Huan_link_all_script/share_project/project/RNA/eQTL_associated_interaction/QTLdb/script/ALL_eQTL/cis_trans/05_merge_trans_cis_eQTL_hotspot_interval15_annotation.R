library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(parallel)
library(conflicted)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/hotspot/")
org1<-read.table("interval_15_cutoff_7.3_cis_1MB_eQTL_segment_hotspot.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(org1)[1]<-"chr"
colnames(org1)[2]<-"start"
colnames(org1)[3]<-"end"

org1$key <-paste(org1$chr,org1$start,org1$end,sep = ":")
tmp<-org1 #9904240
# tmp<-tmp[,-c(1:4)]
rownames(tmp)<-tmp$key
tmp<-tmp%>%dplyr::select(-key)
# rm(org1)
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/annotation/hotspot/cis_1MB")

files<- list.files("~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/annotation/hotspot/cis_1MB","bed.gz")
ProcessBedGz <- function(file = NULL){
# for (file in file2){
  name <- str_replace(file,"_interval_15_cutoff_7.3_cis_1MB_eQTL_segment_hotspot.bed.gz","")  #replace
  # name <- str_replace(name,".bed.gz","")  #replace
  # name <- str_replace(name,"_narrow_peak","")  #replace
  # name <- str_replace(name,"_loops_1_2","")  #replace
  print(name)
  if (name != "hic"){
    print("klk")
    org2<-read.table(file,header = F,sep = "\t") %>% as.data.frame()
    org2$name <-paste(org2$V4,org2$V5,org2$V6,sep=":")
    # colnames(org2)[9]<-name
    org2$key1 <-paste(org2$V1,org2$V2,org2$V3,sep=":")
    org2<-org2[!duplicated(org2$key1),]
    rownames(org2)<-org2$key1
    org2<-org2[,-c(1:6)]
    colnames(org2)[1]<-name
    #--------------------
  } else { #hi-C
    org2<-read.table(file,header = F,sep = "\t") %>% as.data.frame()
    org2<-org2[,-c(7:10)]
    org2$name <-paste(org2$V4,org2$V5,org2$V6,sep=":")
    # colnames(org2)[9]<-name
    org2$key1 <-paste(org2$V1,org2$V2,org2$V3,sep=":")
    org2<-org2[!duplicated(org2$key1),]
    rownames(org2)<-org2$key1
    org2<-org2[,-c(1:6)]
    colnames(org2)[1]<-name
  }
  gc() #remove space
  print("finish")
  return(org2)
}

org2 <- mclapply(files, ProcessBedGz, mc.cores = 20)

for(i in 1:length(org2)){
  tmp<-cbind(tmp,org2[[i]][rownames(tmp),])
  j=i+3
  tmp[,j][!is.na(tmp[,j])]=1
  tmp[,j][is.na(tmp[,j])]=0
  # tmp<-tmp[,-length(tmp)]
  tmp <-tmp%>%dplyr::select(-key1)
  print(i)
}
print("aaa")



for(i in 8:length(org2)){
  print("start")
  tmp<-cbind(tmp,org2[[i]][rownames(tmp),])
  j=i+3
  tmp[,j][!is.na(tmp[,j])]=1
  tmp[,j][is.na(tmp[,j])]=0
  # tmp<-tmp[,-length(tmp)]
  tmp <-tmp%>%dplyr::select(-key1)
  print(i)
}
print("aaa")
tmp[,3][!is.na(tmp[,3])]=1
tmp[,3][is.na(tmp[,3])]=0









for (file in files){
    name <- str_replace(file,"QTLbase_all_eQTL_","")  #replace
    name <- str_replace(name,".bed.gz","")  #replace
    name <- str_replace(name,"_narrow_peak","")  #replace
    print(name)
    if (name != "hic_loops_1" && name != "hic_loops_2" && name !="phastCons100way" ){
        print("klk")
        org2<-read.table(file,header = F,sep = "\t") %>% as.data.frame()
        org2$name <-paste(org2$V6,org2$V7,org2$V8,sep=":")
        colnames(org2)[9]<-name
        org2$key <-paste(org2$V1,org2$V2,sep=":")
        org2<-org2[,-c(1:8)] #remove top 8 column
        tmp <-left_join(tmp,org2,by = "key")
    } else if(name == "hic_loops_1" || name == "hic_loops_2" ){ #hi-C
        org2<-read.table(file,header = F,sep = "\t") %>% as.data.frame()
        org2$name <-paste(org2$V6,org2$V7,org2$V8,org2$V12,sep=":")
        colnames(org2)[13]<-name
        org2$key <-paste(org2$V1,org2$V2,sep=":")
        org2<-org2[,-c(1:12)] #remove top 12 column
        tmp <-left_join(tmp,org2,by = "key")
    }else{
        org2<-read.table(file,header = F,sep = "\t") %>% as.data.frame()
        colnames(org2)[9]<-"conservation"
        org2$key <-paste(org2$V1,org2$V2,sep=":")
        org2<-org2[,-c(1:8)]
        tmp <-left_join(tmp,org2,by = "key")
    }
}
print("finish")

tmp<-tmp[,-c(1:4)]
aa$y[!is.na(aa$y)]=1
aa$d[is.na(aa$d)]=0


gz1 <- gzfile("01_part_merge_annotation1.txt.gz", "w")
write.table(tmp,gz1,row.names = F, col.names = T,quote =F,sep="\t")
close(gz1)
# write.table(tmp,"01_part_merge_annotation.txt",row.names = F, col.names = T,quote =F,sep="\t")


result<-lm (emplambda ~ ctcf + RBP+circRNA+cistrome_CHROMATIN_Accessibility+cistrome_HISTONE+cistrome_TFBS+enhancer+xx1+xx2+lncRNA+conservation+promoter, data= tmp1)
summary(result)




all_info<-tmp


DT <- data.table(x = letters[1:5], y = 1:5)
DT[5,1]<-"NA"

DT2 <-data.table(x = letters[1:6], c = 1:6,d="NA")

aa<-left_join(DT2,DT,by= "x")
aa$y[!is.na(aa$y)]=1
aa$d[is.na(aa$d)]=0