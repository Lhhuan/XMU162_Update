library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(parallel)
library(conflicted)

types= c("cis_1MB","cis_10MB","trans_1MB","trans_10MB")
groups = c("hotspot","non_hotspot")
# for (type in types){
# ProcessBedGz2 <- function(type = NULL){
#   for (group in groups){
# type= "cis_1MB"
type= "cis_10MB"
group = "hotspot"
print(group)
print(type)
#---------------------prepare org background file

setwd(paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans",group,sep="/"))
org1<-read.table(paste("interval_18_cutoff_7.3_",type,"_eQTL_segment_",group,".bed.gz",sep=""),header = F,sep = "\t") %>% as.data.frame()
colnames(org1)[1]<-"chr"
colnames(org1)[2]<-"start"
colnames(org1)[3]<-"end"
org1$key <-paste(org1$chr,org1$start,org1$end,sep = ":")
tmp<-org1 #9904240
# tmp<-tmp[,-c(1:4)]
rownames(tmp)<-tmp$key
tmp<-tmp%>%dplyr::select(-key) 
#-------------------------------------
setwd(paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_18/annotation",group,type,sep="/"))


files<- list.files(paste("~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_18/annotation",group,type,sep="/"),"bed.gz")
ProcessBedGz <- function(file = NULL){
# for (file in file2){
  # name <- str_replace(file,paste("_interval_18_cutoff_7.3_",type,"_eQTL_segment_",group,".bed.gz"),"")  #replace
  name <- str_replace(file,"_interval_18_cutoff_7.3_cis_1MB_eQTL_segment_hotspot.bed.gz","")  #replace
  name <- str_replace(name,"_interval_18_cutoff_7.3_cis_10MB_eQTL_segment_hotspot.bed.gz","")  #replace
  name <- str_replace(name,"_interval_18_cutoff_7.3_trans_10MB_eQTL_segment_hotspot.bed.gz","")  #replace
  name <- str_replace(name,"_interval_18_cutoff_7.3_trans_1MB_eQTL_segment_hotspot.bed.gz","")  #replace
  name <- str_replace(name,"_interval_18_cutoff_7.3_trans_1MB_eQTL_segment_non_hotspot.bed.gz","")  #replace
  name <- str_replace(name,"_interval_18_cutoff_7.3_trans_10MB_eQTL_segment_non_hotspot.bed.gz","")  #replace
  name <- str_replace(name,"_interval_18_cutoff_7.3_cis_10MB_eQTL_segment_non_hotspot.bed.gz","")  #replace
  name <- str_replace(name,"_interval_18_cutoff_7.3_cis_1MB_eQTL_segment_non_hotspot.bed.gz","")  #replace
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
org2 <- mclapply(files, ProcessBedGz, mc.cores = 11)

for(i in 1:length(org2)){
  tmp<-cbind(tmp,org2[[i]][rownames(tmp),])
  j=i+3
  tmp[,j][!is.na(tmp[,j])]=1
  tmp[,j][is.na(tmp[,j])]=0
  # tmp<-tmp[,-length(tmp)]
  tmp <-tmp%>%dplyr::select(-key1)
  print(i)
}

#------------------
# gz1 <- gzfile(paste("~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/annotation_out/",group,"/",type,"/01_annotation_merge_",group,"_",type,".txt.gz",sep = ""), "w")
gz1 <- gzfile(paste("~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_18/annotation_out/05_annotation_merge_",group,"_",type,".txt.gz",sep = ""), "w")
write.table(tmp,gz1,row.names = F, col.names = T,quote =F,sep="\t")
close(gz1)
# rm(org2)
print("finish_org2")

    #----------------
#   }
# }

# mclapply(types, ProcessBedGz2, mc.cores = 5)

