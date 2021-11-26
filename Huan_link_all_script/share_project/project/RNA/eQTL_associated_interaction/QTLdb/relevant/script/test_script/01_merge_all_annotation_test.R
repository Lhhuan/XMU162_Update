library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)


setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/")
org1<-read.table("QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL_normalized_sorted.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(org1)[1]<-"chr"
colnames(org1)[2]<-"start"
colnames(org1)[3]<-"end"
colnames(org1)[4]<-"pos"
colnames(org1)[5]<-"emplambda"
org1$key <-paste(org1$chr,org1$start,sep = ":")
tmp<-org1 #9904240
tmp<-tmp[,-c(1:4)]
rm(org1)
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/relevant/output/annotation")

files<- list.files("~/project/RNA/eQTL_associated_interaction/QTLbase/relevant/output/annotation","bed.gz")
# files<-c(files[1],files[2])
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