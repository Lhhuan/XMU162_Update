library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)


#-----------------------------------------new
setwd("~/project/RNA/eQTL_associated_interaction/QTLbase/MCL/output/")
qtl1s = c("eQTL","mQTL")
rs<-data.frame()
for(qtl in qtl1s){
        #--------------top
        file_name<-paste("../../output/ALL_",qtl,"/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_",qtl,".txt.gz",sep="")
        org<-read.table(file_name,header=T,sep="\t")%>%as.data.frame()%>%arrange(desc(emplambda))
        top<-filter(org,emplambda>0.05)
        org_pos <-top%>%dplyr::select(t,chr)
        final_org_pos<-paste(org_pos$chr,org_pos$t,sep="_")%>%as.data.frame()
        colnames(final_org_pos)[1]<-"SNP"
        final_org_pos$SNP <-as.character(final_org_pos$SNP)
        #---------------------------------------------all xQTL
        all_qtl_file_name <-paste("/share/data0/QTLbase/data/",qtl,".txt.gz",sep="")
        # all_qtl <-fread(all_qtl_file_name,header=T,sep="\t")%>%as.data.frame()
        all_qtl<-read.table(all_qtl_file_name,header = T,sep = "\t") %>% as.data.frame()
        all_qtl_pos <-all_qtl%>%dplyr::select(SNP_chr,SNP_pos,Pvalue)
        all_qtl_pos$SNP <- paste(all_qtl_pos$SNP_chr,all_qtl_pos$SNP_pos,sep="_")
        all_qtl_pos<-all_qtl_pos[,-c(1,2)]
        colnames(all_qtl_pos)[1] <-"P"
        #--------------------------- find the min p of xQTL
        final_top<-inner_join(final_org_pos,all_qtl_pos,by="SNP")
        final_top1<-final_top%>%group_by(SNP)%>%summarise(min(P))%>%as.data.frame()
        rs<-bind_rows(rs,final_top1)
}
#--------------top
qtl<-"miQTL"
file_name<-paste("../../output/ALL_",qtl,"/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_",qtl,".txt.gz",sep="")
org<-read.table(file_name,header=T,sep="\t")%>%as.data.frame()%>%arrange(desc(emplambda))
top<-filter(org,emplambda>0.02)
org_pos <-top%>%dplyr::select(t,chr)
final_org_pos<-paste(org_pos$chr,org_pos$t,sep="_")%>%as.data.frame()
colnames(final_org_pos)[1]<-"SNP"
final_org_pos$SNP <-as.character(final_org_pos$SNP)
#---------------------------------------------all
all_qtl_file_name <-paste("/share/data0/QTLbase/data/",qtl,".txt.gz",sep="")
# all_qtl <-fread(all_qtl_file_name,header=T,sep="\t")%>%as.data.frame()
all_qtl<-read.table(all_qtl_file_name,header = T,sep = "\t") %>% as.data.frame()
all_qtl_pos <-all_qtl%>%dplyr::select(SNP_chr,SNP_pos,Pvalue)
all_qtl_pos$SNP <- paste(all_qtl_pos$SNP_chr,all_qtl_pos$SNP_pos,sep="_")
all_qtl_pos<-all_qtl_pos[,-c(1,2)]
colnames(all_qtl_pos)[1] <-"P"
#---------------------------inner_join
final_top<-inner_join(final_org_pos,all_qtl_pos,by="SNP")
final_top1<-final_top%>%group_by(SNP)%>%summarise(min(P))%>%as.data.frame()
final_rs<-bind_rows(rs,final_top1)
#------------------------------------all_inner_join,find the min p of xQTL
colnames(final_rs)[2] <-"P"
final_rs1<-final_rs%>%group_by(SNP)%>%summarise(min(P))%>%as.data.frame()
#----------------------------------------------------------
# final_rs$pos <-paste(final_rs$chr,final_rs$t,sep="_")
# unique_pos<-unique(final_rs$pos)

# tmp <- str_split_fixed(unique_pos,pattern='_',n=2)
# output<-data.frame(chr= tmp[,1],pos=tmp[,2])
# # bb<-as.data.frame(unique_pos)
# output<-data.frame(SNP=unique_pos)
write.table(final_rs1,"01_top_QTL_eQTL_mQTL_miQTL.txt",row.names = F, col.names = T,quote =F,sep="\t")#把表存下来




#---------------------------------------

# qtls = c("eQTL","riboQTL","reQTL","miQTL","metaQTL","lncRNAQTL","edQTL","cerQTL","pQTL","caQTL","sQTL","hQTL","mQTL")

# setwd("~/project/RNA/eQTL_associated_interaction/QTLbase/MCL/output/")
# qtl1s = c("eQTL","mQTL")
# rs<-data.frame()
# for(qtl in qtl1s){
#         file_name<-paste("../../output/ALL_",qtl,"/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_",qtl,".txt.gz",sep="")
#         org<-fread(file_name,header=T,sep="\t")%>%as.data.frame()%>%arrange(desc(emplambda))
#         top<-filter(org,emplambda>0.05)
#         rs<-bind_rows(rs,top)
# }
# qtl<-"miQTL"
# file_name<-paste("../../output/ALL_",qtl,"/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_",qtl,".txt.gz",sep="")
# org<-fread(file_name,header=T,sep="\t")%>%as.data.frame()%>%arrange(desc(emplambda))
# top<-filter(org,emplambda>0.02)
# final_rs<-bind_rows(rs,top)
# final_rs$pos <-paste(final_rs$chr,final_rs$t,sep="_")
# unique_pos<-unique(final_rs$pos)

# tmp <- str_split_fixed(unique_pos,pattern='_',n=2)
# output<-data.frame(chr= tmp[,1],pos=tmp[,2])
# # bb<-as.data.frame(unique_pos)
# output<-data.frame(SNP=unique_pos)
# write.table(output,"01_top_QTL_eQTL_mQTL_miQTL.txt",row.names = F, col.names = T,quote =F,sep="\t")#把表存下来