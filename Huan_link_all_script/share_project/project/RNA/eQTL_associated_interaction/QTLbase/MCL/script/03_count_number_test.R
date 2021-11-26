#------------compute-0-6 /state/partition1/huan

library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(parallel)

# setwd("/state/partition1/huan")
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/MCL/output/")
LD<-read.table("02_hotspot_in_1kg_phase3_v5_eur_1MB_100_0.2.ld.gz",header = T,sep = "") %>% as.data.frame()
# tmp<-str_split_fixed(LD$SNP_A,pattern='_',n=2)
# LD$chr <-tmp[,1]
# LD$unique_bind<-paste(LD$chr,LD$bind,sep="_")
# unique_bind<-unique(LD$unique_bind)
# length(unique_bind) #766
nrow(LD) #342079971
#---------------------count snp
SNP_A <-unique(LD$SNP_A)
SNP_A<-data.frame(snp = SNP_A)
SNP_B <-unique(LD$SNP_B)
SNP_B<-data.frame(snp = SNP_B)
all_snp <-rbind(SNP_A,SNP_B)
n<-nrow(all_snp) #10089877
print("aaa")
#---------------------------
library(MCL)
rs<-data.frame()
ProcessBedGz <- function(bind1 = NULL){
# for(j in (1:length(unique_bind))){
    LD_I <-filter(LD,unique_bind ==bind1)
    # print(j)
    LD_I_need <-LD_I%>%dplyr::select(SNP_A,SNP_B,R2)
    #-------------------------
    LD_I_need_M <-as.matrix(LD_I_need)
    name<- as.vector(LD_I_need_M[,1:2])%>%unique()
    LD_I_M  <- matrix(NA, nc = length(name), nr = length(name))
    rownames(LD_I_M) <- name; colnames(LD_I_M) <- name
    for(i in 1:nrow(LD_I_need)){
        LD_I_M[LD_I_need_M[i,1], LD_I_need_M[i,2]] <- as.numeric(LD_I_need_M[i,3])
    }
    LD_I_M[is.na(LD_I_M)] <- 0
    result <-mcl(LD_I_M, addLoops = T, expansion = 2, inflation = 1.4)
    #-----------------------
    re_mc <-as.data.frame(result$Cluster)
    re_mc$snp <- name
    tmp2<-str_split_fixed(bind1,pattern='_',n=2)
    bind2= tmp2[1,2]
    re_mc$band <-bind2
    colnames(re_mc)[1]<-"Cluster_number"
      gc() #remove space
    print("finish")
    return(re_mc)
    
    # rs<-bind_rows(rs,re_mc)
}

# unique_bind_head <-unique_bind[746:766]
# unique_bind_head <-unique_bind[1:4]
# unique_bind_head <-unique_bind[5:10]
rm(re_mc)
unique_bind_head <-unique_bind[11:15]

re_mc <- mclapply(unique_bind_head, ProcessBedGz, mc.cores = 30)
print("finish")

final_result<-data.frame()

for(i in 1:length(re_mc)){
  final_result<-bind_rows(final_result,re_mc[[i]])
  print(i)
}

# gz1 <- gzfile("MCL_result_top_4.txt.gz", "w")
gz1 <- gzfile("MCL_result_top_5_10.txt.gz", "w")
#write.csv(df1, gz1)
write.table(final_result,gz1,row.names = F, col.names = T,quote =F,sep="\t")
close(gz1)
    
    LD_I_M <-dcast(LD_I_need,SNP_A~SNP_B)
    rownames(LD_I_M)<-LD_I_M$Var1
    LD_I_M<-LD_I_M[,-1]
    LD_I_M[is.na(LD_I_M)] <- 0
    LD_I_M<-as.matrix(LD_I_M)
    result <-mcl(LD_I_M, addLoops = T, expansion = 2, inflation = 1.4)

#----------------
#-------------------test
aa<-head(LD_I_need,10)
LD_I_M_a <-dcast(aa,SNP_A~SNP_B)
rownames(LD_I_M_a)<-LD_I_M_a$SNP_A
LD_I_M_a<-LD_I_M_a[,-1]



#-------------------------------
a <- as.matrix(aa)
name <- as.vector(a[,1:2])%>%unique()

m <- data.frame(SNP_A = name,
                SNP_B = name,
                R2 = NA,
                stringsAsFactors = FALSE)
m <- matrix(NA, nc = length(name), nr = length(name))
rownames(m) <- name; colnames(m) <- name

for(i in 1:nrow(aa)){
  m[a[i,1], a[i,2]] <- as.numeric(a[i,3])
}
m[is.na(m)] <- 0
result_m <-mcl(m, addLoops = T, expansion = 2, inflation = 1.4)
re_mc <-as.data.frame(result_m$Cluster)
re_mc$snp <- name
tmp2<-str_split_fixed(unique_bind[j],pattern='_',n=2)
band= tmp2[1,2]
re_mc$band <-band
colnames(re_mc)[1]<-"Cluster_number"

rs<-data.frame()
rs<-bind_rows(rs,re_mc)



i=1 #1404453
LD_I <-filter(LD,CHR_A==i)
#--------------
SNP_A <-unique(LD_I$SNP_A)
SNP_A<-data.frame(snp = SNP_A)
SNP_B <-unique(LD_I$SNP_B)
SNP_B<-data.frame(snp = SNP_B)
all_snp <-rbind(SNP_A,SNP_B)
n<-nrow(all_snp) 

LD_M <-dcast(LD,SNP_A~SNP_B)

difference<-LD$CHR_A - LD$CHR_B
unique_difference <-unique(difference)
#0
LD<-LD[,-c(1,2,4,5)]


#-----------------------make matrix for MCL
LD_M <-dcast(LD,SNP_A~SNP_B)
rownames(LD_M)<-LD_M$Var1
LD_M<-LD_M[,-1]
LD_M[is.na(LD_M)] <- 0
#---------------------start MCL
library(MCL)
result <-mcl(mat, addLoops = T, expansion = 2, inflation = 1.4)

#-------------






i=1
aa<-filter(all_QTL,chr==i)
pos<-aa$t
r <- dist(pos,diag = TRUE,method =  "manhattan",upper =T)







qtls = c("eQTL","sQTL","eQTL","mQTL","pQTL","hQTL","caQTL")
for(qtl in qtls){
    filename<-paste("./ALL_eQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_",qtl,".txt.gz",sep="")
    qtl<-fread("./ALL_eQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_eQTL.txt.gz",header = F,sep = "\t") %>% as.data.frame()
}

ALL_QTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_QTL.txt.gz
eQTL<-fread("./ALL_eQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_eQTL.txt.gz",header = F,sep = "\t") %>% as.data.frame()
mQTL<-fread("./ALL_mQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_mQTL.txt.gz",header = F,sep = "\t") %>% as.data.frame()



#-----------------------矩阵转换
sss<-dcast(cc,Var1~Var2)
rownames(sss)<-sss$Var1
sss<-sss[,-1]
sss[is.na(sss)] <- 0