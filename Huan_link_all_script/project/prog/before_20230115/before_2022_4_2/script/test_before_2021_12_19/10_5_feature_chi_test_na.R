
setwd("/home/huanhuan/project/prog/output/")
library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)
library(cowplot)
library(ggpubr)
library(Hmisc)
library(reshape2)
load("/home/huanhuan/project/prog/data/01_add_age_raw_pfs_os_filter_grade.Rdata")
dat$new_pod_total <-dat$pod_total
dat$new_pod_total[grep("1|2",dat$new_pod_total)] <- 0
dat$new_pod_total[grep("3|4",dat$new_pod_total)] <- 1
# dat <-filter(dat,grade=="3a")



#-------------
dat$age_raw_c=NA
dat$age_raw_c[dat$age_raw<=60]=0
dat$age_raw_c[dat$age_raw>60]=1

dat$stage_c=NA
dat$stage_c[dat$stage<3]=0
dat$stage_c[dat$stage>=3]=1

dat$LN_num_c =NA
dat$LN_num_c[dat$LN_num <=4]= 0
dat$LN_num_c[dat$LN_num >4]= 1


dat$extend_num_c <-NA
dat$extend_num_c[dat$extend_num>0]=1
dat$extend_num_c[dat$extend_num==0]=0

dat$BM_c <-dat$BM
dat$spleen_c <-dat$spleen

dat$B2mg_c <- NA
dat$B2mg_c[dat$B2mg<=2.7]=0
dat$B2mg_c[dat$B2mg > 2.7]=1

dat$LDH_c =NA
dat$LDH_c[dat$LDH <=245]= 0
dat$LDH_c[dat$LDH >245]= 1 #----

dat$Lym_Mono <- dat$Lym/dat$Mono
dat$Lym_Mono_c <-NA
dat$Lym_Mono_c[dat$Lym_Mono <=10]=0
dat$Lym_Mono_c[dat$Lym_Mono >10 ]=1

dat$Ki.67_c <-NA 
dat$Ki.67_c[dat$Ki.67 <=20.00]=0
dat$Ki.67_c[dat$Ki.67 >20.00 ]=1

dat$SUVmax_c <-NA
dat$SUVmax_c[dat$SUVmax <=2]=0
dat$SUVmax_c[dat$SUVmax >2 ]=1

dat$SPD_c <-NA
dat$SPD_c[dat$SPD <=20]=0
dat$SPD_c[dat$SPD >20]=1

dat$HGB_c =NA
dat$HGB_c[dat$HGB <120]=1
dat$HGB_c[dat$HGB >=120]=0
#------------------

grade_num <-group_by(dat,grade)%>%summarise(count=n())%>%as.data.frame()
grade_num$class <-"grade_num"
aaa <-dat[,c("grade","age_raw_c","stage_c","LN_num_c","extend_num_c","BM_c","spleen_c","B2mg_c","LDH_c","Lym_Mono_c","Ki.67_c","SUVmax_c","SPD_c","HGB_c")]

count_n<-function(i=NULL){
  cr =filter(aaa,aaa[,i]>0)%>%group_by(grade)%>%summarise(count=n())%>%as.data.frame()
  cr$class=colnames(aaa)[i]
  return(cr)
}
tmp <-lapply(c(2:ncol(aaa)),count_n)
re <-do.call(rbind,tmp)
stage1_2 <- filter(aaa,stage_c==0)%>%group_by(grade)%>%summarise(count=n())%>%as.data.frame()
stage1_2$class <-"Stage(I-II)"

re2 <-bind_rows(grade_num,re,stage1_2)
re2 <-filter(re2,grade!="3b")
re2$class <-gsub("stage_c","Stage(III-IV)",re2$class)

bbb <-dcast(re2[,c(1,3,2)],grade~class)%>%t()
colnames(bbb)<-bbb[1,]
bbb<-bbb[-1,]
ccc <-apply(bbb, 2,as.numeric)%>%as.data.frame()
rownames(ccc)<-rownames(bbb)
colnames(ccc) <-c("Num_0","Num_3a")
ccc$Perc_0 <- round(ccc$Num_0/ccc[5,1]*100,2)
ccc$Perc_3a <- round(ccc$Num_3a/ccc[5,2]*100,2)
rownames(ccc)<-c("Age >60 years","B2M >2.7 mg/L","BM","Extra sites >0","grade_num","HGB <120g/L","Ki.67 >20%","LDH >245U/L","Lymph nodes >4","Lym/Mono >10","SPD >20mm","Spleen","Stage(I-II)","Stage(III-IV)","SUVmax >2")
ddd <-ccc[c("grade_num","Age >60 years","Stage(I-II)","Stage(III-IV)","Lymph nodes >4","Extra sites >0","BM","Spleen","B2M >2.7 mg/L","LDH >245U/L",
            "Lym/Mono >10","Ki.67 >20%","SUVmax >2","SPD >20mm","HGB <120g/L"),c(1,3,2,4)]
library(tibble)
ddd <-add_column(ddd, Variable = rownames(ddd), .before = 1)
write.table(ddd,"10_5_feature_distribution_in_grade_0_3a_NA.txt",row.names = F,col.names = T,quote = F,sep = "\t")
ddd_r <-ddd


chi_t <-function(i){
  n3a1 <-ddd[i,4]
  n3a0 <-ddd[1,4]-n3a1
  n01 <-ddd[i,2]
  n00 <-ddd[1,2]-n01
  tab<-as.table(cbind(c(n3a1,n01), c(n3a0,n00)))
  dimnames(tab) <- list(c("3a", "0"),c("1", "0"))
  tab_Xsqtest <- chisq.test(tab)
  x=tab_Xsqtest$statistic
  names(x)<-NULL
  X_result =data.frame(P=round(tab_Xsqtest$p.value,4),Chi=round(x,2),class=rownames(ddd)[i])
  return(X_result)
}

pstar <-function(i){
  p <-i
  # if(p <= 0.0001){
  #   pval_star = "****"
  if(p <= 0.001){
    pval_star = "***"
  }else if(p <= 0.01){
    pval_star = "**"
  }else if(p <= 0.05){
    pval_star = "*"
  }else{
    pval_star = "ns"
  }
  return(pval_star)
}

#--------3a,0
chi_re <-lapply(c(2:nrow(ddd)),chi_t)
f_chi <-do.call(rbind,chi_re)
pp <-lapply(f_chi[c(1:nrow(f_chi)),1],pstar)
pp_r<-unlist(pp)
f_chi$P_star <- pp_r
ddd$class <-rownames(ddd)
f3a <-left_join(ddd,f_chi,by="class")
f3a<-f3a[,-6]
write.table(f3a,"10_5_feature_0_3a_chi_NA.txt",col.names = T,row.names = F,quote=F,sep = "\t")
