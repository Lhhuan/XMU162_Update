
setwd("D:\\Huan_R\\prog\\data")
library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)
library(cowplot)
library(ggpubr)
library(Hmisc)
library(reshape2)
load("01_add_age_raw_pfs_os_filter_grade.Rdata")
org_grade <-dat
#--------------------
org_grade$birth <-as.Date(org_grade$birth)
org_grade$diagnosis <-as.Date(org_grade$diagnosis)
org_grade$age_count <-difftime(org_grade$diagnosis,org_grade$birth,units = c("days"))
org_grade$age_count <-as.numeric(org_grade$age_count)
org_grade$age_count <-org_grade$age_count/365
#-------------------
org_grade$Mono <- as.numeric(org_grade$Mono)
org_grade$Lym <-as.numeric(org_grade$Lym)
org_grade$Lym_Mono <-org_grade$Lym/org_grade$Mono
# hist(org_grade$Lym_Mono)
org_f_t <- org_grade %>%select(grade,age_count,Ki.67,LN_num,LDH,extend_num,spleen,BM,Lym_Mono,B2mg,stage,SPD,SUVmax)
org_f_t[,c(2:ncol(org_f_t))] <-apply(org_f_t[,c(2:ncol(org_f_t))],2,as.numeric)
grade_num <-group_by(org_f_t,grade)%>%summarise(count=n())%>%as.data.frame()
grade_num$class <-"grade_num"
age_more <-filter(org_f_t,age_count>60)%>%group_by(grade)%>%summarise(count=n())%>%as.data.frame()
age_more$class <-"age>60"
ki67_m <-filter(org_f_t,Ki.67>=20)%>%group_by(grade)%>%summarise(count=n())%>%as.data.frame()
ki67_m$class <-"ki67>=20"
LN_num_m <-filter(org_f_t,LN_num>=6)%>%group_by(grade)%>%summarise(count=n())%>%as.data.frame()
LN_num_m$class <-"LN_num>=6"
extend_num_m <-filter(org_f_t,extend_num>0)%>%group_by(grade)%>%summarise(count=n())%>%as.data.frame()
extend_num_m$class <-"extend_num"
stage1_m <-filter(org_f_t,stage>2)%>%group_by(grade)%>%summarise(count=n())%>%as.data.frame()
stage1_m$class <-"stage:3-4"
stage2_m <-filter(org_f_t,stage<3)%>%group_by(grade)%>%summarise(count=n())%>%as.data.frame()
stage2_m$class <-"stage:1-2"
# ki.67 30,20
# LN_num 6
# extend_num 0
suv_m <-filter(org_f_t,SUVmax>2)%>%group_by(grade)%>%summarise(count=n())%>%as.data.frame()
suv_m$class <- "SUVmax>2"
LDH_m <-filter(org_f_t,LDH>300)%>%group_by(grade)%>%summarise(count=n())%>%as.data.frame()
LDH_m$class <- "LDH>300"
spleen_m <-filter(org_f_t,spleen==1)%>%group_by(grade)%>%summarise(count=n())%>%as.data.frame()
spleen_m$class <- "spleen"
BM_m <-filter(org_f_t,BM==1)%>%group_by(grade)%>%summarise(count=n())%>%as.data.frame()
BM_m$class <- "BM"
Lym_Mono_m <-filter(org_f_t,Lym_Mono>10)%>%group_by(grade)%>%summarise(count=n())%>%as.data.frame()
Lym_Mono_m$class <- "Lym/Mono>10"
B2mg_m <-filter(org_f_t,B2mg >3.4)%>%group_by(grade)%>%summarise(count=n())%>%as.data.frame()
B2mg_m$class <-"B2mg>3.4"
spd_m <-filter(org_f_t,SPD >0)%>%group_by(grade)%>%summarise(count=n())%>%as.data.frame()
spd_m$class <-"SPD>0"
# stage
# suv
# spd
aaa <-bind_rows(grade_num,age_more,ki67_m,LN_num_m,extend_num_m,stage1_m,stage2_m,LDH_m,spleen_m,BM_m,Lym_Mono_m,B2mg_m,suv_m,spd_m)
bbb <-dcast(aaa[,c(1,3,2)],grade~class)%>%t()
colnames(bbb)<-bbb[1,]
bbb<-bbb[-1,]
ccc <-apply(bbb, 2,as.numeric)%>%as.data.frame()
rownames(ccc)<-rownames(bbb)
colnames(ccc) <-c("Num_0","Num_3a","Num_3b")
ccc$Perc_0 <- round(ccc$Num_0/ccc[5,1]*100,2)
ccc$Perc_3a <- round(ccc$Num_3a/ccc[5,2]*100,2)
ccc$Perc_3b <- round(ccc$Num_3b/ccc[5,3]*100,2)
# ddd <-ccc[c(5,1:4,6:nrow(ccc)),c(1,4,2,5,3,6)]
ddd <-ccc[c("grade_num","age>60","stage:1-2","stage:3-4","LN_num>=6","extend_num","BM","spleen","B2mg>3.4","LDH>300",
            "Lym/Mono>10","ki67>=20","SUVmax>2","SPD>0"),c(1,4,2,5,3,6)]
rownames(ddd)<-c("grade_num","age(>60)","stage(I-II)","stage(III-IV)","LN_num(>=6)","extend_num","BM","spleen","B2mg(>3.4)","LDH>300",
                 "Lym/Mono(>10)","ki67(>=20)","SUVmax(>2)","SPD(>0)")

ddd1<-ddd
ddd1$variable <-rownames(ddd1)
ddd1 <-ddd1[,c(7,1:6)]
write.table(ddd1,"02_feature_distribution_in_grade_0_3a_3b.txt",row.names = F,col.names = T,quote = F,sep = "\t")
ddd_r <-ddd
#--------3a,0
ddd_0_3a <-ddd[,1:4]

chi_3a_0 <-function(i){
  n3a1 <-ddd[i,3]
  n3a0 <-ddd[1,3]-n3a1
  n01 <-ddd[i,1]
  n00 <-ddd[1,1]-n01
  tab<-as.table(cbind(c(n3a1,n01), c(n3a0,n00)))
  dimnames(tab) <- list(c("3a", "0"),c("1", "0"))
  tab_Xsqtest <- chisq.test(tab)
  x=tab_Xsqtest$statistic
  names(x)<-NULL
  x<-round(x,2)
  X_result =data.frame(P_value=tab_Xsqtest$p.value,Chi_squared=x,class=rownames(ddd)[i])
  return(X_result)
}
chi_re_3a0 <-lapply(c(2:nrow(ddd)),chi_3a_0)
f_chi_3a0 <-do.call(rbind,chi_re_3a0)
# f_chi_3a0$pval_star <-NA

pstar <-function(i){
  p_value <-i
  # if(p_value <= 0.0001){
  #   pval_star = "****"
  if(p_value <= 0.001){
    pval_star = "***"
  }else if(p_value <= 0.01){
    pval_star = "**"
  }else if(p_value <= 0.05){
    pval_star = "*"
  }else{
    pval_star = "ns"
  }
  return(pval_star)
}
pp <-lapply(f_chi_3a0[c(1:nrow(f_chi_3a0)),1],pstar)
pp_r <-do.call(rbind,pp)
colnames(pp_r)<-"P_star"
f_chi_3a0<-bind_cols(f_chi_3a0,pp_r)
ddd_0_3a$class <-rownames(ddd_0_3a)
f_0_3a <-left_join(ddd_0_3a,f_chi_3a0,by="class")
f_0_3a <-f_0_3a[,c(5,1:4,6,8,7)]
write.table(f_0_3a,"02_feature_0_3a_chi.txt",col.names = T,row.names = F,quote=F,sep = "\t")

#----------------3a,other

eee <-ddd_r
eee$num_0_3b <- eee$Num_0 +eee$Num_3b
eee$perc_0_3b <- eee$num_0_3b/eee[1,"num_0_3b"]*100
eee <- eee[,c(7,8,3,4,1,2,5,6)]
ddd <-eee
ddd_3a_other <-ddd[,1:4]
chi_re_3a_other <-lapply(c(2:nrow(ddd)),chi_3a_0)
f_chi_3a_other <-do.call(rbind,chi_re_3a_other)
pp <-lapply(f_chi_3a_other[c(1:nrow(f_chi_3a_other)),1],pstar)
pp_r <-do.call(rbind,pp)
#-------------------------
colnames(pp_r)<-"P_star"
f_chi_3a_other<-bind_cols(f_chi_3a_other,pp_r)
ddd_3a_other$class <-rownames(ddd_3a_other)
f_3a_other <-left_join(ddd_3a_other,f_chi_3a_other,by="class")
f_3a_other <-f_3a_other[,c(5,1:4,6,8,7)]
write.table(f_3a_other,"02_feature_other_3a_chi.txt",col.names = T,row.names = F,quote=F,sep = "\t")
#---------------------------------3a,3b
# fff <-ddd_r
ddd_3a_3b <- ddd_r[,c(5,6,3,4)]
ddd<-ddd_r[,c(5,6,3,4)]
chi_re_3a_sb <-lapply(c(2:nrow(ddd)),chi_3a_0)
f_chi_3a_3b <-do.call(rbind,chi_re_3a_sb)
pp <-lapply(f_chi_3a_3b[c(1:nrow(f_chi_3a_3b)),1],pstar)
pp_r <-do.call(rbind,pp)
#-------------------------
colnames(pp_r)<-"P_star"
f_chi_3a_3b<-bind_cols(f_chi_3a_3b,pp_r)
ddd_3a_3b$class <-rownames(ddd_3a_3b)
f_3a_3b <-left_join(ddd_3a_3b,f_chi_3a_3b,by="class")
f_3a_3b <-f_3a_3b[,c(5,3,4,1,2,6,8,7)]
write.table(f_3a_3b,"02_feature_3b_3a_chi.txt",col.names = T,row.names = F,quote=F,sep = "\t")


