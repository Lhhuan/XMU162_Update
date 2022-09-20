library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)
library(cowplot)
library(ggpubr)
library(Hmisc)
library(reshape2)
setwd("/home/huanhuan/project/prog/nomograph/output/")
load("/home/huanhuan/project/prog/data/04_add_age_raw_pfs_os_filter_grade_refine_ldh_b2m.Rdata")
dat <-filter(dat1,grade!="3b")
dat$Lym_Mono <- dat$Lym/dat$Mono
library(regplot)

pro_time <- function(x){
  if(x[1,25] == 0){
    a = 0}else{
      a = difftime(x[1,24],x[1,30],units = "days")}
  return(a)}

re <- lapply(1:nrow(dat), function(i){
  print(i)
  x = dat[i,]
  a <- pro_time(x)
  return(a)
}) 
dat$pro_time_day <- unlist(re)


aa <-dat[is.na(dat$pro_time_day),]



dat_used <-dat[,c('B2MG_re0_train','Lym_Mono','LDH_re0_train','HGB','LN_num','SUVmax','BM','grade','pro_time_day','pro_status')]
colnames(dat_used) <-c("B2MG","Lym_Mono","LDH","HGB","Lymph_nodes",'SUVmax','BM','grade','time','pro_status')
#-----------------3a
dat3a <-filter(dat_used,grade=="3a")
res.cox <- coxph(Surv(time, pro_status)~B2MG+Lym_Mono+LDH+HGB+Lymph_nodes+SUVmax+BM, data = dat3a)
pdf("01_3a_nomogram.pdf")
regplot(res.cox)
dev.off()

#-----------------0

dat0 <-filter(dat_used,grade=="0")
res.cox <- coxph(Surv(time, pro_status)~B2MG+Lym_Mono+LDH+HGB+Lymph_nodes+SUVmax+BM, data = dat0)
regplot(res.cox)
dev.off()
# pdf("01_grade012_nomogram.pdf")
# dev.off()