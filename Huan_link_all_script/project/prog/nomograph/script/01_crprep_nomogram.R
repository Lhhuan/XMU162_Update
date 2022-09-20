library(ggplot2)
library(dplyr)
library(stringr)
library(regplot)
library("survival")
rm(list=ls())


setwd("/home/huanhuan/project/prog/nomograph/output/")

all <-read.csv("/home/huanhuan/project/prog/data/all_data.2022-3-21-new.csv",na.strings = "")

# setwd("/home/huanhuan/project/prog/nomograph/script/")
load("./dat_filter.Rdata")
features = c('stage','Bsym','LN6','BM','spleen','extend_num','ECOG','B2MG_re0_train','LDH_re0_train',"pfs_month_new","pro_status","grade")
dat_used <-dat_filter[,features]

# for (i in 1:9){
#   dat_used[,i] <- as.factor(dat_used[,i])
#   print(i)
# }

colnames(dat_used)[1:10] <-c('Stage','Bsym','LoDLIN_6cm','BM','Spleen','Extra_sites','ECOG','B2MG','LDH','Time')
#-----------------3a
dat3a <-filter(dat_used,grade=="3a")
# res.cox <- coxph(Surv(pfs_month_new, pro_status == 1)~stage+Bsym+LN6+BM+spleen+extend_num+ECOG+B2MG+LDH, data = dat3a)
res.cox <- coxph(Surv(Time, pro_status == 1)~as.factor(Stage)+as.factor(Bsym)+as.factor(LoDLIN_6cm)+as.factor(BM)+as.factor(Spleen)+
                   as.factor(Extra_sites)+as.factor(ECOG)+as.factor(B2MG)+as.factor(LDH), data = dat3a)

regplot(res.cox,plots=c("density","spikes"), observation=FALSE,points=TRUE,dencol="#84B1ED",spkcol="#84B1ED",
        failtime = c(12, 24),title="3A")
#-----------------0
dat0 <-filter(dat_used,grade=="0")
res.cox <- coxph(Surv(Time, pro_status == 1)~Stage+Bsym+LoDLIN_6cm+BM+Spleen+Extra_sites+ECOG+B2MG+LDH, data = dat0)
regplot(res.cox,plots=c("density","spikes"), observation=FALSE,points=TRUE,dencol="#84B1ED",spkcol="#84B1ED",
        failtime = c(12, 24),title="Grade0-2")


