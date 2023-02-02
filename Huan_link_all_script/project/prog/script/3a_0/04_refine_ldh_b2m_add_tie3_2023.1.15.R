setwd("D:\\prog\\data")
library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)
library(ggpubr)
load("04_add_age_raw_pfs_os_filter_grade_refine_ldh_b2m.Rdata")
d1 <-dat1[,c("No","pod12","pod24","pod36","pod48","pod_total")]
d1_f <- d1[complete.cases(d1),]
d1_na <-setdiff(d1,d1_f)
d1_na <-left_join(d1_na,dat1[,c('No','pro_status')],by="No")
unique(d1_na$pro_status)
#=================================由此判断dat1 




load("01_tie3_add_pfs_os_filter_20230115.Rdata")
tie3 <-tie3%>%dplyr::select(-c(followup_status))
tie3$BM_extend <-NA
tie3$birth <-NA

tie3 <-tie3[,c("No","age_raw","grade","Ki.67","stage","Bsym","LN_num","BM","spleen",
               "extend_num","LN6","SUVmax","SPD","B2mg","LDH","HGB","Mono","Lym","pfs_month",
               "OS_month","followup","dead_time","dead","pro_time","pro_status","pod12","pod24","pod36","pod48",
               "diagnosis","ECOG","pfs_month_new","os_month_new","hospital","birth","pod_total","BM_extend","B2mg0","LDH0")]


tie3$No <-as.character(tie3$No)
tie3$SPD <-as.character(tie3$SPD)
tie3$pfs_month <-as.integer(tie3$pfs_month)
tie3$OS_month <-as.integer(tie3$OS_month)
tie3$LDH_re01 <-NA
tie3$B2MG_re01 <-NA

tie3$LDH_upper <-NA
tie3$B2M_upper <-NA
tie3$LDH_upper[tie3$hospital=="北协和"]=248 #0-248U/L
tie3$B2M_upper[tie3$hospital=="北协和"]=2.5 #0.8-2.9mg/L #1.09-2.53
tie3$LDH_upper[tie3$hospital=="北医三院"]= 240 #U/L
tie3$B2M_upper[tie3$hospital=="北医三院"]= 3 #mg/L
tie3$LDH_upper[tie3$hospital=="哈尔滨肿瘤医院"]= 245 #U/L
tie3$B2M_upper[tie3$hospital=="哈尔滨肿瘤医院"]= 3 #mg/L

tie3$new_pod_total <-tie3$pod_total
tie3$new_pod_total[grep("1|2",tie3$new_pod_total)] <- 0
tie3$new_pod_total[grep("3|4",tie3$new_pod_total)] <- 1
tie3$LDH_re0 <-NA
tie3$B2MG_re0 <-NA

tie3$LDH_re0[tie3$LDH >tie3$LDH_upper]=1
tie3$LDH_re0[tie3$LDH <=tie3$LDH_upper]=0

tie3$LDH_re0_train <-tie3$LDH_re0
tie3$LDH_re0[is.na(tie3$LDH_re0)]=tie3$LDH0[is.na(tie3$LDH_re0)]
tie3$LDH_re0[is.na(tie3$LDH_re0)&tie3$new_pod_total==0]=0
tie3$LDH_re0[is.na(tie3$LDH_re0)&tie3$new_pod_total==1]=1

tie3$B2MG_re0[tie3$B2mg >tie3$B2M_upper]=1
tie3$B2MG_re0[tie3$B2mg <=tie3$B2M_upper]=0
tie3$B2MG_re0[is.na(tie3$B2MG_re0)]=tie3$B2MG_re01[is.na(tie3$B2MG_re0)]
tie3$B2MG_re0_train <-tie3$B2MG_re0
tie3$B2MG_re0[is.na(tie3$B2MG_re0)]=tie3$B2mg0[is.na(tie3$B2MG_re0)]
tie3$B2MG_re0[is.na(tie3$B2MG_re0)&tie3$new_pod_total==0]=0
tie3$B2MG_re0[is.na(tie3$B2MG_re0)&tie3$new_pod_total==1]=1


dat2 <-bind_rows(dat1,tie3)
save(dat2,file="04_add_tie3_refine_ldh_b2m_20230115.Rdata")



aaa <-tie3[,c("pod12","pod24","pod36","pod48","pod_total")]
aaa <-filter(aaa,pod_total==0)
aaa_f <- aaa[complete.cases(aaa),]
aaa_na <-setdiff(aaa,aaa_f)



