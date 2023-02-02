setwd("D:\\prog\\data")
library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)
library(ggpubr)
data <- read.csv("D:\\prog\\data\\before_2022_4_2\\data_tie1_utf8.csv",na.strings = "")
colnames(data)[78:79] <- c("pro_time","pro_status")
load("D:\\prog\\data\\before_2022_4_2\\01_add_age_raw_pfs_os_filter_grade.Rdata")
dat1 <-left_join(dat,data[,c("No","hospital")],by="No")

dat_error1 <-filter(dat1,pro_status==1&pod48==0)
dat1_1 <- setdiff(dat1,dat_error1)

for(i in 1:nrow(dat_error1)){
  a= as.numeric(difftime(dat_error1[i,57],dat_error1[i,38],units = "days")/30/12) #pro_time, diagnosis
  dat_error1[i,"pod_total_year"]=a
  if(a >4){
    dat_error1[i,c("pod12","pod24","pod36","pod48")] <-c(0,0,0,0)
  }else if(a >3 &a <=4){
    dat_error1[i,c("pod12","pod24","pod36","pod48")] <-c(0,0,0,1)
  }else if(a >2 &a <=3){
    dat_error1[i,c("pod12","pod24","pod36","pod48")] <-c(0,0,1,1)
  }else if(a >1 &a <=2){
    dat_error1[i,c("pod12","pod24","pod36","pod48")] <-c(0,1,1,1)
  }
  else if(a >0 &a <=1){
    dat_error1[i,c("pod12","pod24","pod36","pod48")] <-c(1,1,1,1)
  }
}
dat_error1 <-filter(dat_error1,pod_total_year >0)

dat_error1$pro_status <-as.integer(dat_error1$pro_status)
dat1 <-bind_rows(dat1_1,dat_error1[,1:62])
#==============================================

unique(dat1$hospital)
dat1$hospital[is.na(dat1$hospital)]="中山大学肿瘤防治中心"
peking <-read.csv("D:\\prog\\data\\before_2022_4_2\\B2M_LDH_extract\\peking_utf8.csv",na.strings = "")
peking <-peking[,c("No","LDHg","b2MG2")]
ruijin <-read.csv("D:\\prog\\data\\before_2022_4_2\\B2M_LDH_extract\\ruijin_utf8.csv",na.strings = "")
colnames(ruijin)=ruijin[1,]
ruijin<-ruijin[-1,]
ruijin <-ruijin[,c("No","LDH","B2-MG")]
tianzhong <-read.csv("D:\\prog\\data\\before_2022_4_2\\B2M_LDH_extract\\tianzhong_utf8.csv",na.strings = "")
colnames(tianzhong)=tianzhong[1,]
tianzhong <-tianzhong[2:16,c("No","LDH","B2-MG")]
colnames(peking) <-c("No","LDH_re01","B2MG_re01")
colnames(ruijin) <-c("No","LDH_re01","B2MG_re01")
colnames(tianzhong) <-c("No","LDH_re01","B2MG_re01")

re0 <-rbind(peking,ruijin,tianzhong)
re0$LDH_re01[re0$LDH_re01=="#N/A"]=NA
re0$B2MG_re01[re0$B2MG_re01=="#N/A"]=NA
re0$LDH_re01 <-as.numeric(re0$LDH_re01)
re0$B2MG_re01 <-as.numeric(re0$B2MG_re01)
dat1<-left_join(dat1,re0,by="No")
dat1 <-dat1[,c("No","age_raw","grade","Ki.67","stage","Bsym","LN_num","BM","spleen",
               "extend_num","LN6","SUVmax","SPD","B2mg","LDH","HGB","Mono","Lym","pfs_month",
               "OS_month","followup","dead_time","dead","pro_time","pro_status","pod12","pod24","pod36","pod48",
               "diagnosis","ECOG","pfs_month_new","os_month_new","hospital","birth","pod_total","BM_extend",
               "B2mg0","LDH0","LDH_re01","B2MG_re01")]



load("01_tie2_add_pfs_os_filter.Rdata")
tie2 <-tie2%>%dplyr::select(-c(followup_status,extend_site_out_node))
tie2 <-tie2[,c("No","age_raw","grade","Ki.67","stage","Bsym","LN_num","BM","spleen",
               "extend_num","LN6","SUVmax","SPD","B2mg","LDH","HGB","Mono","Lym","pfs_month",
               "OS_month","followup","dead_time","dead","pro_time","pro_status","pod12","pod24","pod36","pod48",
               "diagnosis","ECOG","pfs_month_new","os_month_new","hospital","birth","pod_total","BM_extend","B2mg0","LDH0")]
tie2$No <-as.character(tie2$No)
tie2$SPD <-as.character(tie2$SPD)
tie2$pfs_month <-as.integer(tie2$pfs_month)
tie2$OS_month <-as.integer(tie2$OS_month)
tie2$LDH_re01 <-NA
tie2$B2MG_re01 <-NA

dat1 <-bind_rows(dat1,tie2)

dat1$LDH_upper <-NA
dat1$B2M_upper <-NA
dat1$LDH_upper[dat1$hospital=="中国医学科学院血液病医院"]=248 #0-248U/L
dat1$B2M_upper[dat1$hospital=="中国医学科学院血液病医院"]=2.5 #0.8-2.9mg/L #1.09-2.53
dat1$LDH_upper[dat1$hospital=="山东大学齐鲁医院"]=230 #120-230U/L
dat1$B2M_upper[dat1$hospital=="山东大学齐鲁医院"]= 1.8 #1.0-3.0 mg/L
dat1$LDH_upper[dat1$hospital=="厦门大学附属第一医院"]=245
dat1$B2M_upper[dat1$hospital=="厦门大学附属第一医院"]=2.7
dat1$LDH_upper[dat1$hospital=="江苏省人民医院"]=250
dat1$B2M_upper[dat1$hospital=="江苏省人民医院"]=3
dat1$LDH_upper[dat1$hospital=="浙江省肿瘤医院"]=240 #0-240U/L
dat1$B2M_upper[dat1$hospital=="浙江省肿瘤医院"]=3 #1-3mg/L
dat1$LDH_upper[dat1$hospital=="中山大学肿瘤防治中心"]= 245
dat1$B2M_upper[dat1$hospital=="中山大学肿瘤防治中心"]=  2.7#2.5
dat1$LDH_upper[dat1$hospital=="山西省肿瘤医院"]= 250#120-250 U/L
dat1$B2M_upper[dat1$hospital=="山西省肿瘤医院"]=  3 #1.3-3.0mg/L
dat1$LDH_upper[dat1$hospital=="大连医科大学附属第二医院"]= 250#135-226U/L
dat1$B2M_upper[dat1$hospital=="大连医科大学附属第二医院"]= 3 #1-3mg/L
dat1$LDH_upper[dat1$hospital=="福建省立医院"]= 245#109~245 U/L
dat1$B2M_upper[dat1$hospital=="福建省立医院"]= 2.7 #1.3~2.7 μg/ml 
dat1$LDH_upper[dat1$hospital=="广东省人民医院"]= 250
dat1$B2M_upper[dat1$hospital=="广东省人民医院"]= 2.8
dat1$LDH_upper[dat1$hospital=="苏大附一"]= 250 #120-250
dat1$B2M_upper[dat1$hospital=="苏大附一"]= 2.7
dat1$LDH_upper[dat1$hospital=="武汉协和"]= 250
dat1$B2M_upper[dat1$hospital=="武汉协和"]= 3
# aa <-unique(dat1[,c("hospital","LDH_upper")])
# table(aa$LDH_upper)
# aa <-unique(dat1[,c("hospital","B2M_upper")])
# table(aa$B2M_upper)
dat1$new_pod_total <-dat1$pod_total
dat1$new_pod_total[grep("1|2",dat1$new_pod_total)] <- 0
dat1$new_pod_total[grep("3|4",dat1$new_pod_total)] <- 1
dat1$LDH_re0 <-NA
dat1$B2MG_re0 <-NA

dat1$LDH_re0[dat1$LDH >dat1$LDH_upper]=1
dat1$LDH_re0[dat1$LDH <=dat1$LDH_upper]=0
dat1$LDH_re0[is.na(dat1$LDH_re0)]=dat1$LDH_re01[is.na(dat1$LDH_re0)]
dat1$LDH_re0_train <-dat1$LDH_re0
dat1$LDH_re0[is.na(dat1$LDH_re0)]=dat1$LDH0[is.na(dat1$LDH_re0)]
dat1$LDH_re0[is.na(dat1$LDH_re0)&dat1$new_pod_total==0]=0
dat1$LDH_re0[is.na(dat1$LDH_re0)&dat1$new_pod_total==1]=1

dat1$B2MG_re0[dat1$B2mg >dat1$B2M_upper]=1
dat1$B2MG_re0[dat1$B2mg <=dat1$B2M_upper]=0
dat1$B2MG_re0[is.na(dat1$B2MG_re0)]=dat1$B2MG_re01[is.na(dat1$B2MG_re0)]
dat1$B2MG_re0_train <-dat1$B2MG_re0
dat1$B2MG_re0[is.na(dat1$B2MG_re0)]=dat1$B2mg0[is.na(dat1$B2MG_re0)]
dat1$B2MG_re0[is.na(dat1$B2MG_re0)&dat1$new_pod_total==0]=0
dat1$B2MG_re0[is.na(dat1$B2MG_re0)&dat1$new_pod_total==1]=1
dat1$LDH_upper[dat1$hospital=="瑞金"]= 180 #91-180U/L #LDH和B2mg为NA
dat1$B2M_upper[dat1$hospital=="瑞金"]= 2.2 #0-2.2ug/ml#LDH和B2mg为NA
dat1$LDH_upper[dat1$hospital=="北京大学"]= 240 #U/L
dat1$B2M_upper[dat1$hospital=="北京大学"]= 3 #mg/L

dat1$LDH_upper[dat1$hospital=="天津市肿瘤医院"]= 247 #U/L #LDH和B2mg为NA
dat1$B2M_upper[dat1$hospital=="天津市肿瘤医院"]= 5.5 #mg/L #LDH和B2mg为NA
# dat1 <-dat1%>%select(-c(hospital,LDH_re01,B2MG_re01))

# dat1<- dat1[which(!is.na(dat1$new_pod_total)),]
save(dat1,file="04_add_age_raw_pfs_os_filter_grade_refine_ldh_b2m.Rdata")

