setwd("D:\\Huan_R\\prog\\data")
library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)
library(cowplot)
library(ggpubr)
library(Hmisc)
data <- read.csv("data.csv",na.strings = "")
colnames(data)[78:79] <- c("pro_time","pro_status")
load("01_add_age_raw_pfs_os_filter_grade.Rdata")
dat1 <-left_join(dat,data[,c("No","hospital")],by="No")
unique(dat1$hospital)
dat1$hospital[is.na(dat1$hospital)]="中山大学肿瘤防治中心"
# dat1$hospital_E <-NA
# dat1$hospital_E[dat1$hospital=="福建省立医院"]="Fujian Provincial Hospital"
# dat1$hospital_E[dat1$hospital=="厦门大学附属第一医院"]="The first affiliated hospital of xiamen university"
# dat1$hospital_E[dat1$hospital=="中山大学肿瘤防治中心"]="Sun Yat-sen University Cancer Center"
# dat1$hospital_E[dat1$hospital=="广东省人民医院"]="Guangdong General Hospital"
# dat1$hospital_E[dat1$hospital=="北京大学"]="Peking University"
# dat1$hospital_E[dat1$hospital=="浙江省肿瘤医院"]="Zhejiang Cancer Hospital"
# dat1$hospital_E[dat1$hospital=="瑞金"]="Ruijin"
# dat1$hospital_E[dat1$hospital=="天津市肿瘤医院"]="Tianjin Medical University cancer hospital"
# dat1$hospital_E[dat1$hospital=="大连医科大学附属第二医院"]="The second hospital of Dalian medical university"
# dat1$hospital_E[dat1$hospital=="山西省肿瘤医院"]="Shanxi provincial cancer hospital"
# dat1$hospital_E[dat1$hospital=="江苏省人民医院"]="Jiangsu provincial hospital"
# dat1$hospital_E[dat1$hospital=="山东大学齐鲁医院"]="Qilu hospital of shandong university"
# dat1$hospital_E[dat1$hospital=="中国医学科学院血液病医院"]="Institute of Hematology & Blood Disease Hospital,Chinese Academy of Medical Sciences"
peking <-read.csv("D:\\Huan_R\\prog\\data\\B2M_LDH_extract\\peking.csv",na.strings = "")
peking <-peking[,c("No","LDHg","b2MG2")]
ruijin <-read.csv("D:\\Huan_R\\prog\\data\\B2M_LDH_extract\\ruijin.csv",na.strings = "")
colnames(ruijin)=ruijin[1,]
ruijin<-ruijin[-1,]
ruijin <-ruijin[,c("No","LDH","B2-MG")]
tianzhong <-read.csv("D:\\Huan_R\\prog\\data\\B2M_LDH_extract\\tianzhong.csv",na.strings = "")
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
dat1$LDH_upper[dat1$hospital=="山西省肿瘤医院"]= 240#114-240 U/L
dat1$B2M_upper[dat1$hospital=="山西省肿瘤医院"]=  2 #mg/L
dat1$LDH_upper[dat1$hospital=="大连医科大学附属第二医院"]= 250#135-226U/L
dat1$B2M_upper[dat1$hospital=="大连医科大学附属第二医院"]= 3 #1-3mg/L
dat1$LDH_upper[dat1$hospital=="福建省立医院"]= 245#109~245 U/L
dat1$B2M_upper[dat1$hospital=="福建省立医院"]= 2.7 #1.3~2.7 μg/ml 
dat1$LDH_upper[dat1$hospital=="广东省人民医院"]= 250
dat1$B2M_upper[dat1$hospital=="广东省人民医院"]= 2.8

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
save(dat1,file="04_add_age_raw_pfs_os_filter_grade_refine_ldh_b2m.Rdata")


# dat1$LDH_upper[dat1$hospital=="瑞金"]= 180 #91-180U/L #LDH和B2mg为NA
# dat1$B2M_upper[dat1$hospital=="瑞金"]= 2.2 #0-2.2ug/ml#LDH和B2mg为NA
# dat1$LDH_upper[dat1$hospital=="北京大学"]= 240 #U/L
# dat1$B2M_upper[dat1$hospital=="北京大学"]= 3 #mg/L

# dat1$LDH_upper[dat1$hospital=="天津市肿瘤医院"]= 247 #U/L #LDH和B2mg为NA
# dat1$B2M_upper[dat1$hospital=="天津市肿瘤医院"]= 5.5 #mg/L #LDH和B2mg为NA