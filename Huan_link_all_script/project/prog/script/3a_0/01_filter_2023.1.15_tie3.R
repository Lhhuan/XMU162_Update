setwd("D:\\prog\\data")
library(ggplot2)
library(dplyr)
library(stringr)
options(stringsAsFactors = FALSE)
all <-read.csv("data_FL3189-23.1.4_UTF8.CSV",na.strings = "")

all <-all[,c("No","age_raw","grade","Ki.67","stage","Bsym","LN_num","BM","spleen",
             "extend_num","LN6","SUVmax","SPD","B2mg","LDH","HGB","Mono","Lym","pfs_month",
             "OS_month","followup","dead_time","dead","pro_time","pro_status","pod12","pod24","pod36","pod48","followup_status","diagnosis","ECOG","LDH245","BMG3","hospital")]

colnames(all)[33:34] <-c("B2mg0","LDH0")
all$No[grep("1620",all$No)] <- 1620
all$No <-as.numeric(all$No)
dat <-filter(all,No>2561)
aa <- filter(all,No<=2561)
dat$pod12[grep("NA",dat$pod12)] <- NA
dat$pod24[grep("NA",dat$pod24)] <- NA
dat$pod36[grep("NA",dat$pod36)] <- NA
dat$pod48[grep("NA",dat$pod48)] <- NA
# dat1<- dat%>%drop_na(pro_status,pod12,pod24,pod36,pod48)
dat <-dat %>%filter(!(is.na(pro_status) & is.na(pod12) & is.na(pod24) & is.na(pod36) & is.na(pod48)))##删掉全部是"哈尔滨肿瘤医院"
# dat_na <-setdiff(dat,dat1) #全部是"哈尔滨肿瘤医院"
dat$diagnosis<-as.Date(dat$diagnosis)
dat$pro_time<-as.Date(dat$pro_time)

dat_error1 <-filter(dat,pro_status==0&pod12==1)
dat_error2 <-filter(dat,pro_status==0&pod48==1)
dat_error1 <-bind_rows(dat_error1,dat_error2)%>%unique()

dat1 <- setdiff(dat,dat_error1)

dat_error1_1<- filter(dat_error1,is.na(pro_time))
dat_error1 <- setdiff(dat_error1,dat_error1_1)
for(i in 1:nrow(dat_error1_1)){
  dat_error1_1[i,c("pod12","pod24","pod36","pod48")] <-c(0,0,0,0)
}

dat_error1$pro_status=1 #pro_time exists

for(i in 1:nrow(dat_error1)){
  a= as.numeric(difftime(dat_error1[i,24],dat_error1[i,31],units = "days")/30/12)
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
dat_error1$pro_status <-as.character(dat_error1$pro_status)
dat <-bind_rows(dat1,dat_error1[,1:35],dat_error1_1[,1:35])

dat_error3 <- filter(dat,pro_status==0)
dat[which(dat$pro_status==0),c("pod12","pod24","pod36","pod48")]<-c(0,0,0,0)

# dat_error4 <- filter(dat,pro_status!=0&pro_status!=1&!(is.na(pro_status)))
dat <- filter(dat,pro_status==0|pro_status==1|is.na(pro_status)) #rm pro_status 中的日期
dat_error4 <- filter(dat,is.na(pro_status)&!is.na(pro_time))
dat1 <- setdiff(dat,dat_error4)
dat_error4$pro_status=1
for(i in 1:nrow(dat_error4)){
  a= as.numeric(difftime(dat_error4[i,24],dat_error4[i,31],units = "days")/30/12)
  dat_error4[i,"pod_total_year"]=a
  if(a >4){
    dat_error4[i,c("pod12","pod24","pod36","pod48")] <-c(0,0,0,0)
  }else if(a >3 &a <=4){
    dat_error4[i,c("pod12","pod24","pod36","pod48")] <-c(0,0,0,1)
  }else if(a >2 &a <=3){
    dat_error4[i,c("pod12","pod24","pod36","pod48")] <-c(0,0,1,1)
  }else if(a >1 &a <=2){
    dat_error4[i,c("pod12","pod24","pod36","pod48")] <-c(0,1,1,1)
  }
  else if(a >0 &a <=1){
    dat_error4[i,c("pod12","pod24","pod36","pod48")] <-c(1,1,1,1)
  }
}

dat_error4 <-filter(dat_error4,pod_total_year >0)
dat_error4$pro_status <-as.character(dat_error4$pro_status)
dat<-bind_rows(dat1,dat_error4[,1:35])

dat <- filter(dat,!(is.na(pod12)|is.na(pod24)))
dat$pro_status[grep("0",dat$pod48)] <- 0
dat$pro_status[grep("1",dat$pod48)] <- 1

dat_error5 <- filter(dat,pod12==1&pod24==0)
dat1 <-setdiff(dat,dat_error5)

dat_error5$pro_status=1
for(i in 1:nrow(dat_error5)){
  a= as.numeric(difftime(dat_error5[i,24],dat_error5[i,31],units = "days")/30/12)
  dat_error5[i,"pod_total_year"]=a
  if(a >4){
    dat_error5[i,c("pod12","pod24","pod36","pod48")] <-c(0,0,0,0)
  }else if(a >3 &a <=4){
    dat_error5[i,c("pod12","pod24","pod36","pod48")] <-c(0,0,0,1)
  }else if(a >2 &a <=3){
    dat_error5[i,c("pod12","pod24","pod36","pod48")] <-c(0,0,1,1)
  }else if(a >1 &a <=2){
    dat_error5[i,c("pod12","pod24","pod36","pod48")] <-c(0,1,1,1)
  }
  else if(a >0 &a <=1){
    dat_error5[i,c("pod12","pod24","pod36","pod48")] <-c(1,1,1,1)
  }
}
dat_error5$pro_status <-as.character(dat_error5$pro_status)
dat<-bind_rows(dat1,dat_error5[,1:35])

unique(dat[,c("pro_status","pod12","pod24","pod36","pod48")])
#==========================================================
# dat$grade_ori <-dat$grade
dat$grade[grep("低级别",dat$grade)] <- 0
dat$grade[grep("1-2级",dat$grade)] <- 0
dat$grade[grep("转化",dat$grade)] <- "3b"
dat$grade[grep("3A",dat$grade)] <- "3a"
dat$grade[grep("3$",dat$grade)] <-"3a"
dat$grade[grep("NA",dat$grade)] <- NA
dat$grade[dat$grade == 1] <- 0
dat$grade[dat$grade == 2] <- 0

dat_3b <-filter(dat,grade=="3b")
dat <-setdiff(dat,dat_3b)

#---------------------------------------
dat[,26] <-as.numeric(dat[,26])
dat[,27] <-as.numeric(dat[,27])
dat[,28] <-as.numeric(dat[,28])
dat[,29] <-as.numeric(dat[,29])

dat$pod_total <-rowSums(dat[,26:29]) #c("pod12","pod24","pod36","pod48")

dat$Ki.67[grep("NA",dat$Ki.67)] <- NA
dat$Ki.67[grep("20-30",dat$Ki.67)] <- 25
dat$Ki.67[grep("30-60",dat$Ki.67)] <- 45
dat$Ki.67[grep("30-40",dat$Ki.67)] <- 35
dat$Ki.67[grep("40-50",dat$Ki.67)] <- 45
dat$Ki.67[grep("60-70",dat$Ki.67)] <- 65
dat$Ki.67[grep("40-60",dat$Ki.67)] <- 50
dat$Ki.67[grep("15-40",dat$Ki.67)] <- 27.5
dat$Ki.67[grep("40-90",dat$Ki.67)] <- 45
dat$Ki.67[grep("25-50",dat$Ki.67)] <- 37.5
dat$Ki.67[grep("10月20日",dat$Ki.67)] <- 15

dat$Ki.67 <-as.numeric(dat$Ki.67)
dat$stage[grep("NA",dat$stage)] <- NA
dat$stage <-as.numeric(dat$stage)
dat$Bsym[grep("NA",dat$Bsym)] <- NA
dat$Bsym <-as.numeric(dat$Bsym)

dat$LN_num[grep("NA",dat$LN_num)] <- NA
dat$LN_num[grep("多发",dat$LN_num)] <- 25 #-----多发
dat$LN_num[grep("肠系膜淋巴结",dat$LN_num)] <- 25
dat$LN_num[grep("单一",dat$LN_num)] <- 25
dat$LN_num <-as.numeric(dat$LN_num)

dat$BM[grep("NA",dat$BM)] <- NA
dat$BM <-as.numeric(dat$BM)

dat$spleen[grep("NA",dat$spleen)] <- NA
dat$spleen <-as.numeric(dat$spleen)

dat$extend_num[grep("NA",dat$extend_num)] <- NA
dat$extend_num <-as.numeric(dat$extend_num)


dat$LN6[grep("NA",dat$LN6)] <- NA
dat$LN6 <-as.numeric(dat$LN6)

dat$SUVmax[grep("NA|na",dat$SUVmax)] <- NA
dat$SUVmax <-as.numeric(dat$SUVmax)

dat$SPD[grep("NA|q",dat$SPD)] <- NA
dat$SPD <-as.numeric(dat$SPD)

dat$B2mg[grep("NA",dat$B2mg)]<- NA
dat$B2mg <-as.numeric(dat$B2mg)

dat$LDH[grep("NA",dat$LDH)]<- NA
dat$LDH <-as.numeric(dat$LDH)
dat$HGB[grep("NA",dat$HGB)]<- NA
dat$HGB <-as.numeric(dat$HGB)

dat$Mono[grep("NA",dat$Mono)]<- NA
dat$Mono <-as.numeric(dat$Mono)

dat$Lym[grep("NA",dat$Lym)]<- NA
dat$Lym <-as.numeric(dat$Lym)

dat$followup[grep("NA",ignore.case = TRUE, dat$followup)] <- NA
dat$followup<-as.Date(dat$followup)
dat$dead_time[grep("NA",ignore.case = TRUE, dat$dead_time)] <- NA
dat$dead_time[grep("17-Mar",dat$dead_time)]<- NA

dat$dead_time<-as.Date(dat$dead_time)

dat$dead[grep("NA",dat$dead)]<- 0 #dat$dead_time[grep("NA",dat$dead)] NA
dat$dead[is.na(dat$dead)] <-0 
dat$dead[grep("/",dat$dead)]<- 0 
dat$dead<-as.numeric(dat$dead)

dat$pro_status<-as.numeric(dat$pro_status)

dat$pfs_month[grep("NA",dat$pfs_month)] <- NA
dat$pfs_month<-as.numeric(dat$pfs_month)
dat$OS_month[grep("NA",dat$OS_month)] <- NA
dat$OS_month<-as.numeric(dat$OS_month)


dat$ECOG[grep("NA",dat$ECOG)] <- NA
dat$ECOG <-as.integer(dat$ECOG)
dat$LDH0[grep("NA",dat$LDH0)] <- NA
dat$LDH0 <-as.numeric(dat$LDH0)
dat$B2mg0[grep("NA",dat$B2mg0)] <- NA
dat$B2mg0 <-as.numeric(dat$B2mg0)

# dat$diagnosis<-as.Date(dat$diagnosis)
save(dat,file = "tie3_filter.Rdata")
write.csv(dat,"tie3_rawdata_filter20230115.csv",row.names = F)

pfs <- function(x){
  if(x[1,25] == 0){
    a = difftime(x[1,21],x[1,31],units = "days")}else{
      a = difftime(x[1,24],x[1,31],units = "days")}
  return(a)}

# pfs <- function(x){
#   if(x[1,"pro_status"] == 0){
#     a = difftime(x[1,'followup'],x[1,"disgnosis"],units = "days")}else{
#       a = difftime(x[1,"pro_time"],x[1,"disgnosis"],units = "days")}
#   return(a)}

re <- lapply(1:nrow(dat), function(i){
  print(i)
  x = dat[i,]
  a <- round(pfs(x)/30)
  return(a)
}) 
dat$pfs_month_new <- unlist(re)

os <- function(x){
  if(x[1,23] == 0){
    a = difftime(x[1,21],x[1,31],units = "days")}else{
      a = difftime(x[1,22],x[1,31],units = "days")}
  return(a)}

re1 <- lapply(1:nrow(dat), function(i){
  print(i)
  x = dat[i,]
  a <- round(os(x)/30)
  return(a)
}) 
dat$os_month_new <- unlist(re1)


tie3<-filter(dat,hospital !="哈尔滨肿瘤医院")
# tie3 <-dat
save(tie3,file="01_tie3_add_pfs_os_filter_20230115.Rdata")