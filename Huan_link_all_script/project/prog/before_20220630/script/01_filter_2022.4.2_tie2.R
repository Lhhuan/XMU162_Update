setwd("D:\\Huan_R\\prog\\data")
library(ggplot2)
library(dplyr)
library(stringr)
options(stringsAsFactors = FALSE)
all <-read.csv("all_data.2022-3-21-new.csv",na.strings = "")
all <-all[,c("No","age_raw","grade","Ki.67","stage","Bsym","LN_num","BM","spleen",
             "extend_num","LN6","SUVmax","SPD","B2mg","LDH","HGB","Mono","Lym","pfs_month",
             "OS_month","followup","dead_time","dead","pro_time","pro_status","pod12","pod24","pod36","pod48","followup_status","diagnosis","ECOG","LDH245","BMG3","hospital","birth","extend_site_out_node")]

colnames(all)[33:34] <-c("B2mg0","LDH0")
all$No[grep("1620",all$No)] <- 1620
all$No <-as.numeric(all$No)
dat <-filter(all,No>2370)
aa <- filter(all,No<=2370)

# dat$grade_ori <-dat$grade
dat$grade[grep("低级别",dat$grade)] <- 0
dat$grade[grep("转化",dat$grade)] <- "3b"
dat$grade[grep("3B",dat$grade)] <- "3b"
dat$grade[grep("高级别",dat$grade)] <- "3b"
dat$grade[grep("滤泡弥漫混合型",dat$grade)] <- "3b"
dat$grade[grep("滤泡弥漫型",dat$grade)] <- "3b"
dat$grade[grep("滤泡中心母细胞型",dat$grade)] <- "3b"
dat$grade[grep("3A",dat$grade)] <- "3a"
dat$grade[grep("3$",dat$grade)] <-"3a"
dat$grade[grep("NA",dat$grade)] <- NA
dat$grade[dat$grade == 1] <- 0
dat$grade[dat$grade == 2] <- 0
#---------------------------------------
# h <-filter(dat,grade=="3a")
# data <- read.csv("filter3_3a_exclude_bigB.2022-3-21-new.csv",na.strings = "")
# data <-filter(data,No>2370)
# diffno <-setdiff(h$No,data$No)
# add <-filter(h,No%in%diffno)
# #--------------------------
dat$pod12[grep("NA",dat$pod12)] <- NA
dat$pod24[grep("NA",dat$pod24)] <- NA
dat$pod36[grep("NA",dat$pod36)] <- NA
dat$pod48[grep("NA",dat$pod48)] <- NA
dat[,26:29] <-lapply(26:29,function(i){as.numeric(dat[,i])})

dat$pod_total <-rowSums(dat[,26:29]) #c("pod12","pod24","pod36","pod48")
dat$Ki.67[grep("NA",dat$Ki.67)] <- NA
dat$Ki.67[grep("50-60",dat$Ki.67)] <- 55
dat$Ki.67[grep("10-40",dat$Ki.67)] <- 30
dat$Ki.67[grep("10月20日",dat$Ki.67)] <- 15
dat$Ki.67[grep("10月30日",dat$Ki.67)] <- 20
dat$Ki.67[grep("20-30",dat$Ki.67)] <- 25
dat$Ki.67[grep("20-40",dat$Ki.67)] <- 30
dat$Ki.67[grep("20-60",dat$Ki.67)] <- 40
dat$Ki.67[grep("30-40",dat$Ki.67)] <- 35
dat$Ki.67[grep("40-60",dat$Ki.67)] <- 50
dat$Ki.67[grep("50-60",dat$Ki.67)] <- 55
dat$Ki.67[grep("50-70",dat$Ki.67)] <- 60
dat$Ki.67[grep("5月10日",dat$Ki.67)] <- 7.5
dat$Ki.67[grep("7-50",dat$Ki.67)] <- 28.5
dat$Ki.67 <-as.numeric(dat$Ki.67)
dat$stage[grep("NA",dat$stage)] <- NA
dat$stage <-as.numeric(dat$stage)
dat$Bsym[grep("NA",dat$Bsym)] <- NA
dat$Bsym <-as.numeric(dat$Bsym)

dat$LN_num[grep("NA",dat$LN_num)] <- NA
dat$LN_num[grep("多发",dat$LN_num)] <- 25 #-----多发
dat$LN_num[grep("肠系膜 腹腔及腹膜后多发肿大淋巴结，FDG代谢异常增高",dat$LN_num)] <- 25
dat$LN_num <-as.numeric(dat$LN_num)

dat$BM[grep("NA",dat$BM)] <- NA
dat$BM <-as.numeric(dat$BM)

dat$spleen[grep("NA",dat$spleen)] <- NA
dat$spleen <-as.numeric(dat$spleen)

dat$extend_num[grep("NA",dat$extend_num)] <- NA
dat$extend_num <-as.numeric(dat$extend_num)


dat$LN6[grep("NA",dat$LN6)] <- NA
dat$LN6 <-as.numeric(dat$LN6)

dat$SUVmax[grep("NA",dat$SUVmax)] <- NA
dat$SUVmax <-as.numeric(dat$SUVmax)

dat$SPD[grep("NA",dat$SPD)] <- NA
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
dat$Lym[grep("3,48",dat$Lym)]<-3.48
dat$Lym <-as.numeric(dat$Lym)

dat$followup[grep("NA",ignore.case = TRUE, dat$followup)] <- NA
dat$followup<-as.Date(dat$followup)
dat$dead_time[grep("NA",ignore.case = TRUE, dat$dead_time)] <- NA
dat$dead_time<-as.Date(dat$dead_time)

dat$dead[grep("NA",dat$dead)]<- 0 #dat$dead_time[grep("NA",dat$dead)] NA
dat$dead<-as.numeric(dat$dead)

dat$pro_time<-as.Date(dat$pro_time)
dat$pro_status[grep("NA",dat$pro_status)]<- 0 #dat$pro_time[grep("NA",dat$pro_status)]
dat$pro_status<-as.numeric(dat$pro_status)


dat$pfs_month[grep("NA",dat$pfs_month)] <- NA
dat$OS_month[grep("NA",dat$OS_month)] <- NA
dat$OS_month[grep("/",dat$OS_month)] <- NA
dat$BM_extend <-dat$extend_site_out_node
dat$BM_extend[grep("NA",dat$BM_extend)] <- NA
dat$BM_extend[grep("1|胃|脾|肺|肝|舌|胰|皮肤|肌肉|口|咽|鼻|肠|肋骨|髂骨|腮腺|腹腔|肾|其他|左|/|扁桃体|甲状腺|,|呼吸",dat$extend_site_out_node)]<- 1
dat$BM_extend[grep("无",dat$extend_site_out_node)]<- 0
dat$BM_extend[grep("骨髓",dat$BM_extend)]<- 0
dat$BM_extend[grep("骨",dat$BM_extend)]<- 1
dat$BM_extend<-as.numeric(dat$BM_extend)
dat$ECOG[grep("NA",dat$ECOG)] <- NA
dat$ECOG <-as.integer(dat$ECOG)
dat$LDH0[grep("NA",dat$LDH0)] <- NA
dat$LDH0 <-as.numeric(dat$LDH0)
dat$B2mg0[grep("NA",dat$B2mg0)] <- NA
dat$B2mg0 <-as.numeric(dat$B2mg0)

dat$diagnosis<-as.Date(dat$diagnosis)
save(dat,file = "tie2_filter.Rdata")
write.csv(dat,"tie2_rawdata_filter.csv",row.names = F)

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


# dat$new_pod_total <-dat$pod_total
# dat$new_pod_total[grep("1|2",dat$new_pod_total)] <- 0
# dat$new_pod_total[grep("3|4",dat$new_pod_total)] <- 1
# dat$
dat$pro_status[is.na(dat$pod_total)]
dat$pod_total[is.na(dat$pod_total)]=0
tie2<-dat

save(tie2,file="01_tie2_add_pfs_os_filter.Rdata")


