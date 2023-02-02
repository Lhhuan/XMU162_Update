setwd("/mnt/c/Users/chenq/Desktop/FL/")
library(ggplot2)
library(dplyr)
library(stringr)
options(stringsAsFactors = FALSE)
data <- read.csv("./data/data.csv",na.strings = "")
data <- data[1:2370,]
dat <- data[,-c(2,8,9,10,12,14,17,18,20,33,35,51,52,56,60,69,78,79)]
dat$pod_total <- rowSums(dat[,2:5])
dat$grade[grep("低级别",dat$grade)] <- 0
dat$grade[grep("转化",dat$grade)] <- "3b"
dat$grade[grep("3A",dat$grade)] <- "3a"
dat$grade[grep(" ",dat$grade)] <- NA
dat$grade[grep("儿童",dat$grade)] <- NA
dat$grade[grep("NA",dat$grade)] <- NA
dat$grade[dat$grade == 1] <- 0
dat$grade[dat$grade == 2] <- 0

dat$Ki.67[grep("NA",dat$Ki.67)] <- NA
dat$stage[grep("NA",dat$stage)] <- NA
dat$Bsym[grep("NA",dat$Bsym)] <- NA
dat$LN_num[grep("NA",dat$LN_num)] <- NA
dat$site0[grep("NA",dat$site0)] <- NA
dat$extend[grep("侵及脾 侵及骨髓",dat$extend)] <- 2
dat$extend[grep("不详",dat$extend)] <- NA
dat$extend[grep(" ",dat$extend)] <- NA
dat$extend[grep("是|有",dat$extend)] <- NA
dat$extend[grep("NA",dat$extend)] <- NA

dat$BM[grep("NA",dat$BM)] <- NA
dat$spleen[grep("NA",dat$spleen)] <- NA
dat <- select(dat,!extend_site)
dat$extend_num[grep("NA",dat$extend_num)] <- NA
dat$BM_extend[grep("NA",dat$BM_extend)] <- NA
dat$LN6[grep("NA",dat$LN6)] <- NA
dat$SUVmax[grep("NA|无|/",dat$SUVmax)] <- NA
dat$SUVmax[grep("已切除",dat$SUVmax)] <- 1.2
dat$SUVmax[grep("1.2-1.3",dat$SUVmax)] <- 1.25
dat$SUVmax[grep("-",dat$SUVmax)] <- NA
dat$SPD[grep("/|NA|外院",dat$SPD)] <- NA
dat$SPD[grep(" ",dat$SPD)]<- NA
dat$primapi[grep("NA",dat$primapi)]<- NA
dat$X150b2mg_ldh[grep("NA",dat$X150b2mg_ldh)]<- NA
dat$b2mg_LDH[grep("NA",dat$b2mg_LDH)]<- NA
dat$ECOG[grep("NA",dat$ECOG)]<- NA
dat$B2mg[grep("NA",dat$B2mg)]<- NA
dat$B2mg0[grep("NA",dat$B2mg0)]<- NA
dat$LDH[grep("NA",dat$LDH)]<- NA
dat$LDH0[grep("NA",dat$LDH0)]<- NA
dat$WBC[grep("NA",dat$WBC)]<- NA
dat$HGB[grep("NA",dat$HGB)]<- NA
dat$HGB0[grep("NA",dat$HGB0)]<- NA
dat$PLT[grep("NA",dat$PLT)]<- NA
dat$Mono[grep("NA",dat$Mono)]<- NA
dat$Lym[grep("NA",dat$Lym)]<- NA
dat$interm_res[grep("NA|不详|未评估|无|N/A",dat$interm_res)]<- NA
dat$interm_res[grep("PR SD SD",dat$interm_res)]<- "SD"
dat$interm_res[grep("CR|CMR|Cru",dat$interm_res)]<- "CR"
dat$interm_res[grep("PR",dat$interm_res)]<- "PR"
dat$interm_res[grep("SD",dat$interm_res)]<- "SD"
dat$interm_res[grep(" ",dat$interm_res)] <- NA

dat$end_res[grep("CR|CMR|Cru|cr",dat$end_res)]<- "CR"
dat$end_res[grep("PR",dat$end_res)]<- "PR"
dat$end_res[grep("NA|还在化疗",dat$end_res)]<- NA
dat$end_res[grep(" ",dat$end_res)] <- NA

dat$response[grep("CR",dat$response)] <- "CR"
dat$response[grep("PD",dat$response)] <- "PD"
dat$response[grep("PR",dat$response)] <- "PR"
dat$response[grep("NA",dat$response)] <- NA
dat$response[grep(" ",dat$response)] <- NA
dat$CR[grep(" ",dat$CR)] <- NA
dat$CR[grep("NA",dat$CR)] <- NA

dat$Deauville[grep("NA|未作|未做|台湾|无",dat$Deauville)] <- NA
dat$Deauville[grep("李",dat$Deauville)] <- 0
dat$Deauville[grep("1-2分|1(第四)|1分|1月1日|1月2日|2(2)|2（2程化疗后）|2(4程)|2(3)|2（4程化疗后）",dat$Deauville)] <- 0
dat$Deauville[grep("1（4程",dat$Deauville)] <-0
dat$Deauville[grep("0（6疗程后）|1(4)|1（4程化疗后|1(6) |1（6程|1(第四)|2(4程)",dat$Deauville)]  <-0
dat$Deauville[grep("1 1",dat$Deauville)] <- 0
dat$Deauville[grep("\\1(第四)",dat$Deauville)] <- 0
dat$Deauville[grep("程\\)",dat$Deauville)] <- 0
dat$Deauville[grep("2，1",dat$Deauville)] <- 0
dat$Deauville[grep("3月",dat$Deauville)] <- 0
dat$Deauville[grep("2,",dat$Deauville)] <- 0
dat$Deauville[grep("4,|3\\.",dat$Deauville)] <- 0
dat$Deauville[grep("4，|2月|2分",dat$Deauville)] <- 0
dat$Deauville[grep("2、",dat$Deauville)] <- 0

dat$Deauville[grep("201504",dat$Deauville)] <- 9##2
dat$Deauville[grep("201",dat$Deauville)] <- 0
dat$Deauville[grep("1\\(",dat$Deauville)] <- 0
dat$Deauville[grep("2\\(",dat$Deauville)] <- 0
dat$Deauville[grep("2（",dat$Deauville)] <- 0
dat$Deauville[grep("3，",dat$Deauville)] <- 0
dat$Deauville[grep("200",dat$Deauville)] <- 0
dat$Deauville[grep("0",dat$Deauville)] <- 0
dat$Deauville[grep("5\\\\3\\\\2",dat$Deauville)] <- 0
dat$Deauville[grep("5 4 3",dat$Deauville)] <- 8 #1
dat$Deauville[grep("5\\\\3",dat$Deauville)] <- 8
dat$Deauville[grep("5",dat$Deauville)] <- 9
dat$Deauville[grep("\\\\1|\\\\2",dat$Deauville)] <- 0
dat$Deauville[grep("3（4程化疗后）|4:3分|3(4)|4 3",dat$Deauville)] <- 8
dat$Deauville[grep("\\\\4",dat$Deauville)] <- 9
dat$Deauville[grep("\\\\3",dat$Deauville)] <- 8
dat$Deauville[grep("3\\(4)",dat$Deauville)] <- 8
dat$Deauville[grep("4",dat$Deauville)] <- 9
dat$Deauville[grep("3分",dat$Deauville)] <- 8
dat$Deauville[grep("1",dat$Deauville)] <- 0
dat$Deauville[grep("3",dat$Deauville)] <- 8
dat$Deauville[grep("2",dat$Deauville)] <- 0
dat$Deauville[grep("8",dat$Deauville)] <- 1
dat$Deauville[grep("9",dat$Deauville)] <- 2

dat$Rmaintain[grep("利妥昔|CD20|美罗华",dat$Rmaintain)]  <- 1
dat$Rmaintain[grep("否|无",dat$Rmaintain)] <- 0
dat$Rmaintain[grep("是|LMZ002C|Y",dat$Rmaintain)] <- 1
dat$Rmaintain[grep("待定|未定|N/A|咨询|建议|在化疗中|NA|R|R-EPOCH|正在治疗",dat$Rmaintain)] <- NA
dat$Rmaintain[grep("N",dat$Rmaintain)] <- 0
dat$Rmaintain[grep("来那度胺",dat$Rmaintain)] <- 0
dat$Rmaintain[grep("0",dat$Rmaintain)] <- 0
dat$Rmaintain[grep("1",dat$Rmaintain)] <- 1
dat$Rmaintain[grep("后",dat$Rmaintain)] <- NA

dat$pro_LDH[grep("NA",dat$pro_LDH)] <- NA
dat$pro_ECOG[grep("NA",dat$ECOG)] <- NA
dat$relapse_time[grep("NA",dat$relapse_time)] <- NA
dat$progress_time[grep(" ",dat$progress_time)] <- NA
dat$progress_time[grep("　",dat$progress_time)] <- NA

dat <- dat[,-53]

dat$relapse_res[grep("^/|NA|未",dat$relapse_res)] <- NA
dat$relapse_res[grep("PR|第二期为pr,化疗结束无影像",dat$relapse_res)] <- "PR"
dat$relapse_res[grep("不|无|可",dat$relapse_res)] <- NA
dat$relapse_res[grep("PD",dat$relapse_res)] <- "PD"
dat$relapse_res[grep("SD",ignore.case = TRUE, dat$relapse_res)] <- "SD"
dat$relapse_res[grep("PR",ignore.case = TRUE, dat$relapse_res)] <- "PR"
dat$relapse_res[grep("CR",ignore.case = TRUE, dat$relapse_res)] <- "CR"
dat$relapse_res[grep("^N$|0",ignore.case = TRUE, dat$relapse_res)] <- NA
dat$relapse_res[grep("NR",ignore.case = TRUE, dat$relapse_res)] <- "SD"

dat$followup[grep("^202/", dat$followup)] <- "2020/01/17"
dat <- dat[-grep("245",ignore.case = TRUE, dat$followup),]
dat$followup[grep("NA|未找到住院号",ignore.case = TRUE, dat$followup)] <- NA
dat$followup <-str_replace(dat$followup,"（患者死亡，死亡日期没说，挂断电话）","")
dat$followup[grep("末次",ignore.case = TRUE, dat$followup)] <- "2019/01"
dat$followup[grep("\\.",ignore.case = TRUE, dat$followup)] <- gsub("\\.","\\/",dat$followup[grep("\\.",ignore.case = TRUE, dat$followup)])
dat$followup[grep("N/A",ignore.case = TRUE, dat$followup)] <- NA
dat$followup[grep("#",ignore.case = TRUE, dat$followup)] <- NA
dat <- dat[-grep("-",ignore.case = TRUE, dat$pfs),]
dat <- dat[-grep("1460",dat$pfs_month),]
dat <- dat[-grep("#VALUE!",dat$pfs_month),]
dat <- dat[-grep("1378",dat$pfs_month),]
dat <- dat[-grep("-",dat$OS_month),]
dat$dead_time[grep("NA",ignore.case = TRUE, dat$dead_time)] <- NA
dat <- dat[,-55]

list <- data.frame(old = c("CR","PR","SD","PD"),new = 1:4)

dat$interm_res <- sapply(dat$interm_res,function(x){
  list$new[match(x,list$old)]})
dat$end_res <- sapply(dat$end_res,function(x){
  list$new[match(x,list$old)]})
dat$response <- sapply(dat$response,function(x){
  list$new[match(x,list$old)]})
dat$relapse_res <- sapply(dat$relapse_res,function(x){
  list$new[match(x,list$old)]})



dat$trans[is.na(dat$trans)] <- 0
dat$trans[is.na(dat$grade)] <- NA
dat$trans_time[is.na(dat$trans_time)] <- 0
dat$trans_time[is.na(dat$trans)] <- NA

dat$relapes[is.na(dat$relapes)] <- 0
dat$relapes[is.na(dat$grade)] <- NA
dat$relapse_time[is.na(dat$relapse_time)] <- 0
dat$relapse_time[is.na(dat$relapes)] <- NA
dat$relapse_res[is.na(dat$relapse_res)] <- 0
dat$relapse_res[is.na(dat$relapes)] <- NA


save(dat,file = "./data/step1.Rdata")



data <- dat
f <- function(x){sum(is.na(data[x]))}
re <- list()
for (i in 1:ncol(data)) {
  print(i)
  re[[i]] <- f(i)
}
result <- as.data.frame(do.call(rbind,re))
colnames(result) <- "na_counts"
result$na_ratio <- result$na_counts/nrow(data)
rownames(result) <- colnames(data)
result$feature <- rownames(result)
result$group <- 1
e <- ggplot(result, aes(x=group, y = na_ratio))
 e + geom_violin(aes(fill = group),trim = TRUE)+geom_point(size=1)+theme_bw(base_size = 18)
 
dat <-dat%>%select(!result$feature[result$na_ratio >0.75])
 
 #---------------------------
# data_after <- read.csv("./data/rawdata.csv")
data1 <-dat

f1 <- function(x){sum(is.na(data1[x,]))}
re1 <- list()
for (i in 1:nrow(data1)) {
   print(i)
   re1[[i]] <- f1(i)
 }
result1 <- as.data.frame(do.call(rbind,re1))
colnames(result1) <- "na_counts"
result1$na_ratio <- result1$na_counts/ncol(data1)
result1$feature <- data1$No
result1$group <- 1
e <- ggplot(result1, aes(x=group, y = na_ratio))
e + geom_violin(aes(fill = group),trim = TRUE)+geom_boxplot(width = 0.2)+geom_point(size=1)+theme_bw(base_size = 18)


#cut 172samples
num <-result1$feature[result1$na_ratio <= quantile(result1$na_ratio,0.75)]
data2 <- filter(data1,No %in%num)
f1 <- function(x){sum(is.na(data2[x,]))}
re1 <- list()
for (i in 1:nrow(data2)) {
  print(i)
  re1[[i]] <- f1(i)
}
result1 <- as.data.frame(do.call(rbind,re1))
colnames(result1) <- "na_counts"
result1$na_ratio <- result1$na_counts/ncol(data2)
result1$feature <- data2$No
result1$group <- 1
e <- ggplot(result1, aes(x=group, y = na_ratio))
e + geom_violin(aes(fill = group),trim = TRUE)+geom_boxplot(width = 0.2)+geom_point(size=1)+theme_bw(base_size = 18)
f2 <- function(x){sum(is.na(data2[x]))}
re2 <- list()
for (i in 1:ncol(data2)) {
  print(i)
  re2[[i]] <- f2(i)
}
result2 <- as.data.frame(do.call(rbind,re2))
colnames(result2) <- "na_counts"
rownames(result2) <- colnames(data2)
result2$na_ratio <- result2$na_counts/nrow(data2)
result2$feature <- rownames(result2)
result2$group <- 1
e <- ggplot(result2, aes(x=group, y = na_ratio))
e + geom_violin(aes(fill = group),trim = TRUE)+geom_boxplot(width = 0.2)+geom_point(size=1)+theme_bw()
# dat_filter <- data2[,-c(2,11,12,25,27,30:37,39)]
dat_filter <-data2
write.csv(dat_filter,"./data/rawdata_filter.csv",row.names = F)
save(dat_filter,file ="./data/step2_filter_sample_and_features.Rdata")


# options(stringsAsFactors = F)
# dat <- read.csv("./data/rawdata_filter.csv")
# dat[grep(0.1,dat$Ki67),]$Ki67 <- 10
# dat[grep("^0.3",dat$Ki67),]$Ki67 <- 30
# dat[grep("^20-30",dat$Ki67),]$Ki67 <- 30
# dat[grep("^25-50",dat$Ki67),]$Ki67 <- 50
# dat[grep("^30-40",dat$Ki67),]$Ki67 <- 40
# dat[grep("^4050",dat$Ki67),]$Ki67 <- 50
# dat[grep("\\+\\+",dat$Ki67),]$Ki67 <- 65
# dat[grep("\\+",dat$Ki67),]$Ki67 <- 35
# total <- read.csv("./data/FLTOTAl.csv")
# total <- total[,c(1,9)]
# dat <- merge(dat,total,by="NO")
# dat <- dat[,-5]
# dat <- dat[!(is.na(dat$Ki67) & is.na(dat$SUVmax)),]
# 
# dat <- dat[!(is.na(dat$stage) & is.na(dat$Lymph.node)),]
# dat <- dat[!(is.na(dat$extranodal) & is.na(dat$spleen) & is.na(dat$BM)),]
# dat <- dat[!(is.na(dat$LDH.U.L.) & is.na(dat$B2.MG.ng.L.)),]
# dat <- dat[!(is.na(dat$ECOG)),]
# dat <- select(dat,!HSCT)
# save(dat,file="./data/step1.Rdata")
#  fill <- function(x,y){
#    re <- lapply(which(is.na(y)),function(i){
#      print(i)
#      range01 <- function(x){(x-min(x))/(max(x)-min(x))}
#      yy <- y[!is.na(y)]
#      a <- quantile(yy,range01(rank(x,na.last = TRUE))[i])
#      return(a)
#    })
#   }
# dat$SUVmax[is.na(dat$SUVmax)] <- unlist(fill(dat$Ki67,dat$SUVmax))
