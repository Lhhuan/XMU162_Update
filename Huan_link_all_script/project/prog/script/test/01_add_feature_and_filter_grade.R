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
dat <- read.csv("rawdata_filter.csv")
dat <- merge(dat,data[,c("No","birth","dead_time","pro_time","pro_status","age_raw")],by="No")
dat$diagnosis<-as.Date(dat$diagnosis)
dat$dead_time<-as.Date(dat$dead_time)
dat$pro_time<-as.Date(dat$pro_time)
dat$followup<-as.Date(dat$followup)
# pfs <- function(x){
#   ifelse(x[1,58] == 0,a = difftime(x[1,50],x[1,38],units = "days"),a = difftime(x[1,57],x[1,38],units = "days"))
#   return(a)}

pfs <- function(x){
  if(x[1,58] == 0){
    a = difftime(x[1,50],x[1,38],units = "days")}else{
      a = difftime(x[1,57],x[1,38],units = "days")}
  return(a)}

re <- lapply(1:nrow(dat), function(i){
  print(i)
  x = dat[i,]
  a <- round(pfs(x)/30)
  return(a)
}) 
dat$pfs_month_new <- unlist(re)
  
os <- function(x){
  if(x[1,51] == 0){
    a = difftime(x[1,50],x[1,38],units = "days")}else{
      a = difftime(x[1,56],x[1,38],units = "days")}
  return(a)}

re1 <- lapply(1:nrow(dat), function(i){
  print(i)
  x = dat[i,]
  a <- round(os(x)/30)
  return(a)
}) 
dat$os_month_new <- unlist(re1)

dat$grade[grep("3$",dat$grade)] <-"3a"
dat <-dat[!is.na(dat$grade),]
dat <-dat[-grep("-7",dat$os_month_new),]
save(dat,file="01_add_age_raw_pfs_os_filter_grade.Rdata")

library("survival")
library("survminer")

#-----------------
# dat$dead <-dat$dead +1
# dat$pro_status <-dat$pro_status +1
fit <- survfit(Surv(os_month_new, dead) ~ grade, data=dat)
pdf("./figure/01_survival/01_survival_os.pdf")
p1 <- ggsurvplot(fit,
                 legend.title="Grade",
                 title="Os",
                 xlab = "Time (Months)")
print(p1) 
dev.off()


fit <- survfit(Surv(pfs_month_new, pro_status) ~ grade, data=dat)
pdf("./figure/01_survival/01_survival_pfs.pdf")
p1 <- ggsurvplot(fit,
                 legend.title="Grade",
                 title="Pfs",
                 xlab = " Time (Months)")
print(p1) 
dev.off()

