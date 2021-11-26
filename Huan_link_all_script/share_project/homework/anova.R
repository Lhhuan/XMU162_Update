library(readxl)
library(knitr)
library(car)

setwd("/home/huanhuan/homework/")
org<-read_excel("HOME.xlsx")
org2<-read.table("home.txt",header=T)
con<-org2$con
time<-org2$time
type<-org2$type
fit<-aov(con ~ time*type)
summary(fit)
fit2<-aov()