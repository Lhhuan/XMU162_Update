library("survival")
library("survminer")
library(randomForest)
library(stringr)
library(dplyr)
library(Hmisc)

library("survival")
library("survminer")
dat <- read.csv("/home/huanhuan/project/prog/data/rawdata_filter.csv",sep=",")
dat <- read.csv("/home/huanhuan/project/prog/data/data.csv",sep=",")
# test$FLIPI2_re[grep("1",test$FLIPI2)] <- 0
# test$FLIPI2_re[grep("2",test$FLIPI2)] <- 1
# test$FLIPI2_re[grep("3|4|5",test$FLIPI2)] <- 2



# fit <- survfit(Surv(pfs_month_new, pro_status) ~ FLIPI1_re, data=dat)
# pdf("./figure/test/aLL_FLIPI1_pfs.pdf")
# p1 <- ggsurvplot(fit,
#                   pval = TRUE,
#                  legend.title="Pod total new",
#                  title="3a: Pfs",
#                  xlab = " Time (Months)",
#                  xlim=c(0,120))
# print(p1) 
# dev.off()
setwd("/home/huanhuan/project/prog/script/hh/")
fit <- survfit(Surv(OS_month, dead) ~ FLIPI1, data=dat)
pdf("./figure/test/aLL_FLIPI1_os_ori.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = " Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()


setwd("/home/huanhuan/project/prog/script/hh/")
fit <- survfit(Surv(OS_month, dead) ~ FLIPI2, data=dat)
pdf("./figure/test/aLL_FLIPI2_os_ori.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = " Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()