library("survival")
library("survminer")
library(randomForest)
library(stringr)
library(dplyr)
library(Hmisc)

setwd("/home/huanhuan/project/prog/script/hh")
load("/home/huanhuan/project/prog/data/01_add_age_raw_pfs_os_filter_grade.Rdata")

dat$FLIPI1_re <- dat$FLIPI1
dat$FLIPI2_re <-dat$FLIPI2
#--------------------------------------------------------------------------------1
# test$FLIPI1_re[grep("1",test$FLIPI1)] <- 0
# # test$FLIPI1_re[grep("2",test$FLIPI1)] <- 1
# test$FLIPI1_re[grep("2|3|4|5",test$FLIPI1)] <- 1

dat$FLIPI1_re[grep("0|1",dat$FLIPI1)] <- "Low"
dat$FLIPI1_re[grep("2",dat$FLIPI1)] <- "Intermediate"
dat$FLIPI1_re[grep("3|4|5",dat$FLIPI1)] <- "High"

fit <- survfit(Surv(pfs_month_new, pro_status) ~ FLIPI1_re, data=dat)
pdf("./figure/test/aLL_FLIPI1_pfs.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)",
                 xlim=c(0,120))
print(p1) 
dev.off()

fit <- survfit(Surv(OS_month, dead) ~ FLIPI1, data=dat)
pdf("./figure/test/aLL_FLIPI1_os.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = " Time (Months)",
                 xlim=c(0,200))
print(p1) 
dev.off()
#-----------------
dat$FLIPI2_re[grep("0|1",dat$FLIPI2)] <- "Low"
dat$FLIPI2_re[grep("2",dat$FLIPI2)] <- "Intermediate"
dat$FLIPI2_re[grep("3|4|5",dat$FLIPI2)] <- "High"

fit <- survfit(Surv(pfs_month_new, pro_status) ~ FLIPI2_re, data=dat)
pdf("./figure/test/aLL_FLIPI2_pfs.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)",
                 xlim=c(0,120))
print(p1) 
dev.off()

fit <- survfit(Surv(os_month_new, dead) ~ FLIPI2_re, data=dat)
pdf("./figure/test/aLL_FLIPI2_os.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = " Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()

an3 <- filter(dat,grade=="3a")
dat3 <- read.csv("/home/huanhuan/project/prog/output/01_add_age_raw_pfs_os_filter_grade_mergrPod_3A.txt",sep="\t")
setdiff(dat3$No,an3$No)
setdiff(an3$No,dat3$No)
filter(an3,No==232)