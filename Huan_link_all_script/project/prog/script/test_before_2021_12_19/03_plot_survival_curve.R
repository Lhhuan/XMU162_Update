library("survival")
library("survminer")
library(randomForest)
library(stringr)
library(dplyr)
library(Hmisc)

setwd("/home/huanhuan/project/prog/output/")
load("../data/01_add_age_raw_pfs_os_filter_grade.Rdata")
dat$new_pod_total <-dat$pod_total
dat$new_pod_total[grep("1|2",dat$new_pod_total)] <- 0
dat$new_pod_total[grep("3|4",dat$new_pod_total)] <- 1
dat$SPD <-as.numeric(dat$SPD)

dat$pod_total3 <-dat$pod_total
dat$pod_total3[grep("1|2",dat$pod_total3)] <- 0
# datk$diff <-abs(datk$Ki.67 -as.numeric(datk$Ki67))

#-------------------------------
prog <- filter(dat,grade=="3a")
fit <- survfit(Surv(pfs_month_new, pro_status) ~ new_pod_total, data=prog)
pdf("./figure/03_ori_survival_3a_new_pod_total_pfs.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod_total",
                 title="3a: Pfs",
                 xlab = " Time (Months)")
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ new_pod_total, data=prog)
pdf("./figure/03_ori_survival_3a_new_pod_total_Os.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total",
                 title="3a: Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()
#-----------------------------------

#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ pod_total, data=prog)
pdf("./figure/03_ori_survival_3a_pod_total_Os.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total",
                 title="3a: Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()
#-----------
fit <- survfit(Surv(os_month_new, dead) ~ pod_total3, data=prog)
pdf("./figure/03_ori_survival_3a_pod_total3_Os.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total",
                 title="3a: Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()
#------------------
fit <- survfit(Surv(pfs_month_new, pro_status) ~ pod_total3, data=prog)
pdf("./figure/03_ori_survival_3a_pod_total3_pfs.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total",
                 title="3a: pfs",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()
#---------------------------------------------------
prog <- filter(dat,grade=="3b")
fit <- survfit(Surv(os_month_new, dead) ~ new_pod_total, data=prog)
pdf("./figure/03_ori_survival_3b_new_pod_total_Os.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total",
                 title="3b: Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()
#----------------------------
prog <- filter(dat,grade=="0")
fit <- survfit(Surv(pfs_month_new, pro_status) ~ new_pod_total, data=prog)
pdf("./figure/03_ori_survival_0_new_pod_total_pfs.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod_total",
                 title="0: Pfs",
                 xlab = " Time (Months)")
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ new_pod_total, data=prog)
pdf("./figure/03_ori_survival_0_new_pod_total_Os.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total",
                 title="0: Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()