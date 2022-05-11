library("survival")
library("survminer")
# library(randomForest)
library(stringr)
library(dplyr)
library(Hmisc)
setwd("/home/huanhuan/project/prog/output/")
dat <- read.csv("06_XGBClassifier_predict_3a_not_fill_CV_num_class.txt",na.strings = "",sep="\t")

dat$pred_new_pod_total <-dat$predict_pod_total
dat$pred_new_pod_total[grep("1|2",dat$pred_new_pod_total)] <- 0
dat$pred_new_pod_total[grep("3|4",dat$pred_new_pod_total)] <- 1

fit <- survfit(Surv(pfs_month_new, pro_status) ~ pred_new_pod_total, data=dat)
pdf("./figure/07_pred_new_pod_total_3a_pfs_CV_c.pdf")
p1 <- ggsurvplot(fit,
                pval = TRUE,
                 legend.title="Pod_total",
                 title="Pfs",
                 xlab = " Time (Months)")
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ pred_new_pod_total, data=dat)
pdf("./figure/07_pred_new_pod_total_3a_Os_CV_C.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total",
                 title="Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()

#-------------------predict_pod_total
fit <- survfit(Surv(pfs_month_new, pro_status) ~ predict_pod_total, data=dat)
pdf("./figure/07_pred_pod_total_3a_pfs_CV_c.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod_total",
                 title="Pfs",
                 xlab = " Time (Months)")
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ predict_pod_total, data=dat)
pdf("./figure/07_pred_pod_total_3a_Os_CV_c.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total",
                 title="Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()

#------
fit <- survfit(Surv(pfs_month_new, pro_status) ~ FLIPI1, data=dat)
pdf("./figure/07_FLIPI1_3a_pfs_CV.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod_total",
                 title="Pfs",
                 xlab = " Time (Months)")
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ FLIPI1, data=dat)
pdf("./figure/07_FLIPI1_3a_Os_CV.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total",
                 title="Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()

