library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)
library(cowplot)
library(ggpubr)
library(Hmisc)
library("survival")
library("survminer")

setwd("/home/huanhuan/project/prog/output/")
dat <- read.csv("01_add_age_raw_pfs_os_filter_grade_mergrPod_3A.txt",sep="\t")

dat$age_s=NA
dat$age_s[dat$age_raw<=60]=0
dat$age_s[dat$age_raw>60]=1

dat$stage_s=NA
dat$stage_s[dat$stage<3]=0
dat$stage_s[dat$stage>=3]=1
dat$stage_s[is.na(dat$stage_s) &dat$new_pod_total==0] =0
dat$stage_s[is.na(dat$stage_s) &dat$new_pod_total==1] =1


dat$HGB_s =NA
dat$HGB_s[dat$HGB <120]=1
dat$HGB_s[dat$HGB >=120]=0

dat$LDH_s =NA
dat$LDH_s[dat$LDH <=245]= 0
dat$LDH_s[dat$LDH >245]= 1 #----
dat$LDH_s[is.na(dat$LDH_s)]=dat$LDH0[is.na(dat$LDH_s)]
dat$LDH_s[is.na(dat$LDH_s)&dat$new_pod_total==0]=0
dat$LDH_s[is.na(dat$LDH_s)&dat$new_pod_total==1]=1


dat$LN_num_s =NA
dat$LN_num_s[dat$LN_num <=4]= 0
dat$LN_num_s[dat$LN_num >4]= 1
dat$LN_num_s[is.na(dat$LN_num_s)&dat$stage_s==1]=1
dat$LN_num_s[is.na(dat$LN_num_s)&dat$stage_s==0]=0

dat$FLIPI1_count <- base::rowSums(dat[,c("age_s","stage_s","HGB_s","LDH_s","LN_num_s")], na.rm = TRUE)

# dat$diff_f1 <-dat$FLIPI1_count -dat$FLIPI1
# dat$diff_f2 <-dat$FLIPI1_count -dat$FLIPI2
#-------------------------------
dat$FLIPI1_count_re =dat$FLIPI1_count
dat$FLIPI1_count_re[grep("0|1",dat$FLIPI1_count)] <- "Low"
dat$FLIPI1_count_re[grep("2",dat$FLIPI1_count)] <- "Intermediate"
dat$FLIPI1_count_re[grep("3|4|5",dat$FLIPI1_count)] <- "High"


fit <- survfit(Surv(pfs_month_new, pro_status) ~ FLIPI1_count_re, data=dat)
pdf("./figure/07_3a_FLIPI1_count_pfs.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="FLIPI1",
                 title="3a: Pfs",
                 xlab = " Time (Months)",
                 xlim=c(0,120))
print(p1) 
dev.off()

fit <- survfit(Surv(os_month_new, dead) ~ FLIPI1_count_re, data=dat)
pdf("./figure/07_3a_FLIPI1_count_os.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="FLIPI1",
                 title="3a: Os",
                 xlab = " Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()
#-----------------------------------

#---------------------------------------------flipi2
dat$B2mg_s <- NA
dat$B2mg_s[dat$B2mg<=2.7]=0
dat$B2mg_s[dat$B2mg > 2.7]=1
dat$B2mg_s[is.na(dat$B2mg_s)]=dat$B2mg0[is.na(dat$B2mg_s)]
dat$B2mg_s[is.na(dat$B2mg_s)&dat$new_pod_total==1]=1
dat$B2mg_s[is.na(dat$B2mg_s)&dat$new_pod_total==0]=0


dat$LN6_s <- dat$LN6
dat$LN6_s[is.na(dat$LN6_s)&dat$new_pod_total==1]=1
dat$LN6_s[is.na(dat$LN6_s)&dat$new_pod_total==0]=0

dat$BM_s <- dat$BM
dat$BM_s[is.na(dat$BM_s)&dat$new_pod_total==1]=1
dat$BM_s[is.na(dat$BM_s)&dat$new_pod_total==0]=0
#HGB_s, age_s
dat$FLIPI2_count <- base::rowSums(dat[,c("B2mg_s","LN6_s","BM_s","HGB_s","age_s")], na.rm = TRUE)

dat$FLIPI2_count_re =dat$FLIPI2_count
dat$FLIPI2_count_re[grep("0",dat$FLIPI2_count)] <- "Low"
dat$FLIPI2_count_re[grep("1|2",dat$FLIPI2_count)] <- "Intermediate"
dat$FLIPI2_count_re[grep("3|4|5",dat$FLIPI2_count)] <- "High"
#--------------------------------------------------------

fit <- survfit(Surv(pfs_month_new, pro_status) ~ FLIPI2_count_re, data=dat)
pdf("./figure/070_3a_FLIPI2_count_pfs.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="FLIPI2",
                 title="3a: Pfs",
                 xlab = " Time (Months)",
                 xlim=c(0,120))
print(p1) 
dev.off()

fit <- survfit(Surv(os_month_new, dead) ~ FLIPI2_count_re, data=dat)
pdf("./figure/070_3a_FLIPI2_count_os.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="FLIPI2",
                 title="3a: Os",
                 xlab = " Time (Months)",
                 xlim=c(0,110))
# print(p1) 
    p1
dev.off()
# ggsave("./figure/070_3a_FLIPI2_count_os.pdf",p1)
#------------------------------------------primapi

dat$primapi_re  <- NA
dat$primapi_re[dat$B2mg>3]= "High"
dat$primapi_re[dat$B2mg<=3 &dat$BM_s==1 ]="Intermediate"
dat$primapi_re[dat$B2mg<=3 &dat$BM_s==0 ]="Low"
dat$primapi_re[is.na(dat$primapi_re) & dat$pod_total==0 ]="Low"
dat$primapi_re[is.na(dat$primapi_re) & dat$pod_total>=3 ]="High"
dat$primapi_re[is.na(dat$primapi_re)]="Intermediate"



save(dat,file ="07_add_model_index_age_raw_pfs_os_filter_grade_mergrPod_3A.Rdata")

fit <- survfit(Surv(pfs_month_new, pro_status) ~ primapi_re, data=dat)
pdf("./figure/070_3a_Primapi_re_count_pfs.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="PRIMA-PI",
                 title="3a: Pfs",
                 xlab = " Time (Months)",
                 xlim=c(0,120))
print(p1) 
dev.off()

fit <- survfit(Surv(os_month_new, dead) ~ primapi_re, data=dat)
pdf("./figure/070_3a_Primapi_re_count_os.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="PRIMA-PI",
                 title="3a: Os",
                 xlab = " Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()
#----------------------------------b2mg
fit <- survfit(Surv(pfs_month_new, pro_status) ~ X150b2mg_ldh, data=dat)
pdf("./figure/070_3a_X150b2mg_ldh_re_count_pfs.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="B2mg_LDH",
                 title="3a: Pfs",
                 xlab = " Time (Months)",
                 xlim=c(0,120))
print(p1) 
dev.off()

fit <- survfit(Surv(os_month_new, dead) ~ X150b2mg_ldh, data=dat)
pdf("./figure/070_3a_X150b2mg_ldh_count_os.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="B2mg_LDH",
                 title="3a: Os",
                 xlab = " Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()




#-----------------------------------------------

#-----------------------------------
