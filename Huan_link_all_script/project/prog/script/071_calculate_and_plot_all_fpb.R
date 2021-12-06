library("survival")
library("survminer")
library(randomForest)
library(stringr)
library(dplyr)
library(Hmisc)


setwd("/home/huanhuan/project/prog/output/")
load("/home/huanhuan/project/prog/data/01_add_age_raw_pfs_os_filter_grade.Rdata")
dat$new_pod_total <-dat$pod_total
dat$new_pod_total[grep("1|2",dat$new_pod_total)] <- 0
dat$new_pod_total[grep("3|4",dat$new_pod_total)] <- 1
# dat <-filter(dat,grade=="3a")
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

dat$diff_f1 <-dat$FLIPI1_count -dat$FLIPI1
dat$diff_f2 <-dat$FLIPI1_count -dat$FLIPI2
#-------------------------------
dat$FLIPI1_count_re =dat$FLIPI1_count
dat$FLIPI1_count_re[grep("0|1",dat$FLIPI1_count)] <- "Low"
dat$FLIPI1_count_re[grep("2",dat$FLIPI1_count)] <- "Intermediate"
dat$FLIPI1_count_re[grep("3|4|5",dat$FLIPI1_count)] <- "High"


fit <- survfit(Surv(pfs_month_new, pro_status) ~ FLIPI1_count_re, data=dat)
pdf("./figure/071_all_FLIPI1_count_pfs.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)",
                 xlim=c(0,120))
print(p1) 
dev.off()

fit <- survfit(Surv(os_month_new, dead) ~ FLIPI1_count_re, data=dat)
# dat<-filter(dat,grade=="3a")
pdf("./figure/071_all_FLIPI1_count_os.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = " Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()