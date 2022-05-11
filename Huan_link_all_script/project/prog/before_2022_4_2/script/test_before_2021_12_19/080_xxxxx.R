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
load("07_add_model_index_age_raw_pfs_os_filter_grade_mergrPod_3A.Rdata")
pod_total0=which(dat$new_pod_total==0)
# set.seed(112231) #/TOP1
# set.seed(1124)
set.seed(168)
test_number0 <-sample(x=pod_total0, round(length(pod_total0)*1/3),replace = F)
pod_total1=which(dat$new_pod_total!=0)
# set.seed(112231)
# set.seed(1124)
set.seed(168)
test_number1 <-sample(x=pod_total1, round(length(pod_total1)*1/3),replace = F)

test_set_number = c(test_number0,test_number1)
train_set_number =setdiff(1:nrow(dat),test_set_number)
#------------------------------------------------------------------new_cutoff
dat$B2mg_c <-dat$B2mg
dat$B2mg_c[dat$B2mg_c <=3.4]=0
dat$B2mg_c[dat$B2mg_c >3.4 ]=1

dat$LN_num_c <-dat$LN_num
dat$LN_num_c[dat$LN_num_c <=10]=0
dat$LN_num_c[dat$LN_num_c >10]=1

dat$LDH_c <-dat$LDH
dat$LDH_c[dat$LDH_c <=270]=0
dat$LDH_c[dat$LDH_c >270]=1

dat$age_raw_c <-dat$age_raw
dat$age_raw_c[dat$age_raw_c <=60]=0
dat$age_raw_c[dat$age_raw_c >60 ]=1

dat$Lym_Mono_c <-dat$Lym_Mono
dat$Lym_Mono_c[dat$Lym_Mono_c <=5]=0
dat$Lym_Mono_c[dat$Lym_Mono_c >5 ]=1

dat$Ki.67_c <-dat$Ki.67
dat$Ki.67_c[dat$Ki.67_c <=70.00]=0
dat$Ki.67_c[dat$Ki.67_c >70.00 ]=1

dat$SPD_c <-dat$SPD
dat$SPD_c[dat$SPD_c <=20]=0
dat$SPD_c[dat$SPD_c >20]=1

dat$SUVmax_c <-dat$SUVmax
dat$SUVmax_c[dat$SUVmax_c <=2]=0
dat$SUVmax_c[dat$SUVmax_c >2 ]=1

dat$HGB_c <-dat$HGB
dat$HGB_c[dat$HGB_c >=120 ]=0
dat$HGB_c[dat$HGB_c <120]=1

#-------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------FLIPI
test=dat[test_set_number,]

fit <- survfit(Surv(pfs_month_new, pro_status) ~ FLIPI1_count_re, data=test)
pdf("./figure/test/069_pre_survival_3a_pod_total_new_sum_pfs_FLIPI1.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)",
                 xlim=c(0,120))
print(p1) 
dev.off()

fit <- survfit(Surv(os_month_new, dead) ~ FLIPI1_count_re, data=test)
pdf("./figure/test/069_pre_survival_3a_pod_total_new_sum_os_FLIPI1.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = " Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()
#--------------------------------------------------------------------------------2
fit <- survfit(Surv(pfs_month_new, pro_status) ~ FLIPI2_count_re, data=test)
pdf("./figure/test/069_pre_survival_3a_pod_total_new_sum_pfs_FLIPI2.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)",
                 xlim=c(0,120))
print(p1) 
dev.off()

fit <- survfit(Surv(os_month_new, dead) ~ FLIPI2_count_re, data=test)
pdf("./figure/test/069_pre_survival_3a_pod_total_new_sum_os_FLIPI2.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = " Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()
#--------------------------------------------------------------------------------------------adjust1 weight_a
dat$LDH_c_f <- dat$LDH_c * 24.53
dat$LN_num_c_f <-dat$LN_num_c *13.21
dat$Lym_Mono_c_f <-dat$Lym_Mono_c*11.79
dat$SPD_c_f<- dat$SPD_c *11.32
dat$B2mg_c_f <- dat$B2mg_c*9.91
dat$Bsym_c_f <- dat$Bsym*8.96
dat$age_raw_c_f <- dat$age_raw_c*6.60
dat$HGB_c_f<- dat$HGB_c *6.13
dat$Ki.67_c_f<- dat$Ki.67_c *2.83
dat$BM_c_f<-dat$BM*2.83
dat$SUVmax_c_f<- dat$SUVmax_c *1.89
#-------------
dat_s<-dat[,c("B2mg_c_f","LN_num_c_f","LDH_c_f","age_raw_c_f","Lym_Mono_c_f","HGB_c_f","Ki.67_c_f","SPD_c_f","SUVmax_c_f","Bsym_c_f","BM_c_f")]
dat$predict_Podtotal_new <- base::rowSums(dat_s, na.rm = TRUE)
dat$Predict_Pod_total_f <-dat$predict_Podtotal_new
test=dat[test_set_number,]
train = dat[train_set_number,]
train_new_pod <-table(train$new_pod_total)%>%as.data.frame()
colnames(train_new_pod ) <-c("pod_total_new","number")
train_new_pod$ratio <-train_new_pod$number/sum(train_new_pod$number)


train <-train[order(train$predict_Podtotal_new),]
cutoff_y = round(train[271,"Predict_Pod_total_f"])
test$Predict_Pod_total_f[test$Predict_Pod_total_f<=cutoff_y ]=0
test$Predict_Pod_total_f[test$Predict_Pod_total_f>cutoff_y]=1


library("survival")
library("survminer")

fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=test)
pdf("./figure/test/069_pre_survival_3a_pod_total_new_sum_pfs_adjust1_new_cutoff_y_32.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)",
                 xlim=c(0,120))
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=test)
pdf("./figure/test/069_pre_survival_3a_pod_total_new_sum_Os_adjust1_new_cutoff_y_32.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()

#---------------------------------------------------------Ratio_adjust2
dat$LDH_c_f <- dat$LDH_c * 10
dat$LN_num_c_f <-dat$LN_num_c *6
dat$Lym_Mono_c_f <-dat$Lym_Mono_c*5
dat$SPD_c_f<- dat$SPD_c *5
dat$B2mg_c_f <- dat$B2mg_c*4
dat$Bsym_c_f <- dat$Bsym*4
dat$age_raw_c_f <- dat$age_raw_c*3
dat$HGB_c_f<- dat$HGB_c *3
dat$Ki.67_c_f<- dat$Ki.67_c *1
dat$BM_c_f<-dat$BM*1
dat$SUVmax_c_f<- dat$SUVmax_c *1
#-------------
dat_s<-dat[,c("B2mg_c_f","LN_num_c_f","LDH_c_f","age_raw_c_f","Lym_Mono_c_f","HGB_c_f","Ki.67_c_f","SPD_c_f","SUVmax_c_f","Bsym_c_f","BM_c_f")]
dat$predict_Podtotal_new <- base::rowSums(dat_s, na.rm = TRUE)
dat$Predict_Pod_total_f <-dat$predict_Podtotal_new
test=dat[test_set_number,]
train = dat[train_set_number,]
train_new_pod <-table(train$new_pod_total)%>%as.data.frame()
colnames(train_new_pod ) <-c("pod_total_new","number")
train_new_pod$ratio <-train_new_pod$number/sum(train_new_pod$number)


train <-train[order(train$predict_Podtotal_new),]
cutoff_y = round(train[271,"Predict_Pod_total_f"])
test$Predict_Pod_total_f[test$Predict_Pod_total_f<=cutoff_y ]=0
test$Predict_Pod_total_f[test$Predict_Pod_total_f>cutoff_y]=1


fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=test)
pdf("./figure/test/069_pre_survival_3a_pod_total_new_sum_pfs_adjust2_new_cutoff_y14.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)",
                 xlim=c(0,120))
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=test)
pdf("./figure/test/069_pre_survival_3a_pod_total_new_sum_Os_adjust2_new_cutoff_y14.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()


#------------------------------------------------------------------------------------------------------------------------------------------------------------adjust3
dat$LDH_c_f <- dat$LDH_c * 4
dat$LN_num_c_f <-dat$LN_num_c *3
dat$Lym_Mono_c_f <-dat$Lym_Mono_c*3
dat$SPD_c_f<- dat$SPD_c *3
dat$B2mg_c_f <- dat$B2mg_c*2
dat$Bsym_c_f <- dat$Bsym*2
dat$age_raw_c_f <- dat$age_raw_c*2
dat$HGB_c_f<- dat$HGB_c *2
dat$Ki.67_c_f<- dat$Ki.67_c *1
dat$BM_c_f<-dat$BM*1
dat$SUVmax_c_f<- dat$SUVmax_c *1

dat_s<-dat[,c("B2mg_c_f","LN_num_c_f","LDH_c_f","age_raw_c_f","Lym_Mono_c_f","HGB_c_f","Ki.67_c_f","SPD_c_f","SUVmax_c_f","Bsym_c_f","BM_c_f")]
dat$predict_Podtotal_new <- base::rowSums(dat_s, na.rm = TRUE)
dat$Predict_Pod_total_f <-dat$predict_Podtotal_new
test=dat[test_set_number,]
train = dat[train_set_number,]
train_new_pod <-table(train$new_pod_total)%>%as.data.frame()
colnames(train_new_pod ) <-c("pod_total_new","number")
train_new_pod$ratio <-train_new_pod$number/sum(train_new_pod$number)


train <-train[order(train$predict_Podtotal_new),]
cutoff_y = round(train[271,"Predict_Pod_total_f"])


test$Predict_Pod_total_f[test$Predict_Pod_total_f<=cutoff_y ]=0
test$Predict_Pod_total_f[test$Predict_Pod_total_f>cutoff_y]=1



fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=test)
pdf("./figure/test/069_pre_survival_3a_pod_total_new_sum_pfs_adjust3_new_cutoff_y7.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)",
                 xlim=c(0,120))
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=test)
pdf("./figure/test/069_pre_survival_3a_pod_total_new_sum_Os_adjust3_new_cutoff_y7.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()

#--------------------------------------------------------------------------------------------------------------------------------------------------------------adjust4

dat$LDH_c_f <- dat$LDH_c * 5
dat$LN_num_c_f <-dat$LN_num_c *3
dat$Lym_Mono_c_f <-dat$Lym_Mono_c*3
dat$SPD_c_f<- dat$SPD_c *2
dat$B2mg_c_f <- dat$B2mg_c*2
dat$Bsym_c_f <- dat$Bsym*2
dat$age_raw_c_f <- dat$age_raw_c*1
dat$HGB_c_f<- dat$HGB_c *1
dat$Ki.67_c_f<- dat$Ki.67_c *1
dat$BM_c_f<-dat$BM*1
dat$SUVmax_c_f<- dat$SUVmax_c *1
dat_s<-dat[,c("B2mg_c_f","LN_num_c_f","LDH_c_f","age_raw_c_f","Lym_Mono_c_f","HGB_c_f","Ki.67_c_f","SPD_c_f","SUVmax_c_f","Bsym_c_f","BM_c_f")]
dat$predict_Podtotal_new <- base::rowSums(dat_s, na.rm = TRUE)
dat$Predict_Pod_total_f <-dat$predict_Podtotal_new
test=dat[test_set_number,]
train = dat[train_set_number,]
train_new_pod <-table(train$new_pod_total)%>%as.data.frame()
colnames(train_new_pod ) <-c("pod_total_new","number")
train_new_pod$ratio <-train_new_pod$number/sum(train_new_pod$number)


train <-train[order(train$predict_Podtotal_new),]
cutoff_y = round(train[271,"Predict_Pod_total_f"])
test$Predict_Pod_total_f[test$Predict_Pod_total_f<=cutoff_y ]=0
test$Predict_Pod_total_f[test$Predict_Pod_total_f>cutoff_y]=1
fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=test)
pdf("./figure/test/069_pre_survival_3a_pod_total_new_sum_pfs_adjust4_new_cutoff_y7.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = "Time (Months)",
                 xlim=c(0,120))
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=test)
pdf("./figure/test/069_pre_survival_3a_pod_total_new_sum_Os_adjust4_new_cutoff_y7.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
