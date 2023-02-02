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
load("09_test_train_dataset.Rdata")
dat$primapi_re_n <- NA
dat$primapi_re_n[dat$primapi_re =="Low"]=0
dat$primapi_re_n[dat$primapi_re =="Intermediate"]=1
dat$primapi_re_n[dat$primapi_re =="High"]=2
#------------------------------------------------------------------new_cutoff
dat$Lym_Mono_c <-NA
dat$Lym_Mono_c[dat$Lym_Mono <=10]=0
dat$Lym_Mono_c[dat$Lym_Mono >10 ]=1
dat$Lym_Mono_c[is.na(dat$Lym_Mono_c)&dat$new_pod_total==1]=1
dat$Lym_Mono_c[is.na(dat$Lym_Mono_c)&dat$new_pod_total==0]=0

dat$B2mg_c <-NA
dat$B2mg_c[dat$B2mg <=2.7]=0
dat$B2mg_c[dat$B2mg >2.7 ]=1
dat$B2mg_c[is.na(dat$B2mg_c)&dat$new_pod_total==1]=1
dat$B2mg_c[is.na(dat$B2mg_c)&dat$new_pod_total==0]=0

dat$SPD_c <-NA
dat$SPD_c[dat$SPD <=20]=0
dat$SPD_c[dat$SPD >20]=1
dat$SPD_c[is.na(dat$SPD_c)&dat$new_pod_total==1]=1
dat$SPD_c[is.na(dat$SPD_c)&dat$new_pod_total==0]=0

dat$LDH_c <-NA
dat$LDH_c[dat$LDH <=245]=0
dat$LDH_c[dat$LDH >245]=1
dat$LDH_c[is.na(dat$LDH_c)&dat$new_pod_total==1]=1
dat$LDH_c[is.na(dat$LDH_c)&dat$new_pod_total==0]=0

dat$LN_num_c <-NA
dat$LN_num_c[dat$LN_num <=4]=0
dat$LN_num_c[dat$LN_num >4]=1
dat$LN_num_c[is.na(dat$LN_num_c)&dat$new_pod_total==1]=1
dat$LN_num_c[is.na(dat$LN_num_c)&dat$new_pod_total==0]=0
#----------
dat$SUVmax_c <-NA
dat$SUVmax_c[dat$SUVmax <=2]=0
dat$SUVmax_c[dat$SUVmax >2 ]=1
dat$SUVmax_c[is.na(dat$SUVmax_c)&dat$new_pod_total==1]=1
dat$SUVmax_c[is.na(dat$SUVmax_c)&dat$new_pod_total==0]=0

dat$age_raw_c <-dat$age_s

dat$Ki.67_c <-NA 
dat$Ki.67_c[dat$Ki.67 <=20.00]=0
dat$Ki.67_c[dat$Ki.67 >20.00 ]=1
dat$Ki.67_c[is.na(dat$Ki.67_c)&dat$new_pod_total==1]=1
dat$Ki.67_c[is.na(dat$Ki.67_c)&dat$new_pod_total==0]=0

#------
dat$HGB_c =dat$HGB_s
# dat$HGB_c <-NA
# dat$HGB_c[dat$HGB >=120 ]=0
# dat$HGB_c[dat$HGB <120]=1
# dat$HGB_c[is.na(dat$HGB_c)&dat$new_pod_total==1]=1
# dat$HGB_c[is.na(dat$HGB_c)&dat$new_pod_total==0]=0


dat$BM_c <-dat$BM
dat$BM_c[is.na(dat$BM_c)&dat$new_pod_total==1]=1
dat$BM_c[is.na(dat$BM_c)&dat$new_pod_total==0]=0

dat$Bsym_c <- dat$Bsym
dat$Bsym_c[is.na(dat$Bsym_c)&dat$new_pod_total==1]=1
dat$Bsym_c[is.na(dat$Bsym_c)&dat$new_pod_total==0]=0
#-------------------------------------------------------------------------------------------------------------
test <-dat[test_set_number,]
im <-read.table("/home/huanhuan/project/prog/output/092_feature_ratio.txt",header=T,sep="\t")

# test$sum <-base::rowSums(test[,c("Lym_Mono","B2mg","SPD","LDH","LN_num","SUVmax","Ki.67","HGB","BM","Bsym")])
# base::rowSums(dat_s, na.rm = TRUE)
#-----------------------------------------------adjust1
test$Lym_Mono_c_f <-test$Lym_Mono_c*21.40
test$B2mg_c_f <- test$B2mg_c*18.78
test$SPD_c_f<- test$SPD_c *17.47
test$LDH_c_f <- test$LDH_c * 15.28
test$HGB_c_f<- test$HGB_c *6.99
test$Bsym_c_f <- test$Bsym_c*6.55
test$LN_num_c_f <-test$LN_num_c *6.11
test$BM_c_f<-test$BM_c*3.06
test$SUVmax_c_f<- test$SUVmax_c *2.62
test$age_raw_c_f <- test$age_raw_c*1.31
test$Ki.67_c_f<- test$Ki.67_c *0.44

test$sum_score<-base::rowSums(test[,c("B2mg_c_f","LN_num_c_f","LDH_c_f","age_raw_c_f","Lym_Mono_c_f","HGB_c_f","Ki.67_c_f","SPD_c_f","SUVmax_c_f","Bsym_c_f","BM_c_f")])
test[is.na(test$sum_score),]
test$adjust1 <- test$sum_score
# adjust1 <-test
#-----------------------------------------------adjust2
test$Lym_Mono_c_f <-test$Lym_Mono_c*10
test$B2mg_c_f <- test$B2mg_c*9
test$SPD_c_f<- test$SPD_c *8
test$LDH_c_f <- test$LDH_c * 7
test$HGB_c_f<- test$HGB_c *3
test$Bsym_c_f <- test$Bsym_c*3
test$LN_num_c_f <-test$LN_num_c *3
test$BM_c_f<-test$BM_c*1
test$SUVmax_c_f<- test$SUVmax_c *1
test$age_raw_c_f <- test$age_raw_c*1
test$Ki.67_c_f<- test$Ki.67_c *1

test$sum_score<-base::rowSums(test[,c("B2mg_c_f","LN_num_c_f","LDH_c_f","age_raw_c_f","Lym_Mono_c_f","HGB_c_f","Ki.67_c_f","SPD_c_f","SUVmax_c_f","Bsym_c_f","BM_c_f")])
# test[is.na($sum_score),]
test$adjust2 <- test$sum_score
#-----------------------------------------------adjust3
test$Lym_Mono_c_f <-test$Lym_Mono_c*5
test$B2mg_c_f <- test$B2mg_c*5
test$SPD_c_f<- test$SPD_c *4
test$LDH_c_f <- test$LDH_c * 4
test$HGB_c_f<- test$HGB_c *2
test$Bsym_c_f <- test$Bsym_c*2
test$LN_num_c_f <-test$LN_num_c *2
test$BM_c_f<-test$BM_c*1
test$SUVmax_c_f<- test$SUVmax_c *1
test$age_raw_c_f <- test$age_raw_c*1
test$Ki.67_c_f<- test$Ki.67_c *1

test$sum_score<-base::rowSums(test[,c("B2mg_c_f","LN_num_c_f","LDH_c_f","age_raw_c_f","Lym_Mono_c_f","HGB_c_f","Ki.67_c_f","SPD_c_f","SUVmax_c_f","Bsym_c_f","BM_c_f")])
# test[is.na($sum_score),]
test$adjust3 <- test$sum_score
#-----------------------------------------------adjust4
test$Lym_Mono_c_f <-test$Lym_Mono_c*3
test$B2mg_c_f <- test$B2mg_c*3
test$SPD_c_f<- test$SPD_c *3
test$LDH_c_f <- test$LDH_c * 3
test$HGB_c_f<- test$HGB_c *2
test$Bsym_c_f <- test$Bsym_c*2
test$LN_num_c_f <-test$LN_num_c *2
test$BM_c_f<-test$BM_c*1
test$SUVmax_c_f<- test$SUVmax_c *1
test$age_raw_c_f <- test$age_raw_c*1
test$Ki.67_c_f<- test$Ki.67_c *1

test$sum_score<-base::rowSums(test[,c("B2mg_c_f","LN_num_c_f","LDH_c_f","age_raw_c_f","Lym_Mono_c_f","HGB_c_f","Ki.67_c_f","SPD_c_f","SUVmax_c_f","Bsym_c_f","BM_c_f")])
# test[is.na($sum_score),]
test$adjust4 <- test$sum_score
#----------------------------------------------
# setwd("/share/Projects/huanhuan/project/prog/")

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.background = element_blank(), 
                                                  axis.title.y = element_text(size = 10),
                                                  axis.title.x = element_text(size = 10),
                                                  axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5),legend.title=element_blank())

mycolor <-c("#F57F17","#673AB7","#9C27B0","#E53935","#827717","#1B5E20","#006064","#01579B")
#--------------------------------------------------------------------------8
library(plotROC)
longtest <-melt_roc(test,"new_pod_total",c("FLIPI1_count","FLIPI2_count","primapi_re_n","X150b2mg_ldh","adjust1","adjust2","adjust3","adjust4"))
longtest$name <-gsub("_count","",longtest$name)
longtest$name <-gsub("primapi_re_n","PRIMA-PI",longtest$name)
longtest$name <-gsub("X150b2mg_ldh","LDH+B2mg",longtest$name)
library(Hmisc)
longtest$name <- capitalize(longtest$name)
pdf("./figure/093_cutoff1_roc8.pdf",width=7,height=6)
p <-ggplot(longtest,aes(d=D,m=M,color=name))+geom_roc() + scale_color_manual(values=mycolor)+#,color=name aes(color=mycolor)
  p_theme+xlab("FPR")+ylab("TPR") #+scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
auc  <- calc_auc(p)$AUC
names(auc) <-levels(factor(longtest[,"name"]))
auc <-sort(auc,decreasing=T)
p1 <- p+annotate("text",x=0.85,y=0.15,
  label=paste(names(auc)[1]," AUC: ",round(auc[1],3),"\n",
              names(auc)[2]," AUC: ",round(auc[2],3),"\n",
              names(auc)[3]," AUC: ",round(auc[3],3),"\n",
              names(auc)[4]," AUC: ",round(auc[4],3),"\n",
              names(auc)[5]," AUC: ",round(auc[5],3),"\n",
              names(auc)[6]," AUC: ",round(auc[6],3),"\n",
              names(auc)[7]," AUC: ",round(auc[7],3),"\n",
              names(auc)[8]," AUC: ",round(auc[8],3),"\n"),
              size=3)
p1
dev.off()
#------------------------------------------------------
#--------------------------------------------------------------------------adjust4
mycolor <-c("#E53935","#827717","#1B5E20","#006064","#01579B")
library(plotROC)
longtest <-melt_roc(test,"new_pod_total",c("FLIPI1_count","FLIPI2_count","primapi_re_n","X150b2mg_ldh","adjust4"))
longtest$name <-gsub("_count","",longtest$name)
longtest$name <-gsub("primapi_re_n","PRIMA-PI",longtest$name)
longtest$name <-gsub("X150b2mg_ldh","LDH+B2mg",longtest$name)
library(Hmisc)
longtest$name <- capitalize(longtest$name)
# setwd("/share/Projects/huanhuan/project/prog/")
pdf("./figure/093_cutoff1_roc_adjust4.pdf",width=7,height=6)
p <-ggplot(longtest,aes(d=D,m=M,color=name))+geom_roc() + scale_color_manual(values=mycolor)+#,color=name aes(color=mycolor)
  p_theme+xlab("FPR")+ylab("TPR") #+scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
auc  <- calc_auc(p)$AUC
names(auc) <-levels(factor(longtest[,"name"]))
auc <-sort(auc,decreasing=T)
p1 <- p+annotate("text",x=0.85,y=0.15,
  label=paste(names(auc)[1]," AUC: ",round(auc[1],3),"\n",
              names(auc)[2]," AUC: ",round(auc[2],3),"\n",
              names(auc)[3]," AUC: ",round(auc[3],3),"\n",
              names(auc)[4]," AUC: ",round(auc[4],3),"\n",
              names(auc)[5]," AUC: ",round(auc[5],3),"\n"),
              size=3)
p1
dev.off()
#--------------------------------------------------------------------------adjust2
library(plotROC)
longtest <-melt_roc(test,"new_pod_total",c("FLIPI1_count","FLIPI2_count","primapi_re_n","X150b2mg_ldh","adjust2"))
longtest$name <-gsub("_count","",longtest$name)
longtest$name <-gsub("primapi_re_n","PRIMA-PI",longtest$name)
longtest$name <-gsub("X150b2mg_ldh","LDH+B2mg",longtest$name)
library(Hmisc)
longtest$name <- capitalize(longtest$name)
# setwd("/share/Projects/huanhuan/project/prog/")
pdf("./figure/093_cutoff1_roc_adjust2.pdf",width=7,height=6)
p <-ggplot(longtest,aes(d=D,m=M,color=name))+geom_roc() + scale_color_manual(values=mycolor)+#,color=name aes(color=mycolor)
  p_theme+xlab("FPR")+ylab("TPR") #+scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
auc  <- calc_auc(p)$AUC
names(auc) <-levels(factor(longtest[,"name"]))
auc <-sort(auc,decreasing=T)
p1 <- p+annotate("text",x=0.85,y=0.15,
  label=paste(names(auc)[1]," AUC: ",round(auc[1],3),"\n",
              names(auc)[2]," AUC: ",round(auc[2],3),"\n",
              names(auc)[3]," AUC: ",round(auc[3],3),"\n",
              names(auc)[4]," AUC: ",round(auc[4],3),"\n",
              names(auc)[5]," AUC: ",round(auc[5],3),"\n"),
              size=3)
p1
dev.off()
#--------------------------------------------------------------------------adjust3
library(plotROC)
longtest <-melt_roc(test,"new_pod_total",c("FLIPI1_count","FLIPI2_count","primapi_re_n","X150b2mg_ldh","adjust3"))
longtest$name <-gsub("_count","",longtest$name)
longtest$name <-gsub("primapi_re_n","PRIMA-PI",longtest$name)
longtest$name <-gsub("X150b2mg_ldh","LDH+B2mg",longtest$name)
library(Hmisc)
longtest$name <- capitalize(longtest$name)

pdf("./figure/093_cutoff1_roc_adjust3.pdf",width=7,height=6)
p <-ggplot(longtest,aes(d=D,m=M,color=name))+geom_roc() + scale_color_manual(values=mycolor)+#,color=name aes(color=mycolor)
  p_theme+xlab("FPR")+ylab("TPR") #+scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
auc  <- calc_auc(p)$AUC
names(auc) <-levels(factor(longtest[,"name"]))
auc <-sort(auc,decreasing=T)
p1 <- p+annotate("text",x=0.85,y=0.15,
  label=paste(names(auc)[1]," AUC: ",round(auc[1],3),"\n",
              names(auc)[2]," AUC: ",round(auc[2],3),"\n",
              names(auc)[3]," AUC: ",round(auc[3],3),"\n",
              names(auc)[4]," AUC: ",round(auc[4],3),"\n",
              names(auc)[5]," AUC: ",round(auc[5],3),"\n"),
              size=3)
p1
dev.off()


auc  <- calc_auc(p)$AUC
# ng <-levels(factor(longtest[,"name"]))
names(auc) <-levels(factor(longtest[,"name"]))
auc <-sort(auc,decreasing=T)

#-------------------------------------------------
#-------------------------------------------------
#------------------------------------------------------------------------FLIPI


custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}

fit <- survfit(Surv(pfs_month_new, pro_status) ~ FLIPI1_count_re, data=test)
pdf("./figure/093_pfs_FLIPI1.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#B71C1C", "#38aa34","#21618C"),
                legend.labs=c("High","Intermediate","Low"), #标签
                pval = TRUE,
                legend.title="FLIPI1",
                title="Progression Free Survival",
                xlab = " Time (Months)",
                xlim=c(0,110),
                ggtheme = custom_theme()
                )
print(p1$plot) 
dev.off()



fit <- survfit(Surv(os_month_new, dead) ~ FLIPI1_count_re, data=test)
pdf("./figure/093_os_FLIPI1.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#B71C1C", "#38aa34","#21618C"),
                legend.labs=c("High","Intermediate","Low"), #标签
                pval = TRUE,
                legend.title="FLIPI1",
                title="Overall survival",
                xlab = " Time (Months)",
                xlim=c(0,110),
                ggtheme = custom_theme()
)
print(p1$plot) 
dev.off()
#--------------------------------------------------------------------------------2
fit <- survfit(Surv(pfs_month_new, pro_status) ~ FLIPI2_count_re, data=test)
pdf("./figure/093_pfs_FLIPI2.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#B71C1C", "#38aa34","#21618C"),
                legend.labs=c("High","Intermediate","Low"), #标签
                pval = TRUE,
                legend.title="FLIPI2",
                title="Progression Free Survival",
                xlab = " Time (Months)",
                xlim=c(0,110),
                ggtheme = custom_theme()
                )
print(p1$plot) 
dev.off()

fit <- survfit(Surv(os_month_new, dead) ~ FLIPI2_count_re, data=test)
pdf("./figure/093_os_FLIPI2.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#B71C1C", "#38aa34","#21618C"),
                legend.labs=c("High","Intermediate","Low"), #标签
                pval = TRUE,
                legend.title="FLIPI1",
                title="Overall survival",
                xlab = " Time (Months)",
                xlim=c(0,110),
                ggtheme = custom_theme()
                )
print(p1$plot) 
dev.off()
#-------------------------------------------------prima_pi
fit <- survfit(Surv(pfs_month_new, pro_status) ~ primapi_re, data=test)
pdf("./figure/093_pfs_PRIMAPI.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#B71C1C", "#38aa34","#21618C"),
                legend.labs=c("High","Intermediate","Low"), #标签
                pval = TRUE,
                legend.title="PRIMA-PI",
                title="Progression Free Survival",
                xlab = " Time (Months)",
                xlim=c(0,110),
                ggtheme = custom_theme()
                )
print(p1$plot) 
dev.off()

fit <- survfit(Surv(os_month_new, dead) ~ primapi_re, data=test)
pdf("./figure/093_os_PRIMAPI.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#B71C1C", "#38aa34","#21618C"),
                legend.labs=c("High","Intermediate","Low"), #标签
                pval = TRUE,
                legend.title="PRIMA-PI",
                title="Overall survival",
                xlab = " Time (Months)",
                xlim=c(0,110),
                ggtheme = custom_theme()
                )
print(p1$plot) 
dev.off()
#--------------------------------------------------------------------------LDH + β2M
fit <- survfit(Surv(pfs_month_new, pro_status) ~ X150b2mg_ldh, data=test)
pdf("./figure/093_pfs_X150b2mg_ldh.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#B71C1C", "#38aa34","#21618C"),
                legend.labs=c("High","Intermediate","Low"), #标签
                pval = TRUE,
                legend.title="LDH + B2M",
                title="Progression Free Survival",
                xlab = " Time (Months)",
                xlim=c(0,110),
                ggtheme = custom_theme()
                )
print(p1$plot) 
dev.off()


fit <- survfit(Surv(os_month_new, dead) ~ X150b2mg_ldh, data=test)
pdf("./figure/093_os_X150b2mg_ldh.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#B71C1C", "#38aa34","#21618C"),
                legend.labs=c("High","Intermediate","Low"), #标签
                pval = TRUE,
                legend.title="LDH + B2M",
                title="Overall survival",
                main ="Overall survival",
                xlab = " Time (Months)",
                xlim=c(0,110),
                ggtheme = custom_theme()
                )
print(p1$plot) 
dev.off()

#--------------------------------------------------------------------------------------------adjust2
cutoff_y=12
test$adjust2_01 <-NA
test$adjust2_01[test$adjust2<=cutoff_y ]=0
test$adjust2_01[test$adjust2>cutoff_y]=1
library("survival")
library("survminer")

fit <- survfit(Surv(pfs_month_new, pro_status) ~ adjust2_01, data=test)
pdf("./figure/093_pod24_pfs_cutoff1_adjust2.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C","#B71C1C"),
                legend.labs=c("Low","High"), #标签
                pval = TRUE,
                legend.title="POD24",
                title="Progression Free Survival",
                xlab = " Time (Months)",
                xlim=c(0,110),
                ggtheme = custom_theme()
                )
print(p1$plot) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ adjust2_01, data=test)
pdf("./figure/093_pod24_Os_cutoff1_adjust2.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C","#B71C1C"),
                legend.labs=c("Low","High"), #标签
                pval = TRUE,
                legend.title="POD24",
                title="Overall survival",
                xlab = " Time (Months)",
                xlim=c(0,110),
                ggtheme = custom_theme()
                )
print(p1$plot) 
dev.off()

#--------------------------------------------------------------------------------------------adjust3
cutoff_y=8
test$adjust3_01 <-NA
test$adjust3_01[test$adjust3<=cutoff_y ]=0
test$adjust3_01[test$adjust3>cutoff_y]=1
library("survival")
library("survminer")

fit <- survfit(Surv(pfs_month_new, pro_status) ~ adjust3_01, data=test)
pdf("./figure/093_pod24_pfs_cutoff1_adjust3.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C","#B71C1C"),
                legend.labs=c("Low","High"), #标签
                pval = TRUE,
                legend.title="POD24",
                title="Progression Free Survival",
                xlab = " Time (Months)",
                xlim=c(0,110),
                ggtheme = custom_theme()
                )
print(p1$plot) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ adjust3_01, data=test)
pdf("./figure/093_pod24_Os_cutoff1_adjust3.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C","#B71C1C"),
                legend.labs=c("Low","High"), #标签
                pval = TRUE,
                legend.title="POD24",
                title="Overall survival",
                xlab = " Time (Months)",
                xlim=c(0,110),
                ggtheme = custom_theme()
                )
print(p1$plot) 
dev.off()

#--------------------------------------------------------------------------------------------adjust4
cutoff_y=7
test$adjust4_01 <-NA
test$adjust4_01[test$adjust4<=cutoff_y ]=0
test$adjust4_01[test$adjust4>cutoff_y]=1
library("survival")
library("survminer")

fit <- survfit(Surv(pfs_month_new, pro_status) ~ adjust4_01, data=test)
pdf("./figure/093_pod24_pfs_cutoff1_adjust4.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C","#B71C1C"),
                legend.labs=c("Low","High"), #标签
                pval = TRUE,
                legend.title="POD24",
                title="Progression Free Survival",
                xlab = " Time (Months)",
                xlim=c(0,110),
                ggtheme = custom_theme()
                )
print(p1$plot) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ adjust4_01, data=test)
pdf("./figure/093_pod24_Os_cutoff1_adjust4.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C","#B71C1C"),
                legend.labs=c("Low","High"), #标签
                pval = TRUE,
                legend.title="POD24",
                title="Overall survival",
                xlab = " Time (Months)",
                xlim=c(0,110),
                ggtheme = custom_theme()
                )
print(p1$plot) 
dev.off()