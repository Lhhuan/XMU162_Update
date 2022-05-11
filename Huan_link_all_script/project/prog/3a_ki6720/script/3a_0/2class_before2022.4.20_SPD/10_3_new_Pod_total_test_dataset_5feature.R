library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)
library(cowplot)
library(ggpubr)
library(Hmisc)
library("survival")
library("survminer")

setwd("/home/huanhuan/project/prog/script/3a_0/output/")
load("07_add_model_index_age_raw_pfs_os_filter_grade_mergrPod_3A_0.Rdata")
load("09_test_train_dataset.Rdata")
dat$primapi_re_n <- NA
dat$primapi_re_n[dat$primapi_re =="Low"]=0
dat$primapi_re_n[dat$primapi_re =="Intermediate"]=1
dat$primapi_re_n[dat$primapi_re =="High"]=2
#------------------------------------------------------------------new_cutoff
dat$LDH_c <-dat$LDH_re0
dat$LN_num_c <-dat$LN_num_s
dat$HGB_c =dat$HGB_s

dat$ECOG_c =NA
dat$ECOG_c[dat$ECOG <=1]=0
dat$ECOG_c[dat$ECOG >1]=1
dat$ECOG_c[is.na(dat$ECOG_c)&dat$new_pod_total==1]=1
dat$ECOG_c[is.na(dat$ECOG_c)&dat$new_pod_total==0]=0
dat$B2M_c <-dat$B2MG_re0

dat$Lym_Mono_c <-NA
dat$Lym_Mono_c[dat$Lym_Mono <=10]=0
dat$Lym_Mono_c[dat$Lym_Mono >10 ]=1
dat$Lym_Mono_c[is.na(dat$Lym_Mono_c)&dat$new_pod_total==1]=1
dat$Lym_Mono_c[is.na(dat$Lym_Mono_c)&dat$new_pod_total==0]=0
dat$BM_c <-dat$BM_s

dat$SPD_c <-NA
dat$SPD_c[dat$SPD <=20]=0
dat$SPD_c[dat$SPD >20]=1
dat$SPD_c[is.na(dat$SPD_c)&dat$new_pod_total==1]=1
dat$SPD_c[is.na(dat$SPD_c)&dat$new_pod_total==0]=0

# dat$SUVmax_c <-NA
# dat$SUVmax_c[dat$SUVmax <=2]=0
# dat$SUVmax_c[dat$SUVmax >2 ]=1
# dat$SUVmax_c[is.na(dat$SUVmax_c)&dat$new_pod_total==1]=1
# dat$SUVmax_c[is.na(dat$SUVmax_c)&dat$new_pod_total==0]=0

# dat$Lym_Mono_c <-NA
# dat$Lym_Mono_c[dat$Lym_Mono <=10]=0
# dat$Lym_Mono_c[dat$Lym_Mono >10 ]=1
# dat$Lym_Mono_c[is.na(dat$Lym_Mono_c)&dat$new_pod_total==1]=1
# dat$Lym_Mono_c[is.na(dat$Lym_Mono_c)&dat$new_pod_total==0]=0



# dat$LN6_c <-dat$LN6_s

# #-------------------------------------------------------------------------------------------------------------
test <-dat[test_set_number,]
# im <-read.table("/share/Projects/huanhuan/project/prog/10_2_feature_ratio.txt",header=T,sep="\t")

# test$sum <-base::rowSums(test[,c("Lym_Mono","B2mg","SPD","LDH","LN_num","SUVmax","Ki.67","HGB","BM","Bsym")])
# base::rowSums(dat_s, na.rm = TRUE)

#-----------------------------------------------adjust2
test$LDH_c_f <- test$LDH_c*11
test$LN_num_c_f <-test$LN_num_c *11
test$HGB_c_f<- test$HGB_c*10
test$B2M_c_f <-test$B2M_c*9
test$ECOG_c_f<- test$ECOG_c *8

test$sum_score<-base::rowSums(test[,c("LDH_c_f","LN_num_c_f","HGB_c_f","B2M_c_f","ECOG_c_f")])
# test[is.na($sum_score),]
test$adjust2 <- test$sum_score
#-----------------------------------------------adjust3
test$LDH_c_f <- test$LDH_c*5
test$LN_num_c_f <-test$LN_num_c *5
test$HGB_c_f<- test$HGB_c*5
test$B2M_c_f <-test$B2M_c*4
test$ECOG_c_f<- test$ECOG_c *4

test$sum_score<-base::rowSums(test[,c("LDH_c_f","LN_num_c_f","HGB_c_f","B2M_c_f","ECOG_c_f")])
# test[is.na($sum_score),]
test$adjust3 <- test$sum_score
#-----------------------------------------------adjust4
test$LDH_c_f <- test$LDH_c*5
test$LN_num_c_f <-test$LN_num_c *5
test$HGB_c_f<- test$HGB_c*5
test$B2M_c_f <-test$B2M_c*5
test$ECOG_c_f<- test$ECOG_c *4

test$sum_score<-base::rowSums(test[,c("LDH_c_f","LN_num_c_f","HGB_c_f","B2M_c_f","ECOG_c_f")])
# test[is.na($sum_score),]
test$adjust4 <- test$sum_score
#--------------------------------------
#-----------------------------------------------adjust5
test$LDH_c_f <- test$LDH_c*6
test$LN_num_c_f <-test$LN_num_c *6
test$HGB_c_f<- test$HGB_c*5
test$B2M_c_f <-test$B2M_c*5
test$ECOG_c_f<- test$ECOG_c *4

test$sum_score<-base::rowSums(test[,c("LDH_c_f","LN_num_c_f","HGB_c_f","B2M_c_f","ECOG_c_f")])
# test[is.na($sum_score),]
test$adjust5 <- test$sum_score

#-----------------------------------------------adjust6
test$LDH_c_f <- test$LDH_c*6
test$LN_num_c_f <-test$LN_num_c *6
test$HGB_c_f<- test$HGB_c*5
test$B2M_c_f <-test$B2M_c*4
test$ECOG_c_f<- test$ECOG_c *4

test$sum_score<-base::rowSums(test[,c("LDH_c_f","LN_num_c_f","HGB_c_f","B2M_c_f","ECOG_c_f")])
# test[is.na($sum_score),]
test$adjust6 <- test$sum_score

#-----------------------------------------------adjust7
test$LDH_c_f <- test$LDH_c*2
test$LN_num_c_f <-test$LN_num_c *2
test$HGB_c_f<- test$HGB_c*2
test$B2M_c_f <-test$B2M_c*1
test$ECOG_c_f<- test$ECOG_c *1

test$sum_score<-base::rowSums(test[,c("LDH_c_f","LN_num_c_f","HGB_c_f","B2M_c_f","ECOG_c_f")])
# test[is.na($sum_score),]
test$adjust7 <- test$sum_score


p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.background = element_rect(color="black",size=1.2),
                                                  axis.title.y = element_text(size = 15),
                                                  axis.title.x = element_text(size = 15),
                                                  axis.text.y = element_text(size = 12,colour = "black"),
                                                  axis.text.x = element_text(size = 12,colour = "black"),
                                                  plot.title = element_text(hjust = 0.5),legend.title=element_blank(),legend.text=element_text(size=9))

mycolor <-c("#673AB7","#9C27B0","#E53935","#827717","#1B5E20","#006064","#01579B")
#--------------------------------------------------------------------------8
library(plotROC)
longtest <-melt_roc(test,"new_pod_total",c("FLIPI1_count","FLIPI2_count","primapi_re_n","b2mg_ldh1.5s","adjust6","adjust7","adjust4"))
longtest$name <-gsub("_count","",longtest$name)
longtest$name <-gsub("primapi_re_n","PRIMA-PI",longtest$name)
longtest$name <-gsub("b2mg_ldh1.5s","LDH+B2mg",longtest$name)
library(Hmisc)
longtest$name <- capitalize(longtest$name)
# pdf("./figure/10_3_cutoff1_roc7.pdf",width=7,height=6)
p <-ggplot(longtest,aes(d=D,m=M,color=name))+geom_roc(show.legend = TRUE,labels=TRUE) + scale_color_manual(values=mycolor)+#,color=name aes(color=mycolor)
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
              names(auc)[7]," AUC: ",round(auc[7],3),"\n"),
              size=3.5)
            
ggsave("./figure/10_3_cutoff1_roc7_new.png",p1,dpi=300,width=7,height=5.8)

#------------------------------------------------------
#--------------------------------------------------------------------------adjust5
TP <-length(test$new_pod_total[test$adjust5>=5&test$new_pod_total==1])
FP <-length(test$new_pod_total[test$adjust5>=5&test$new_pod_total==0])
FN <-length(test$new_pod_total[test$adjust5<5&test$new_pod_total==1])
TN <-length(test$new_pod_total[test$adjust5<5&test$new_pod_total==0])
TPR = TP / (TP+FN)
FPR = FP / (FP + TN)
# mycolor <-c("#E53935","#827717","#1B5E20","#006064","#01579B")
mycolor <-c("#827717","#1B5E20","#006064","#E53935","#01579B")
test$Ours <-test$adjust5

library(plotROC)
longtest <-melt_roc(test,"new_pod_total",c("FLIPI1_count","FLIPI2_count","primapi_re_n","b2mg_ldh1.5s","Ours"))
longtest$name <-gsub("_count","",longtest$name)
longtest$name <-gsub("primapi_re_n","PRIMA-PI",longtest$name)
longtest$name <-gsub("b2mg_ldh1.5s","LDH+B2mg",longtest$name)
library(Hmisc)
longtest$name <- capitalize(longtest$name)
# setwd("/share/Projects/huanhuan/project/prog/")
# pdf("./figure/10_3_cutoff1_roc_adjust3.pdf",width=7,height=6)
# p <-ggplot(longtest,aes(d=D,m=M,color=name))+geom_roc(n.cuts = 0,show.legend = TRUE,labels=TRUE) + scale_color_manual(values=mycolor)+#,color=name aes(color=mycolor)
p <-ggplot(longtest,aes(d=D,m=M,color=name))+geom_roc(n.cuts = 0,show.legend = TRUE,labels=TRUE) + scale_color_manual(values=mycolor)+#,color=name aes(color=mycolor)
 theme_bw() +p_theme+xlab("FPR")+ylab("TPR") #+scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
auc  <- calc_auc(p)$AUC
names(auc) <-levels(factor(longtest[,"name"]))
auc <-sort(auc,decreasing=T)
p1 <- p+annotate("text",x=0.85,y=0.15,
  label=paste(names(auc)[1]," AUC: ",round(auc[1],3),"\n",
              names(auc)[2]," AUC: ",round(auc[2],3),"\n",
              names(auc)[3]," AUC: ",round(auc[3],3),"\n",
              names(auc)[4]," AUC: ",round(auc[4],3),"\n",
              names(auc)[5]," AUC: ",round(auc[5],3),"\n"),
              size=3.5)
p1 <-p1 +annotate("text",x=0.43,y=0.82,label= "5",size=5,color="#E53935" )+geom_point(x=FPR,y=TPR ,size=4,color="#E53935" )

  p2 <-p1 +scale_x_continuous(limits= c(0, 1), breaks= seq(0,1,0.2))+ #expand = c(0, 0)
  scale_y_continuous(limits= c(0, 1), breaks= seq(0,1,0.2)) #+theme(panel.background = element_rect(color="black",size=1.2))
ggsave("./figure/10_3_roc_adjust5_NEW.png",p2,dpi=300,width=7,height=5.8)
pdf("aaa.pdf")
print(p2)
dev.off()

#------------------------------------------------------------------------FLIPI


custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}


cutoff_y=2
test$FLIPI1_01 <-NA
test$FLIPI1_01[test$FLIPI1_count<cutoff_y ]=0
test$FLIPI1_01[test$FLIPI1_count>=cutoff_y]=1
fit <- survfit(Surv(pfs_month_new, pro_status) ~ FLIPI1_01, data=test)
# pdf("./figure/093_pfs_FLIPI1.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C","#B71C1C"),
                legend.labs=c("Low risk","High risk"), #标签
                pval = TRUE,
                legend.title="POD24",
                title="PFS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
# print(p1$plot) 

ggsave("./figure/10_3_pfs_FLIPI1_test.png",p1$plot,dpi=300,height=5.2,width=5)
# dev.off()



fit <- survfit(Surv(os_month_new, dead) ~ FLIPI1_01, data=test)
# pdf("./figure/093_os_FLIPI1.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C","#B71C1C"),
                legend.labs=c("Low risk","High risk"), #标签
                pval = TRUE,
                legend.title="POD24",
                title="OS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
# print(p1$plot) 

ggsave("./figure/10_3_OS_FLIPI1_test.png",p1$plot,dpi=300,height=5.2,width=5)
#--------------------------------------------------------------------------------2
cutoff_y=2
test$FLIPI2_01 <-NA
test$FLIPI2_01[test$FLIPI2_count<cutoff_y ]=0
test$FLIPI2_01[test$FLIPI2_count>=cutoff_y]=1

fit <- survfit(Surv(pfs_month_new, pro_status) ~ FLIPI2_01, data=test)
# pdf("./figure/093_pfs_FLIPI1.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C","#B71C1C"),
                legend.labs=c("Low risk","High risk"), #标签
                pval = TRUE,
                legend.title="POD24",
                title="PFS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
# print(p1$plot) 

ggsave("./figure/10_3_pfs_FLIPI2_test.png",p1$plot,dpi=300,height=5.2,width=5)

fit <- survfit(Surv(os_month_new, dead) ~ FLIPI2_01, data=test)
# pdf("./figure/093_os_FLIPI1.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C","#B71C1C"),
                legend.labs=c("Low risk","High risk"), #标签
                pval = TRUE,
                legend.title="POD24",
                title="OS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
# print(p1$plot) 

ggsave("./figure/10_3_OS_FLIPI2_test.png",p1$plot,dpi=300,height=5.2,width=5)
#-------------------------------------------------prima_pi

cutoff_y=1
test$primapi_01 <-NA
test$primapi_01[test$primapi_re_n<cutoff_y ]=0
test$primapi_01[test$primapi_re_n>=cutoff_y]=1

fit <- survfit(Surv(pfs_month_new, pro_status) ~ primapi_01, data=test)
# pdf("./figure/093_pfs_PRIMAPI.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C","#B71C1C"),
                legend.labs=c("Low risk","High risk"), #标签
                pval = TRUE,
                legend.title="POD24",
                title="PFS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
# print(p1$plot) 

ggsave("./figure/10_3_pfs_primapi_test.png",p1$plot,dpi=300,height=5.2,width=5)

fit <- survfit(Surv(os_month_new, dead) ~ primapi_01, data=test)
# pdf("./figure/093_os_PRIMAPI.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C","#B71C1C"),
                legend.labs=c("Low risk","High risk"), #标签
                pval = TRUE,
                legend.title="POD24",
                title="OS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
# print(p1$plot) 

ggsave("./figure/10_3_OS_primapi_test.png",p1$plot,dpi=300,height=5.2,width=5)
#--------------------------------------------------------------------------LDH + β2M
cutoff_y=1
test$b2mldh_01 <-NA
test$b2mldh_01[test$b2mg_ldh1.5s<cutoff_y ]=0
test$b2mldh_01[test$b2mg_ldh1.5s>=cutoff_y]=1
fit <- survfit(Surv(pfs_month_new, pro_status) ~ b2mldh_01, data=test)
# pdf("./figure/093_pfs_b2mg_ldh1.5s.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C","#B71C1C"),
                legend.labs=c("Low risk","High risk"), #标签
                pval = TRUE,
                legend.title="POD24",
                title="PFS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
# print(p1$plot) 

ggsave("./figure/10_3_pfs_b2mg_ldh1.5s_test.png",p1$plot,dpi=300,height=5.2,width=5)


fit <- survfit(Surv(os_month_new, dead) ~ b2mldh_01, data=test)
# pdf("./figure/093_os_b2mg_ldh1.5s.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C","#B71C1C"),
                legend.labs=c("Low risk","High risk"), #标签
                pval = TRUE,
                legend.title="POD24",
                title="OS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
# print(p1$plot) 

ggsave("./figure/10_3_OS_b2mg_ldh1.5s_test.png",p1$plot,dpi=300,height=5.2,width=5)


#--------------------------------------------------------------------------------------------adjust7
# cutoff_y=4
# test$adjust7_01 <-NA
# test$adjust7_01[test$adjust7<cutoff_y ]=0
# test$adjust7_01[test$adjust7>=cutoff_y]=1
# library("survival")
# library("survminer")

# fit <- survfit(Surv(pfs_month_new, pro_status) ~ adjust7_01, data=test)
# # pdf("./figure/10_3_pod24_pfs_cutoff1_adjust3.pdf",height=5.2,width=5)
# p1 <- ggsurvplot(fit,
#                 palette=c("#21618C","#B71C1C"),
#                 legend.labs=c("Low","High"), #标签
#                 pval = TRUE,
#                 legend.title="POD24",
#                 title="Progression Free Survival",
#                 xlab = " Time (Months)",
#                 xlim=c(0,110),
#                 ggtheme = custom_theme()
#                 )
# ggsave("./figure/10_3_pod24_pfs_cutoff1_adjust7_new.png",p1$plot,dpi=300,height=5.2,width=5)
# #------------------------------
# fit <- survfit(Surv(os_month_new, dead) ~ adjust7_01, data=test)
# # pdf("./figure/10_3_pod24_Os_cutoff1_adjust3.pdf",height=5.2,width=5)
# p1 <- ggsurvplot(fit,
#                 palette=c("#21618C","#B71C1C"),
#                 legend.labs=c("Low","High"), #标签
#                 pval = TRUE,
#                 legend.title="POD24",
#                 title="Overall survival",
#                 xlab = " Time (Months)",
#                 xlim=c(0,110),
#                 ggtheme = custom_theme()
#                 )
# ggsave("./figure/10_3_pod24_Os_cutoff1_adjust7_new.png",p1$plot,dpi=300,height=5.2,width=5)
#-----------------
#--------------------------------------------------------------------------------------------adjust5
cutoff_y=5
test$adjust5_01 <-NA
test$adjust5_01[test$adjust5<cutoff_y ]=0
test$adjust5_01[test$adjust5>=cutoff_y]=1
library("survival")
library("survminer")

fit <- survfit(Surv(pfs_month_new, pro_status) ~ adjust5_01, data=test)
# pdf("./figure/10_3_pod24_pfs_cutoff1_adjust5.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C","#B71C1C"),
                legend.labs=c("Low risk","High risk"), #标签
                pval = TRUE,
                legend.title="POD24",
                title="PFS",
                xlab = " Time (Months)",
                xlim=c(0,110),
                ggtheme = custom_theme()
                )
ggsave("./figure/10_3_pod24_pfs_cutoff1_adjust5_new.png",p1$plot,dpi=300,height=5.2,width=5)

fit <- survfit(Surv(os_month_new, dead) ~ adjust5_01, data=test)
# pdf("./figure/10_3_pod24_Os_cutoff1_adjust5.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C","#B71C1C"),
                legend.labs=c("Low risk","High risk"), #标签
                pval = TRUE,
                legend.title="POD24",
                title="OS",
                xlab = " Time (Months)",
                xlim=c(0,110),
                ggtheme = custom_theme()
                )
ggsave("./figure/10_3_pod24_Os_cutoff1_adjust5_new.png",p1$plot,dpi=300,height=5.2,width=5)


#--------------------------------------------------------------------------------------------ori


custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}

fit <- survfit(Surv(pfs_month_new, pro_status) ~ new_pod_total, data=dat)
# pdf("./figure/10_3_pod24_pfs_cutoff1_adjust3.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C","#B71C1C"),
                legend.labs=c("Low risk","High risk"), #标签
                pval = TRUE,
                legend.title="POD24",
                title="PFS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
ggsave("./figure/10_3_pod24_pfs_ori.png",p1$plot,dpi=300,height=5.2,width=5)

fit <- survfit(Surv(os_month_new, dead) ~new_pod_total, data=dat)
# pdf("./figure/10_3_pod24_Os_cutoff1_adjust3.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C","#B71C1C"),
                legend.labs=c("Low risk","High risk"), #标签
                pval = TRUE,
                legend.title="POD24",
                title="OS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
ggsave("./figure/10_3_pod24_os_ori.png",p1$plot,dpi=300,height=5.2,width=5)

#---------------------------------------------------------------------------------
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.background = element_blank(), 
                                                  axis.title.y = element_text(size = 10),
                                                  axis.title.x = element_text(size = 10),
                                                  axis.line = element_line(colour = "black"),legend.position ="none",plot.title = element_text(hjust = 0.5))


pod_total_new  <-table(dat$new_pod_total)%>%as.data.frame()
# pdf("./figure/10_3_distrbution_pod_total_new.pdf",width = 5,height = 5)
p1 <-ggplot(data = pod_total_new , mapping = aes(x=Var1,y=Freq)) + geom_bar(stat = 'identity', fill = "#5eaaa8", width=0.5)+
  p_theme+labs(y="Number of samples",x="Pod total class") +scale_y_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black"))

ggsave("./figure/10_3_distrbution_pod_total_new.png",p1,dpi=300,height=4,width=4)
