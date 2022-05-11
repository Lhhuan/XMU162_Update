setwd("/home/huanhuan/project/prog/output/")
library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)
library(cowplot)
library(ggpubr)
library(Hmisc)
setwd("/home/huanhuan/project/prog/output/")
dat <- read.csv("01_add_age_raw_pfs_os_filter_grade_mergrPod_3A.txt",sep="\t")
load("shap_weight_overlap_feature.Rdata")
dat$pod_total3 <-dat$pod_total
dat$pod_total3[grep("1|2",dat$pod_total3)] <- 0
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.background = element_blank(), 
                                                  axis.title.y = element_text(size = 10),
                                                  axis.title.x = element_text(size = 10),
                                                  axis.line = element_line(colour = "black"),legend.position ="none",plot.title = element_text(hjust = 0.5))

a <-table(dat$pod_total)%>%as.data.frame()
pod_total_new  <-table(dat$new_pod_total)%>%as.data.frame()
colnames(pod_total_new) <-c("new_pod_total","number")
pod_total_new$ratio <-pod_total_new$number/sum(pod_total_new$number)


#----------------------------------------------------------
pod_total3<-table(dat$pod_total3)%>%as.data.frame()
colnames(pod_total3) <-c("pod_total3","number")
pod_total3$ratio <-pod_total3$number/sum(pod_total3$number)

#-------------

overlap_fe <- overlap_f[-c(5,6,14)]
# overlap_fe <- c("B2mg","LN_num","LDH","age_raw","Lym_Mono","HGB","Ki.67","SPD","SUVmax","Bsym","BM")
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.background = element_blank(), 
                                                  axis.title.y = element_text(size = 12),
                                                  axis.title.x = element_text(size = 12),
                                                  axis.line = element_line(colour = "black"),legend.position ="none",plot.title = element_text(hjust = 0.5))

plot_feature_dis <-function(i){
  data <-dat[,overlap_fe[i]]%>%as.data.frame()
  colnames(data) <-"a"
  features <-overlap_fe[i]
  features<-gsub("_"," ",features)
  features <-gsub("Lym Mono","Lym/Mono",features)
  # pdf("test.pdf")
  p1 <-ggplot(data,aes(x=a))+ geom_density(fill="lightblue")+
    p_theme+labs(y="Density",x=features)
  # print(p1)
  # dev.off()
}
plist = lapply(1:9,plot_feature_dis)

a <-table(dat$Bsym)%>%as.data.frame()
plist[[10]] <-ggplot(data = a, mapping = aes(x=Var1,y=Freq)) + geom_bar(stat = 'identity', fill = "lightblue", width=0.5)+
  p_theme+labs(y="Number of samples",x="Bsym") +scale_y_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(color="black",size=8),
    axis.text.x = element_text(color="black",size=8))

a <-table(dat$BM)%>%as.data.frame()
plist[[11]] <-ggplot(data = a, mapping = aes(x=Var1,y=Freq)) + geom_bar(stat = 'identity', fill = "lightblue", width=0.5)+
  p_theme+labs(y="Number of samples",x="BM") +scale_y_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(color="black",size=8),
    axis.text.x = element_text(color="black",size=8))

pdf("./figure/066_density_plot_of_features.pdf",width=12, height=9)
# CombinePlots(plist,ncol=3,nrow=2)
p2<-gridExtra::marrangeGrob(plist,nrow=3,ncol=4,top="Distribution of features")
print(p2)
dev.off()

weight <- read.csv("067_not_fill_weight_overlap_feature_importance_weight.txt",sep="\t")
# weight_a <-weight[which(weight$Feature %in% overlap_fe),]
weight_a <-weight


weight_a$Ratio <- weight_a$Feature_importance/sum(weight_a$Feature_importance)
weight_a$Ratio_adjust1 <- round(weight_a$Ratio *100,2)
weight_a$Ratio_adjust2 <- round(weight_a$Feature_importance /5)
# weight_a$Ratio_adjust2[weight_a$Feature_importance==25] <-3
write.table(weight_a,"069_feature_ratio.txt",col.names=T,row.names=F,quote=F,sep="\t")

c_cutoff <-function(i){
  data <-dat[,overlap_fe[i]]%>%as.data.frame()
  data1 <-data[!is.na(data),]%>%data.frame()
  data1 <-data1[order(data1$.),]%>%data.frame()
  colnames(data1) <-"a"
  cutoff1 = data1[round(nrow(data1)*0.80434783),"a"]
  # cutoff2 = data1[round(nrow(data1)*0.8695652),"a"]
  rs <-data.frame(feature=overlap_fe[i],cutoff1=cutoff1)
}

ff <-lapply(c(1:5,7:9),c_cutoff)
cutoff <-do.call(rbind,ff)

c_cutoff_d <-function(i){
  data <-dat[,overlap_fe[i]]%>%as.data.frame()
  data1 <-data[!is.na(data),]%>%data.frame()
  data1 <-data1[order(-data1$.),]%>%data.frame()
  colnames(data1) <-"a"
  cutoff1 = data1[round(nrow(data1)*0.80434783),"a"]
  # cutoff2 = data1[round(nrow(data1)*0.8695652),"a"]
  rs <-data.frame(feature=overlap_fe[i],cutoff1=cutoff1)
}
ff1 <-lapply(6,c_cutoff_d)
cutoff <-bind_rows(cutoff,ff1)
cutoff <-bind_rows(cutoff, data.frame(feature="Bsym",cutoff1=0))
cutoff <-bind_rows(cutoff, data.frame(feature="BM",cutoff1=0))
write.table(cutoff,"066_cutoff_new_podtotal.txt",col.names=T,row.names=F,quote=F,sep="\t")
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


dat$SUVmax_c <-dat$SUVmax
dat$SUVmax_c[dat$SUVmax_c <=2]=0
dat$SUVmax_c[dat$SUVmax_c >2 ]=1

dat$SPD_c <-dat$SPD
dat$SPD_c[dat$SPD_c <=20]=0
dat$SPD_c[dat$SPD_c >20]=1

dat$Ki.67_c <-dat$Ki.67
dat$Ki.67_c[dat$Ki.67_c <=70.00]=0
dat$Ki.67_c[dat$Ki.67_c >70.00 ]=1

dat$HGB_c <-dat$HGB
dat$HGB_c[dat$HGB_c <=120]=1
dat$HGB_c[dat$HGB_c >120 ]=0
#-------------------------------------------------------------------------------------------------------------
#weight_a
dat$B2mg_c_f <- dat$B2mg_c *0.09905660
dat$LN_num_c_f <-dat$LN_num_c *0.13207547
dat$LDH_c_f <- dat$LDH_c *0.24528302
dat$age_raw_c_f <- dat$age_raw_c *0.06603774
dat$Lym_Mono_c_f <-dat$Lym_Mono_c *0.11792453
dat$HGB_c_f<- dat$HGB_c *0.06132075
dat$Ki.67_c_f<- dat$Ki.67_c *0.02830189
dat$SPD_c_f<- dat$SPD_c *0.11320755
dat$SUVmax_c_f<- dat$SUVmax_c *0.01886792
dat$Bsym_c_f <- dat$Bsym*0.08962264
dat$BM_c_f<-dat$BM*0.02830189
# dat$Bsym_c_f <- dat$Bsym
# dat$BM_c_f<-dat$BM
#-------------
dat_s<-dat[,c("B2mg_c_f","LN_num_c_f","LDH_c_f","age_raw_c_f","Lym_Mono_c_f","HGB_c_f","Ki.67_c_f","SPD_c_f","SUVmax_c_f","Bsym_c_f","BM_c_f")]


dat$predict_Podtotal_new <- base::rowSums(dat_s, na.rm = TRUE)
dat$Predict_Pod_total_f <-dat$predict_Podtotal_new
dat1 <-dat[order(dat$Predict_Pod_total_f),]
dat1[1:round(nrow(dat)*0.80434783),"Predict_Pod_total_f"]=0
a=round(nrow(dat)*0.80434783) +1
dat1[a:nrow(dat),"Predict_Pod_total_f"]=1
save(dat,file="069_count_predict_pod_total_new.Rdata")

library("survival")
library("survminer")

fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/069_pre_survival_3a_pod_total_new_sum_pfs_new_cutoff.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)")
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/069_pre_survival_3a_pod_total_new_sum_Os_new_cutoff.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = "Time (Months)",
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
dat$Predict_Pod_total_f <-dat$predict_Podtotal_new/100
dat1 <-dat[order(dat$Predict_Pod_total_f),]
dat1[1:round(nrow(dat)*0.80434783),"Predict_Pod_total_f"]=0
a=round(nrow(dat)*0.80434783) +1
dat1[a:nrow(dat),"Predict_Pod_total_f"]=1


library("survival")
library("survminer")

fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/069_pre_survival_3a_pod_total_new_sum_pfs_adjust1_new_cutoff.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)")
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/069_pre_survival_3a_pod_total_new_sum_Os_adjust1_new_cutoff.pdf")
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
# dat$Predict_Pod_total_f <-dat$predict_Podtotal_new/43
dat1 <-dat[order(dat$Predict_Pod_total_f),]
dat1[1:round(nrow(dat)*0.80434783),"Predict_Pod_total_f"]=0
a=round(nrow(dat)*0.80434783) +1
dat1[a:nrow(dat),"Predict_Pod_total_f"]=1

fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/069_pre_survival_3a_pod_total_new_sum_pfs_adjust2_new_cutoff.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)")
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/069_pre_survival_3a_pod_total_new_sum_Os_adjust2_new_cutoff.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()
#---------------------------------------------------------14
dat_s<-dat[,c("B2mg_c_f","LN_num_c_f","LDH_c_f","age_raw_c_f","Lym_Mono_c_f","HGB_c_f","Ki.67_c_f","SPD_c_f","SUVmax_c_f","Bsym_c_f","BM_c_f")]
dat$predict_Podtotal_new <- base::rowSums(dat_s, na.rm = TRUE)
dat$Predict_Pod_total_f <-dat$predict_Podtotal_new
dat1 <-dat[order(dat$Predict_Pod_total_f),]
dat1$Predict_Pod_total_f[dat1$Predict_Pod_total_f<22 ]=0
dat1$Predict_Pod_total_f[dat1$Predict_Pod_total_f>=22 ]=1

fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/069_pre_survival_3a_pod_total_new_sum_pfs_adjust2_new_cutoff_ycutoff_22.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)")
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/069_pre_survival_3a_pod_total_new_sum_Os_adjust2_new_cutoff_ycutoff_22.pdf")
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
# dat$Predict_Pod_total_f <-dat$predict_Podtotal_new/33
dat$Predict_Pod_total_f <-dat$predict_Podtotal_new
dat1 <-dat[order(dat$Predict_Pod_total_f),]
dat1[1:round(nrow(dat)*0.80434783),"Predict_Pod_total_f"]=0
a=round(nrow(dat)*0.80434783) +1
dat1[a:nrow(dat),"Predict_Pod_total_f"]=1 
fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/069_pre_survival_3a_pod_total_new_sum_pfs_adjust3_new_cutoff.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)")
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/069_pre_survival_3a_pod_total_new_sum_Os_adjust3_new_cutoff.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()
#------------------------------------------------------------------------------------------------------------------------------------------------------y cutoff12
dat1 <-dat[order(dat$Predict_Pod_total_f),]
dat1$Predict_Pod_total_f[dat1$Predict_Pod_total_f<8 ]=0
dat1$Predict_Pod_total_f[dat1$Predict_Pod_total_f>=8 &dat1$Predict_Pod_total_f<16 ]=3

dat1$Predict_Pod_total_f[dat1$Predict_Pod_total_f>=16 ]=4
fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/069_pre_survival_3a_pod_total_new_sum_pfs_adjust3_new_cutoff_Y_cutoff12.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)",
                xlim=c(0,120)))
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/069_pre_survival_3a_pod_total_new_sum_Os_adjust3_new_cutoff_Y_cutoff12.pdf")
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
# dat$Predict_Pod_total_f <-dat$predict_Podtotal_new/33
dat$Predict_Pod_total_f <-dat$predict_Podtotal_new
dat1 <-dat[order(dat$Predict_Pod_total_f),]
dat1[1:round(nrow(dat)*0.80434783),"Predict_Pod_total_f"]=0
a=round(nrow(dat)*0.80434783) +1
dat1[a:nrow(dat),"Predict_Pod_total_f"]=1 
fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/069_pre_survival_3a_pod_total_new_sum_pfs_adjust4_new_cutoff.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)")
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/069_pre_survival_3a_pod_total_new_sum_Os_adjust4_new_cutoff.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()
#------------------------------------------------------------------------------------------------------------------------------------------------------y cutoff11
dat1 <-dat[order(dat$Predict_Pod_total_f),]
dat1$Predict_Pod_total_f[dat1$Predict_Pod_total_f<11 ]=0
# dat1$Predict_Pod_total_f[dat1$Predict_Pod_total_f>=10 &dat1$Predict_Pod_total_f<20 ]=3

dat1$Predict_Pod_total_f[dat1$Predict_Pod_total_f>=11 ]=1
fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/069_pre_survival_3a_pod_total_new_sum_pfs_adjust4_new_cutoff_Y_cutoff11.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)")
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/069_pre_survival_3a_pod_total_new_sum_Os_adjust4_new_cutoff_Y_cutoff11.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()




#-------------------------------------------------------------------------------------------------------------------------------------------pod_total3，features cutoff2
pod_total3  <-table(dat$pod_total3)%>%as.data.frame()
colnames(pod_total3) <-c("pod_total3","number")
pod_total3$ratio <-pod_total3$number/sum(pod_total3$number)
#-------
c_cutoff <-function(i){
  data <-dat[,overlap_fe[i]]%>%as.data.frame()
  data1 <-data[!is.na(data),]%>%data.frame()
  data1 <-data1[order(data1$.),]%>%data.frame()
  colnames(data1) <-"a"
  cutoff1 = data1[round(nrow(data1)*0.80434783),"a"]
  cutoff2 = data1[round(nrow(data1)*0.8695652),"a"]
  rs <-data.frame(feature=overlap_fe[i],cutoff1=cutoff1,cutoff2=cutoff2)
}

ff <-lapply(c(1:5,7:9),c_cutoff)
cutoff <-do.call(rbind,ff)
#-------
c_cutoff_d <-function(i){
  data <-dat[,overlap_fe[i]]%>%as.data.frame()
  data1 <-data[!is.na(data),]%>%data.frame()
  data1 <-data1[order(-data1$.),]%>%data.frame()
  colnames(data1) <-"a"
  cutoff1 = data1[round(nrow(data1)*0.80434783),"a"]
  cutoff2 = data1[round(nrow(data1)*0.8695652),"a"]
  rs <-data.frame(feature=overlap_fe[i],cutoff1=cutoff1,cutoff2=cutoff2)
}
ff1 <-lapply(6,c_cutoff_d)
cutoff <-bind_rows(cutoff,ff1)
cutoff <-bind_rows(cutoff, data.frame(feature="Bsym",cutoff1=0,cutoff2=1))
cutoff <-bind_rows(cutoff, data.frame(feature="BM",cutoff1=0,cutoff2=1))
write.table(cutoff,"066_cutoff_podtotal3.txt",col.names=T,row.names=F,quote=F,sep="\t")

#-------------------------------------------
dat$B2mg_c <-dat$B2mg
dat$B2mg_c[dat$B2mg_c <=3.4]=0
dat$B2mg_c[dat$B2mg_c >3.4 & dat$B2mg_c <=3.8 ]=1
dat$B2mg_c[dat$B2mg_c >3.8 ]=2

dat$LN_num_c <-dat$LN_num
dat$LN_num_c[dat$LN_num_c <=6]=0
dat$LN_num_c[dat$LN_num_c >6 & dat$LN_num_c<=11]=1
dat$LN_num_c[dat$LN_num_c >11]=2

dat$LDH_c <-dat$LDH
dat$LDH_c[dat$LDH_c <=270]=0
dat$LDH_c[dat$LDH_c >270 & dat$LDH_c <=326]=1
dat$LDH_c[dat$LDH_c >326]=2

dat$age_raw_c <-dat$age_raw
dat$age_raw_c[dat$age_raw_c <=60]=0
dat$age_raw_c[dat$age_raw_c >60 & dat$age_raw_c <=66]=1
dat$age_raw_c[dat$age_raw_c >66 ]=2

dat$Lym_Mono_c <-dat$Lym_Mono
dat$Lym_Mono_c[dat$Lym_Mono_c <=6]=0
dat$Lym_Mono_c[dat$Lym_Mono_c >6 & dat$Lym_Mono_c <=6.5 ]=1
dat$Lym_Mono_c[dat$Lym_Mono_c >6.5 ]=2

dat$SUVmax_c <-dat$SUVmax
dat$SUVmax_c[dat$SUVmax_c <=2]=0
dat$SUVmax_c[dat$SUVmax_c >2 & dat$SUVmax_c<=20]=1
dat$SUVmax_c[dat$SUVmax_c >20 ]=2

dat$SPD_c <-dat$SPD
dat$SPD_c[dat$SPD_c <=0]=0
dat$SPD_c[dat$SPD_c >0 & dat$SPD_c <=33 ]=1
dat$SPD_c[dat$SPD_c >33]=2

dat$Ki.67_c <-dat$Ki.67
dat$Ki.67_c[dat$Ki.67_c <=20.00]=0
dat$Ki.67_c[dat$Ki.67_c >20.00 & dat$Ki.67_c <=70.00]=1
dat$Ki.67_c[dat$Ki.67_c >70.00 ]=2

dat$HGB_c <-dat$HGB
dat$HGB_c[dat$HGB_c >120]=0
dat$HGB_c[dat$HGB_c >= 111 & dat$HGB_c <120]=1
dat$HGB_c[dat$HGB_c =<111]=2


dat$B2mg_c_f <- dat$B2mg_c*3
dat$LN_num_c_f <-dat$LN_num_c *3
dat$LDH_c_f <- dat$LDH_c * 3
dat$age_raw_c_f <- dat$age_raw_c*3
dat$Lym_Mono_c_f <-dat$Lym_Mono_c*2
dat$HGB_c_f<- dat$HGB_c *2
dat$Ki.67_c_f<- dat$Ki.67_c *2
dat$SPD_c_f<- dat$SPD_c *2
dat$SUVmax_c_f<- dat$SUVmax_c *2
dat$Bsym_c_f <- dat$Bsym*1
dat$BM_c_f<-dat$BM*1
dat_s<-dat[,c("B2mg_c_f","LN_num_c_f","LDH_c_f","age_raw_c_f","Lym_Mono_c_f","HGB_c_f","Ki.67_c_f","SPD_c_f","SUVmax_c_f","Bsym_c_f","BM_c_f")]
dat$predict_Podtotal_new <- base::rowSums(dat_s, na.rm = TRUE)
# dat$Predict_Pod_total_f <-dat$predict_Podtotal_new/33
dat$Predict_Pod_total_f <-dat$predict_Podtotal_new
dat1 <-dat[order(dat$Predict_Pod_total_f),]
dat1[1:round(nrow(dat)*0.80434783),"Predict_Pod_total_f"]=0
a=round(nrow(dat)*0.80434783) +1
dat1[a:nrow(dat),"Predict_Pod_total_f"]=1 
fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/066_pre_survival_3a_pod_total_new_sum_pfs_adjust3_new_cutoff_features3.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)")
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/066_pre_survival_3a_pod_total_new_sum_Os_adjust3_new_cutoff_features3.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()

#--------------------------------------------------------
dat1 <-dat[order(dat$Predict_Pod_total_f),]
dat1$Predict_Pod_total_f[dat1$Predict_Pod_total_f<12 ]=0
dat1$Predict_Pod_total_f[dat1$Predict_Pod_total_f>=12 ]=1

fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/066_pre_survival_3a_pod_total_new_sum_pfs_adjust3_new_cutoff_features3_Y_cutoff12.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)")
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/066_pre_survival_3a_pod_total_new_sum_Os_adjust3_new_cutoff_features3_Y_cutoff12.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()


#------------------------------------------------------------------------------------------------------------------------------------------------------y cutoff12,24
dat1 <-dat[order(dat$Predict_Pod_total_f),]
dat1$Predict_Pod_total_f[dat1$Predict_Pod_total_f<12 ]=0
dat1$Predict_Pod_total_f[dat1$Predict_Pod_total_f>=12 &dat1$Predict_Pod_total_f <24 ]=3
dat1$Predict_Pod_total_f[dat1$Predict_Pod_total_f>=24 ]=4

fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/066_pre_survival_3a_pod_total_new_sum_pfs_adjust3_new_cutoff_Y_features3_cutoff12_24.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Pfs",
                 xlab = " Time (Months)")
print(p1) 
dev.off()
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=dat1)
pdf("./figure/066_pre_survival_3a_pod_total_new_sum_Os_adjust3_new_cutoff_Y_features3_cutoff12_24.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()
