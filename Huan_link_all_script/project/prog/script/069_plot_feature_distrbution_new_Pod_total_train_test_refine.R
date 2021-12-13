setwd("/home/huanhuan/project/prog/output/")
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
load("shap_weight_overlap_feature.Rdata")

pod_total0=which(dat$new_pod_total==0)
set.seed(111)
test_number0 <-sample(x=pod_total0, round(length(pod_total0)*1/3),replace = F)
pod_total1=which(dat$new_pod_total!=0)
set.seed(111)
test_number1 <-sample(x=pod_total1, round(length(pod_total1)*1/3),replace = F)

test_set_number = c(test_number0,test_number1)
train_set_number =setdiff(1:nrow(dat),test_set_number)



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
#-------------------------------
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
cutoff_y = round(train[271,"Predict_Pod_total_f"]) #nrow(train)*0.8041543

test$Predict_Pod_total_f[test$Predict_Pod_total_f<cutoff_y ]=0
test$Predict_Pod_total_f[test$Predict_Pod_total_f>=cutoff_y]=1



fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=test)
pdf("./figure/069_pre_survival_3a_pod_total_new_sum_pfs_adjust3_new_cutoff_y8.pdf")
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
pdf("./figure/069_pre_survival_3a_pod_total_new_sum_Os_adjust3_new_cutoff_y8.pdf")
p1 <- ggsurvplot(fit,
                  pval = TRUE,
                 legend.title="Pod total new",
                 title="3a: Os",
                 xlab = "Time (Months)",
                 xlim=c(0,110))
print(p1) 
dev.off()
