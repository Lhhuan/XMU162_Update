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
# set.seed(178)
# set.seed(188)
set.seed(2)
test_number0 <-sample(x=pod_total0, round(length(pod_total0)*1/3),replace = F)
pod_total1=which(dat$new_pod_total!=0)
# set.seed(112231)
# set.seed(1124)
# set.seed(178)
# set.seed(188)
set.seed(2)
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
dat$FLIPI1_c <-NA  
dat$FLIPI1_c[dat$FLIPI1_count_re=="Intermediate"|dat$FLIPI1_count_re=="Low"]=0
dat$FLIPI1_c[dat$FLIPI1_count_re=="High"]=1
#-----
dat$FLIPI2_c <-NA  
dat$FLIPI2_c[dat$FLIPI2_count_re=="Intermediate"|dat$FLIPI2_count_re=="Low"]=0
dat$FLIPI2_c[dat$FLIPI2_count_re=="High"]=1

dat$primapi_re_c <-NA  
dat$primapi_re_c[dat$primapi_re=="Intermediate"|dat$primapi_re=="Low"]=0
dat$primapi_re_c[dat$primapi_re=="High"]=1
dat$b2mg_ldh_c <-NA  
dat$b2mg_ldh_c[dat$X150b2mg_ldh <2]=0
dat$b2mg_ldh_c[dat$X150b2mg_ldh ==2]=1
#------------------------------------------------------------------------FLIPI
test=dat[test_set_number,]

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}

fit <- survfit(Surv(pfs_month_new, pro_status) ~ FLIPI1_count_re, data=test)
pdf("./figure/08_pre_survival_3a_pod_pfs_FLIPI1.pdf",height=5.2,width=5)
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
pdf("./figure/08_pre_survival_3a_os_FLIPI1.pdf",height=5.2,width=5)
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
pdf("./figure/08_pre_survival_3a_pfs_FLIPI2.pdf",height=5.2,width=5)
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
pdf("./figure/08_pre_survival_3a_os_FLIPI2.pdf",height=5.2,width=5)
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
pdf("./figure/08_pre_survival_3a_pfs_PRIMAPI.pdf",height=5.2,width=5)
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
pdf("./figure/08_pre_survival_3a_os_PRIMAPI.pdf",height=5.2,width=5)
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
pdf("./figure/08_pre_survival_3a_pfs_X150b2mg_ldh.pdf",height=5.2,width=5)
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
pdf("./figure/08_pre_survival_3a_os_X150b2mg_ldh.pdf",height=5.2,width=5)
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
adjust1 <-test[,c("new_pod_total","Predict_Pod_total_f")]

library("survival")
library("survminer")

fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=test)
pdf("./figure/08_pre_survival_3a_pod24_pfs_adjust1.pdf",height=5.2,width=5)
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
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=test)
pdf("./figure/08_pre_survival_3a_pod24_Os_adjust1.pdf",height=5.2,width=5)
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
adjust2 <-test[,c("new_pod_total","Predict_Pod_total_f")]

fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=test)
pdf("./figure/08_pre_survival_3a_pod24_pfs_adjust2.pdf",height=5.2,width=5)
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
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=test)
pdf("./figure/08_pre_survival_3a_pod24_Os_adjust2.pdf",height=5.2,width=5)
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
adjust3 <-test[,c("new_pod_total","Predict_Pod_total_f")]


fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=test)
pdf("./figure/08_pre_survival_3a_pod24_pfs_adjust3.pdf",height=5.2,width=5)
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
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=test)
pdf("./figure/08_pre_survival_3a_pod24_os_adjust3.pdf",height=5.2,width=5)
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
adjust4 <-test[,c("new_pod_total","Predict_Pod_total_f")]

fit <- survfit(Surv(pfs_month_new, pro_status) ~ Predict_Pod_total_f, data=test)
pdf("./figure/08_pre_survival_3a_pod24_pfs_adjust4.pdf",height=5.2,width=5)
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
fit <- survfit(Surv(os_month_new, dead) ~ Predict_Pod_total_f, data=test)
pdf("./figure/08_pre_survival_3a_pod24_os_adjust4.pdf",height=5.2,width=5)
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

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(pROC)
roc_list=NULL
roc_list[[1]] <-roc(test$new_pod_total,test$FLIPI1_c,direction=  "<")
roc_list[[2]] <-roc(test$new_pod_total,test$FLIPI2_c,,direction=  "<")
roc_list[[3]] <-roc(test$new_pod_total,test$primapi_re_c,direction=  "<")
roc_list[[4]] <-roc(test$new_pod_total,test$b2mg_ldh_c,direction=  "<")
roc_list[[5]] <-roc(adjust1$new_pod_total,adjust1$Predict_Pod_total_f,direction=  "<")
roc_list[[6]] <-roc(adjust2$new_pod_total,adjust2$Predict_Pod_total_f,direction=  "<")
roc_list[[7]] <-roc(adjust3$new_pod_total,adjust3$Predict_Pod_total_f,direction=  "<")
roc_list[[8]] <-roc(adjust4$new_pod_total,adjust4$Predict_Pod_total_f,direction=  "<")

# mycolor <-c("#cc561e","#beca5c","#709fb0","#ff9292","#5b6d5b","#c15050","#e84545","#9ede73","#845ec2","#28527a")
# mycolor <-c("#FFEE58","#beca5c","#709fb0","#ff9292","#5b6d5b","#c15050","#e84545","#28527a")
mycolor <-c("#FFEE58","#D4E157","#8BC34A","#4CAF50","#3F51B5","#673AB7","#9C27B0","#E53935")
pdf("./figure/08_model_ROC.pdf")

plot(roc_list[[1]],col=mycolor[1],plot=TRUE,legacy.axes=TRUE)
plot.roc(roc_list[[2]],add=TRUE,col=mycolor[2],plot=TRUE,legacy.axes=TRUE)
plot.roc(roc_list[[3]],add=TRUE,col=mycolor[3],plot=TRUE,legacy.axes=TRUE)
plot.roc(roc_list[[4]],add=TRUE,col=mycolor[4],plot=TRUE,legacy.axes=TRUE)
plot.roc(roc_list[[5]],add=TRUE,col=mycolor[5],plot=TRUE,legacy.axes=TRUE)
plot.roc(roc_list[[6]],add=TRUE,col=mycolor[6],plot=TRUE,legacy.axes=TRUE)
plot.roc(roc_list[[7]],add=TRUE,col=mycolor[7],plot=TRUE,legacy.axes=TRUE)
plot.roc(roc_list[[8]],add=TRUE,col=mycolor[8],plot=TRUE,legacy.axes=TRUE)

value_condition<-c("FLIPI1","FLIPI2","PRIMA-PI","LDH+B2mg","Adjust1","Adjust2","Adjust3","Adjust4")
# value_condition<-c("FLIPI1","FLIPI2","PRIMA-PI","Adjust1","Adjust2","Adjust3","Adjust4")
legend("bottomright", legend=value_condition, col=mycolor, lwd=2,cex = 0.8)
# legend("bottomright", legend=c("AUC of FLIPI1:0.5205","AUC of Breast: 0.621","AUC of Ovary_Fallopian_Tube: 0.61",
#                                "AUC of Skin: 0.307","AUC of CNS_Brain: 0.588", "AUC of Liver: 0.504", 
#                                "AUC of Uterus: 0.486","AUC of Lung: 0.912","AUC of Bladder_Urinary_Tract: 0.408",
#                                "AUC of ALL_tissue: 0.927"
#                                ), col=mycolor, lwd=2,cex = 0.8)

# title(main = "ROC")
dev.off()

library(modEvA)
#--------
pdf("./figure/08_aupr_FLIPI1_c.pdf")
aupr=AUC(obs=test$new_pod_total,pred=test$FLIPI1_c,curve = "PR", simplif=TRUE, main = "PR curve")
dev.off()
#-----
pdf("./figure/08_aupr_FLIPI2_c.pdf")
aupr=AUC(obs=test$new_pod_total,pred=test$FLIPI2_c,curve = "PR", simplif=TRUE, main = "PR curve")
dev.off()
#--------
pdf("./figure/08_aupr_primapi.pdf")
aupr=AUC(obs=test$new_pod_total,pred=test$primapi_re_c,curve = "PR", simplif=TRUE, main = "PR curve")
dev.off()
#---------

#---------
pdf("./figure/08_aupr_adjust1.pdf")
aupr=AUC(obs=adjust1$new_pod_total,pred=adjust1$Predict_Pod_total_f,curve = "PR", simplif=TRUE, main = "PR curve")
dev.off()
#---------
pdf("./figure/08_aupr_adjust2.pdf")
aupr=AUC(obs=adjust2$new_pod_total,pred=adjust2$Predict_Pod_total_f,curve = "PR", simplif=TRUE, main = "PR curve")
dev.off()
#---------
pdf("./figure/08_aupr_adjust3.pdf")
aupr=AUC(obs=adjust3$new_pod_total,pred=adjust3$Predict_Pod_total_f,curve = "PR", simplif=TRUE, main = "PR curve")
dev.off()
#-----------
pdf("./figure/08_aupr_adjust4.pdf")
aupr=AUC(obs=adjust4$new_pod_total,pred=adjust4$Predict_Pod_total_f,curve = "PR", simplif=TRUE, main = "PR curve")
dev.off()
#---------
pdf("./figure/08_aupr_b2mg_ldh.pdf")
aupr=AUC(obs=test$new_pod_total,pred=test$b2mg_ldh_c,curve = "PR", simplif=TRUE, main = "PR curve")
dev.off()

library(ggplot2)
library(plotROC)
p <- ggplot(test, aes(d = new_pod_total, m = FLIPI2_count)) + geom_roc()
p
dev.off()


# pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
# plot(pr)

# i=1
# for(tissue in value_condition){
#   file_name <-paste(tissue,"predict.txt",sep="_")
#   org<-read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
#   roc_list[[i]]<-roc(org$true_value1,org$predict_value1,direction=  "<")
#   print(i)
#   i=i+1
#   print(tissue)
# }