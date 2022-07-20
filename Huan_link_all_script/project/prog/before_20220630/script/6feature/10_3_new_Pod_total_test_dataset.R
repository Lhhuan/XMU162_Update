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
dat$B2M_c <-dat$B2MG_re0
dat$Lym_Mono_c <-NA
dat$Lym_Mono_c[dat$Lym_Mono <=10]=0
dat$Lym_Mono_c[dat$Lym_Mono >10 ]=1
dat$Lym_Mono_c[is.na(dat$Lym_Mono_c)&dat$new_pod_total==1]=1
dat$Lym_Mono_c[is.na(dat$Lym_Mono_c)&dat$new_pod_total==0]=0
dat$LDH_c <-dat$LDH_re0
dat$HGB_c =dat$HGB_s
dat$LN_num_c <-dat$LN_num_s
dat$SUVmax_c <-NA
dat$SUVmax_c[dat$SUVmax <=10]=0
dat$SUVmax_c[dat$SUVmax >10 ]=1
dat$SUVmax_c[is.na(dat$SUVmax_c)&dat$new_pod_total==1]=1
dat$SUVmax_c[is.na(dat$SUVmax_c)&dat$new_pod_total==0]=0

dat$BM_c <-dat$BM_s

dat$extend_num_c <-NA
dat$extend_num_c[dat$extend_num>0]=1
dat$extend_num_c[dat$extend_num==0]=0
dat$extend_num_c[is.na(dat$extend_num_c)&dat$new_pod_total==1]=1
dat$extend_num_c[is.na(dat$extend_num_c)&dat$new_pod_total==0]=0

setwd("/home/huanhuan/project/prog/script/6feature/output/")
# dat$age_c =dat$age_s


# dat$BM_extend_c <-dat$BM_extend
# dat$BM_extend_c[is.na(dat$BM_extend_c)&dat$new_pod_total==1]=1
# dat$BM_extend_c[is.na(dat$BM_extend_c)&dat$new_pod_total==0]=0

# dat$ECOG_c =NA
# dat$ECOG_c[dat$ECOG <=1]=0
# dat$ECOG_c[dat$ECOG >1]=1
# dat$ECOG_c[is.na(dat$ECOG_c)&dat$new_pod_total==1]=1
# dat$ECOG_c[is.na(dat$ECOG_c)&dat$new_pod_total==0]=0

# dat$SPD_c <-NA
# dat$SPD_c[dat$SPD <=20]=0
# dat$SPD_c[dat$SPD >20]=1
# dat$SPD_c[is.na(dat$SPD_c)&dat$new_pod_total==1]=1
# dat$SPD_c[is.na(dat$SPD_c)&dat$new_pod_total==0]=0

# #-------------------------------------------------------------------------------------------------------------
test <-dat[test_set_number,]
# im <-read.table("/share/Projects/huanhuan/project/prog/10_2_feature_ratio.txt",header=T,sep="\t")

# test$sum <-base::rowSums(test[,c("Lym_Mono","B2mg","SPD","LDH","LN_num","SUVmax","Ki.67","HGB","BM","Bsym")])
# base::rowSums(dat_s, na.rm = TRUE)
#-----------------------------------------------adjust2
test$Lym_Mono_c_f <-test$Lym_Mono_c*9
test$LDH_c_f <- test$LDH_c *8
test$B2M_c_f <-test$B2M_c*8
test$LN_num_c_f <-test$LN_num_c *5
test$HGB_c_f<- test$HGB_c *3
test$BM_c_f<- test$BM_c *2
test$sum_score<-base::rowSums(test[,c("Lym_Mono_c_f","LDH_c_f","B2M_c_f","LN_num_c_f","HGB_c_f","BM_c_f")])

test$adjust2 <- test$sum_score

#-----------------------------------------------adjust3
test$Lym_Mono_c_f <-test$Lym_Mono_c*4
test$LDH_c_f <- test$LDH_c *4
test$B2M_c_f <-test$B2M_c*4
test$LN_num_c_f <-test$LN_num_c *3
test$HGB_c_f<- test$HGB_c *2
test$BM_c_f<- test$BM_c *1
test$sum_score<-base::rowSums(test[,c("Lym_Mono_c_f","LDH_c_f","B2M_c_f","LN_num_c_f","HGB_c_f","BM_c_f")])

test$adjust3 <- test$sum_score
#-----------------------------------------------adjust4
test$Lym_Mono_c_f <-test$Lym_Mono_c*2
test$LDH_c_f <- test$LDH_c *2
test$B2M_c_f <-test$B2M_c*2
test$LN_num_c_f <-test$LN_num_c *1
test$HGB_c_f<- test$HGB_c *1
test$BM_c_f<- test$BM_c *1
test$sum_score<-base::rowSums(test[,c("Lym_Mono_c_f","LDH_c_f","B2M_c_f","LN_num_c_f","HGB_c_f","BM_c_f")])

test$adjust4 <- test$sum_score

#----------------------------------------------
test$Ours <-test$adjust4
save(test,file="10_3_test_vaild.Rdata")
# setwd("/share/Projects/huanhuan/project/prog/")

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
longtest <-melt_roc(test,"new_pod_total",c("FLIPI1_count","FLIPI2_count","primapi_re_n","b2mg_ldh1.5s","adjust3","adjust4","adjust2"))
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
#--------------------------------------------------------------------------Ours
test$Ours <-test$adjust4
TP <-length(test$new_pod_total[test$Ours>=3&test$new_pod_total==1])
FP <-length(test$new_pod_total[test$Ours>=3&test$new_pod_total==0])
FN <-length(test$new_pod_total[test$Ours<3&test$new_pod_total==1])
TN <-length(test$new_pod_total[test$Ours<3&test$new_pod_total==0])
TPR = TP / (TP+FN)
FPR = FP / (FP + TN)
mycolor <-c("#827717","#1B5E20","#006064","#E53935","#01579B")


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
# p1 <-p1 +annotate("text",x=0.3,y=TPR,label= "3",size=5,color="#E53935" )+geom_point(x=FPR,y=TPR ,size=4,color="#E53935" )

  p2 <-p1 +scale_x_continuous(limits= c(0, 1), breaks= seq(0,1,0.2))+ #expand = c(0, 0)
  scale_y_continuous(limits= c(0, 1), breaks= seq(0,1,0.2)) #+theme(panel.background = element_rect(color="black",size=1.2))
ggsave("./figure/10_3_roc_adjust5_NEW.png",p2,dpi=300,width=7,height=5.8)
pdf("aaa.pdf")
print(p2)
dev.off()


custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}
#--------------------------------------------------------------------------------------------ours
cutoff_y=3
test$Ours_01 <-NA
test$Ours_01[test$Ours<cutoff_y ]=0
test$Ours_01[test$Ours>=cutoff_y]=1
library("survival")
library("survminer")

fit <- survfit(Surv(pfs_month_new, pro_status) ~ Ours_01, data=test)
# pdf("./figure/10_3_pod24_pfs_cutoff1_adjust3.pdf",height=5.2,width=5)
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
ggsave("./figure/10_3_pod24_pfs_cutoff1_Ours.png",p1$plot,dpi=300,height=5.2,width=5)
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ Ours_01, data=test)
# pdf("./figure/10_3_pod24_Os_cutoff1_adjust3.pdf",height=5.2,width=5)
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
ggsave("./figure/10_3_pod24_Os_cutoff1_Ours.png",p1$plot,dpi=300,height=5.2,width=5)
#-----------------

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
