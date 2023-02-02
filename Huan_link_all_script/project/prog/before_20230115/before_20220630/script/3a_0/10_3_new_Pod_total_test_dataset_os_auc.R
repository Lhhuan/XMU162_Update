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
dat$B2M_c <-dat$B2MG_re0

dat$ECOG_c =NA
dat$ECOG_c[dat$ECOG <=1]=0
dat$ECOG_c[dat$ECOG >1]=1
dat$ECOG_c[is.na(dat$ECOG_c)&dat$new_pod_total==1]=1
dat$ECOG_c[is.na(dat$ECOG_c)&dat$new_pod_total==0]=0


dat$Lym_Mono_c <-NA
dat$Lym_Mono_c[dat$Lym_Mono <=10]=0
dat$Lym_Mono_c[dat$Lym_Mono >10 ]=1
dat$Lym_Mono_c[is.na(dat$Lym_Mono_c)&dat$new_pod_total==1]=1
dat$Lym_Mono_c[is.na(dat$Lym_Mono_c)&dat$new_pod_total==0]=0
dat$BM_c <-dat$BM_s


dat$SUVmax_c <-NA
dat$SUVmax_c[dat$SUVmax <=10]=0
dat$SUVmax_c[dat$SUVmax >10 ]=1
dat$SUVmax_c[is.na(dat$SUVmax_c)&dat$new_pod_total==1]=1
dat$SUVmax_c[is.na(dat$SUVmax_c)&dat$new_pod_total==0]=0

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


#-----------------------------------------------adjust
dat$LDH_c_f <- dat$LDH_c*2
dat$LN_num_c_f <-dat$LN_num_c *2
dat$HGB_c_f<- dat$HGB_c*2
dat$B2M_c_f <-dat$B2M_c*1
dat$ECOG_c_f<- dat$ECOG_c *1
dat$Lym_Mono_c_f <-dat$Lym_Mono_c*1
dat$BM_c_f <-dat$BM_c*1
dat$SUVmax_c_f<- dat$SUVmax_c *1

dat$sum_score<-base::rowSums(dat[,c("LDH_c_f","LN_num_c_f","HGB_c_f","B2M_c_f","ECOG_c_f","Lym_Mono_c_f","BM_c_f","SUVmax_c_f")])
# dat[is.na($sum_score),]
dat$adjust4 <- dat$sum_score


p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.background = element_rect(color="black",size=1.2),
                                                  axis.title.y = element_text(size = 15),
                                                  axis.title.x = element_text(size = 15),
                                                  axis.text.y = element_text(size = 12,colour = "black"),
                                                  axis.text.x = element_text(size = 12,colour = "black"),
                                                  plot.title = element_text(hjust = 0.5),legend.title=element_blank(),legend.text=element_text(size=9))


#------------------------------------------------------
#--------------------------------------------------------------------------adjust
TP <-length(dat$new_pod_total[dat$adjust5>=4&dat$new_pod_total==1])
FP <-length(dat$new_pod_total[dat$adjust5>=4&dat$new_pod_total==0])
FN <-length(dat$new_pod_total[dat$adjust5<4&dat$new_pod_total==1])
TN <-length(dat$new_pod_total[dat$adjust5<4&dat$new_pod_total==0])
TPR = TP / (TP+FN)
FPR = FP / (FP + TN)
# mycolor <-c("#E53935","#827717","#1B5E20","#006064","#01579B")
mycolor <-c("#827717","#1B5E20","#006064","#E53935","#01579B")
dat$Ours <-dat$adjust4

library(plotROC)
longtest <-melt_roc(dat,"dead",c("FLIPI1_count","FLIPI2_count","primapi_re_n","b2mg_ldh1.5s","Ours"))
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
# p1 <-p1 +annotate("text",x=0.39,y=0.75,label= "4",size=5,color="#E53935" )+geom_point(x=FPR,y=TPR ,size=4,color="#E53935" )

  p2 <-p1 +scale_x_continuous(limits= c(0, 1), breaks= seq(0,1,0.2))+ #expand = c(0, 0)
  scale_y_continuous(limits= c(0, 1), breaks= seq(0,1,0.2)) #+theme(panel.background = element_rect(color="black",size=1.2))
ggsave("./figure/10_3_roc_Ours_ALL_dead.png",p2,dpi=300,width=7,height=5.8)
pdf("ALL_os.pdf")
print(p2)
dev.off()
#-----------
test <-dat[test_set_number,]


library(plotROC)
longtest <-melt_roc(test,"dead",c("FLIPI1_count","FLIPI2_count","primapi_re_n","b2mg_ldh1.5s","Ours"))
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
# p1 <-p1 +annotate("text",x=0.39,y=0.75,label= "4",size=5,color="#E53935" )+geom_point(x=FPR,y=TPR ,size=4,color="#E53935" )

  p2 <-p1 +scale_x_continuous(limits= c(0, 1), breaks= seq(0,1,0.2))+ #expand = c(0, 0)
  scale_y_continuous(limits= c(0, 1), breaks= seq(0,1,0.2)) #+theme(panel.background = element_rect(color="black",size=1.2))
ggsave("./figure/10_3_roc_adjust5_test_dead.png",p2,dpi=300,width=7,height=5.8)
pdf("test_os.pdf")
print(p2)
dev.off()


#------------------------------------------------------------------------FLIPI


custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}
#-----------------
#--------------------------------------------------------------------------------------------adjust5
# cutoff_y=4
# test$adjust5_01 <-NA
# test$adjust5_01[test$adjust5<cutoff_y ]=0
# test$adjust5_01[test$adjust5>=cutoff_y]=1
# library("survival")
# library("survminer")

# fit <- survfit(Surv(pfs_month_new, pro_status) ~ adjust5_01, data=test)
# # pdf("./figure/10_3_pod24_pfs_cutoff1_adjust5.pdf",height=5.2,width=5)
# p1 <- ggsurvplot(fit,
#                 palette=c("#21618C","#B71C1C"),
#                 legend.labs=c("Low risk","High risk"), #标签
#                 pval = TRUE,
#                 legend.title="POD24",
#                 title="PFS",
#                 xlab = " Time (Months)",
#                 xlim=c(0,110),
#                 ggtheme = custom_theme()
#                 )
# ggsave("./figure/10_3_pod24_pfs_cutoff1_adjust5_new.png",p1$plot,dpi=300,height=5.2,width=5)

# fit <- survfit(Surv(os_month_new, dead) ~ adjust5_01, data=test)
# # pdf("./figure/10_3_pod24_Os_cutoff1_adjust5.pdf",height=5.2,width=5)
# p1 <- ggsurvplot(fit,
#                 palette=c("#21618C","#B71C1C"),
#                 legend.labs=c("Low risk","High risk"), #标签
#                 pval = TRUE,
#                 legend.title="POD24",
#                 title="OS",
#                 xlab = " Time (Months)",
#                 xlim=c(0,110),
#                 ggtheme = custom_theme()
#                 )
# ggsave("./figure/10_3_pod24_Os_cutoff1_adjust5_new.png",p1$plot,dpi=300,height=5.2,width=5)

