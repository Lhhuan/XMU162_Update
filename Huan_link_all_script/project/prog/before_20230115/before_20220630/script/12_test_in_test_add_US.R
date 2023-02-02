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
load("10_3_ALL_data_OS_valid.Rdata")
all <-dat

load("10_3_test_vaild.Rdata")
load("11_us_for_vaild.Rdata")

dat <- rbind(test[,c("FLIPI1_count","FLIPI2_count","primapi_re_n","Ours","new_pod_total","pro_status","pfs_month_new","dead","os_month_new")],
  US[,c("FLIPI1_count","FLIPI2_count","primapi_re_n","Ours","new_pod_total","pro_status","pfs_month_new","dead","os_month_new")])

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.background = element_rect(color="black",size=1.2),
                                                  axis.title.y = element_text(size = 15),
                                                  axis.title.x = element_text(size = 15),
                                                  axis.text.y = element_text(size = 12,colour = "black"),
                                                  axis.text.x = element_text(size = 12,colour = "black"),
                                                  plot.title = element_text(hjust = 0.5),legend.title=element_blank(),legend.text=element_text(size=9))

mycolor <-c("#827717","#1B5E20","#E53935","#006064","#01579B")
#--------------------------------------------------------------------------8
library(plotROC)
longtest <-melt_roc(dat,"new_pod_total",c("FLIPI1_count","FLIPI2_count","primapi_re_n","Ours"))
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
              names(auc)[4]," AUC: ",round(auc[4],3),"\n"),
              size=3.5)
            
ggsave("./figure/12_pod24_test_in_merge.png",p1,dpi=300,width=7,height=5.8)



TP <-length(dat$new_pod_total[dat$Ours>=3&dat$new_pod_total==1])
FP <-length(dat$new_pod_total[dat$Ours>=3&dat$new_pod_total==0])
FN <-length(dat$new_pod_total[dat$Ours<3&dat$new_pod_total==1])
TN <-length(dat$new_pod_total[dat$Ours<3&dat$new_pod_total==0])
TPR = TP / (TP+FN)
FPR = FP / (FP + TN)
mycolor <-c("#827717","#1B5E20","#E53935","#006064","#01579B")
# test$Ours <-test$adjust5

library(plotROC)
longtest <-melt_roc(dat,"new_pod_total",c("FLIPI1_count","FLIPI2_count","primapi_re_n","Ours"))
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
              names(auc)[4]," AUC: ",round(auc[4],3),"\n"),
              size=3.5)
p1 <-p1 +annotate("text",x=0.38,y=TPR,label= "3",size=5,color="#E53935" )+geom_point(x=FPR,y=TPR ,size=4,color="#E53935" )

  p2 <-p1 +scale_x_continuous(limits= c(0, 1), breaks= seq(0,1,0.2))+ #expand = c(0, 0)
  scale_y_continuous(limits= c(0, 1), breaks= seq(0,1,0.2)) #+theme(panel.background = element_rect(color="black",size=1.2))
ggsave("./figure/12_pod24_test_in_merge_new.png",p2,dpi=300,width=7,height=5.8)
pdf("aaa.pdf")
print(p2)
dev.off()

#------------------------------
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}

cutoff_y=3
dat$Ours_01 <-NA
dat$Ours_01[dat$Ours<cutoff_y ]=0
dat$Ours_01[dat$Ours>=cutoff_y]=1
library("survival")
library("survminer")

fit <- survfit(Surv(pfs_month_new, pro_status) ~ Ours_01, data=dat)
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
ggsave("./figure/12_merge_vaild_pod24_pfs_Ours.png",p1$plot,dpi=300,height=5.2,width=5)
#------------------------------
fit <- survfit(Surv(os_month_new, dead) ~ Ours_01, data=dat)
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
ggsave("./figure/12_merge_vaild_pod24_Os_Ours.png",p1$plot,dpi=300,height=5.2,width=5)






dat <- rbind(all[,c("FLIPI1_count","FLIPI2_count","primapi_re_n","Ours","new_pod_total","pro_status","pfs_month_new","dead","os_month_new")],
  US[,c("FLIPI1_count","FLIPI2_count","primapi_re_n","Ours","new_pod_total","pro_status","pfs_month_new","dead","os_month_new")])


library(plotROC)
longtest <-melt_roc(dat,"dead",c("FLIPI1_count","FLIPI2_count","primapi_re_n","Ours"))
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
              names(auc)[4]," AUC: ",round(auc[4],3),"\n"),
              size=3.5)
# p1 <-p1 +annotate("text",x=0.35,y=TPR,label= "5",size=5,color="#E53935" )+geom_point(x=FPR,y=TPR ,size=4,color="#E53935" )

  p2 <-p1 +scale_x_continuous(limits= c(0, 1), breaks= seq(0,1,0.2))+ #expand = c(0, 0)
  scale_y_continuous(limits= c(0, 1), breaks= seq(0,1,0.2)) #+theme(panel.background = element_rect(color="black",size=1.2))
ggsave("./figure/12_merge_all_us_vaild_OS_NEW.png",p2,dpi=300,width=7,height=5.8)
# pdf("aaa.pdf")
# print(p2)
# dev.off()

