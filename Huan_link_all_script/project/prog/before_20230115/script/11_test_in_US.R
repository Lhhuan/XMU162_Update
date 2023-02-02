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
dat1 <- read.csv("../data/4C_Follicular_lymphoma_3A_79_cases_KHY_04-22-22.csv",na.strings = "")
dat<- dat1[-1,]

colnames(dat)<- c("ID","age_raw","Ki.67","stage","Bsym","LN_num","LN_long","BM","spleen","extend_num","SUVmax","SPD","ECOG","B2mg","B2M_upper","LDH_re0","LDH_upper","HGB","Lym",
"Mono","diagnosis","pro_time","followup","dead","OS_month","dead_time","OS_month1","TTP","pfs_month","pro_status","OS_year","TTP1","PFS_year","FLIPI_BINARY","PRIM_TX")

names <-data.frame(original_name= colnames(dat1),refine_name= colnames(dat))


count_pod <- function(x){
  if(is.na(x[1,22])){
    a = NA}else{
      a = difftime(x[1,22],x[1,21],units = "days")}
  return(a)}

re <- lapply(1:nrow(dat), function(i){
  print(i)
  x = dat[i,]
  a <- round(count_pod(x)/30)
  return(a)
}) 
dat$pod <- unlist(re)
dat$new_pod_total <-NA

dat$new_pod_total[is.na(dat$pod)]=0
dat$new_pod_total[dat$pod >24]=0
dat$new_pod_total[dat$pod <=24]=1

dat$age_s=NA
dat$age_s[dat$age_raw<=60]=0
dat$age_s[dat$age_raw>60]=1

dat$stage_s=NA
dat$stage_s[dat$stage<3]=0
dat$stage_s[dat$stage>=3]=1

dat$HGB_s =NA
dat$HGB_s[grep("≤120",dat$HGB)] <- 1
dat$HGB_s[grep(">120",dat$HGB)] <- 0

dat$LDH_s =dat$LDH_re0

dat$LN_num_s =NA
dat$LN_num_s[grep("≥4",dat$LN_num)] <- 1
dat$LN_num_s[grep("<4|3|1",dat$LN_num)] <- 0

dat$FLIPI1_count <- base::rowSums(dat[,c("age_s","stage_s","HGB_s","LDH_s","LN_num_s")], na.rm = TRUE)

dat$B2MG_re0 <-NA
dat$B2MG_re0[dat$B2mg >dat$B2M_upper]=1
dat$B2MG_re0[dat$B2mg <=dat$B2M_upper]=0
dat$B2mg_s <- dat$B2MG_re0

dat$LN6 <- NA
dat$LN6 [dat$LN_long >=6]=1
dat$LN6 [dat$LN_long <6]=0
dat$LN6_s <- dat$LN6

dat$BM_s <- as.numeric(dat$BM)
dat$FLIPI2_count <- base::rowSums(dat[,c("B2mg_s","LN6_s","BM_s","HGB_s","age_s")], na.rm = TRUE)

dat$primapi_re  <- NA
dat$primapi_re[dat$B2mg>3]= "High"
dat$primapi_re[dat$B2mg<=3 &dat$BM_s==1 ]="Intermediate"
dat$primapi_re[dat$B2mg<=3 &dat$BM_s==0 ]="Low"
dat$primapi_re[is.na(dat$primapi_re) & dat$pod_total==0 ]="Low"
dat$primapi_re[is.na(dat$primapi_re) & dat$pod_total>=3 ]="High"
dat$primapi_re[is.na(dat$primapi_re)]="Intermediate"

dat$primapi_re_n <- NA
dat$primapi_re_n[dat$primapi_re =="Low"]=0
dat$primapi_re_n[dat$primapi_re =="Intermediate"]=1
dat$primapi_re_n[dat$primapi_re =="High"]=2

dat$dead[grep("lived",dat$dead)] <- 0
dat$dead[grep("dead",dat$dead)] <- 1
dat$dead <-as.numeric(dat$dead)
dat$Lym_Mono <- dat$Lym/dat$Mono

dat$os_month_new <- dat$OS_month
dat$pfs_month_new <- dat$pfs_month

save(dat,file="11_US.Rdata")
#------------------------------------------------------------------new_cutoff
dat$LDH_c <-dat$LDH_re0
dat$B2M_c <-dat$B2MG_re0
dat$Lym_Mono_c <-NA
dat$Lym_Mono_c[dat$Lym_Mono <=10]=0
dat$Lym_Mono_c[dat$Lym_Mono >10 ]=1

dat$HGB_c =dat$HGB_s
dat$LN_num_c <-dat$LN_num_s
dat$SUVmax_c <-NA
dat$SUVmax_c[dat$SUVmax <=10]=0
dat$SUVmax_c[dat$SUVmax >10 ]=1


dat$BM_c <-dat$BM_s

dat$extend_num_c <-NA
dat$extend_num_c[dat$extend_num>0]=1
dat$extend_num_c[dat$extend_num==0]=0
# #-------------------------------------------------------------------------------------------------------------
dat$B2M_c_f <-dat$B2M_c*2
dat$Lym_Mono_c_f <-dat$Lym_Mono_c*2
dat$LDH_c_f <- dat$LDH_c *2
dat$HGB_c_f<- dat$HGB_c *1
dat$LN_num_c_f <-dat$LN_num_c *1
dat$SUVmax_c_f<- dat$SUVmax_c *1
dat$BM_c_f<- dat$BM_c *1
dat$sum_score<-base::rowSums(dat[,c("B2M_c_f","Lym_Mono_c_f","LDH_c_f","HGB_c_f","LN_num_c_f","SUVmax_c_f","BM_c_f")])

dat$Ours <- dat$sum_score

US <-dat

save(US,file="11_us_for_vaild.Rdata")

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
            
ggsave("./figure/11_pod24_test_in_US.png",p1,dpi=300,width=7,height=5.8)



# TP <-length(test$new_pod_total[test$adjust5>=5&test$new_pod_total==1])
# FP <-length(test$new_pod_total[test$adjust5>=5&test$new_pod_total==0])
# FN <-length(test$new_pod_total[test$adjust5<5&test$new_pod_total==1])
# TN <-length(test$new_pod_total[test$adjust5<5&test$new_pod_total==0])
# TPR = TP / (TP+FN)
# FPR = FP / (FP + TN)
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
# p1 <-p1 +annotate("text",x=0.35,y=TPR,label= "5",size=5,color="#E53935" )+geom_point(x=FPR,y=TPR ,size=4,color="#E53935" )

  p2 <-p1 +scale_x_continuous(limits= c(0, 1), breaks= seq(0,1,0.2))+ #expand = c(0, 0)
  scale_y_continuous(limits= c(0, 1), breaks= seq(0,1,0.2)) #+theme(panel.background = element_rect(color="black",size=1.2))
ggsave("./figure/11_pod24_test_in_US_NEW.png",p2,dpi=300,width=7,height=5.8)
pdf("aaa.pdf")
print(p2)
dev.off()


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
ggsave("./figure/11_OS_test_in_US_NEW.png",p2,dpi=300,width=7,height=5.8)
# pdf("aaa.pdf")
# print(p2)
# dev.off()

