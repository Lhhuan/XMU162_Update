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

# # dat$BM_c <-dat$BM_s

# # dat$extend_num_c <-NA
# # dat$extend_num_c[dat$extend_num>0]=1
# # dat$extend_num_c[dat$extend_num==0]=0
# # dat$extend_num_c[is.na(dat$extend_num_c)&dat$new_pod_total==1]=1
# # dat$extend_num_c[is.na(dat$extend_num_c)&dat$new_pod_total==0]=0


# dat$age_c =dat$age_s
# #------------------------------------------------------------------------------------------------------------
dat$Lym_Mono_c_f <-dat$Lym_Mono_c*2
dat$LDH_c_f <- dat$LDH_c *1
dat$HGB_c_f<- dat$HGB_c *1
dat$B2M_c_f <-dat$B2M_c*1
dat$SUVmax_c_f <-dat$SUVmax_c*1
dat$LN_num_c_f <-dat$LN_num_c *1

dat$sum_score<-base::rowSums(dat[,c("Lym_Mono_c_f","LDH_c_f","HGB_c_f","B2M_c_f","SUVmax_c_f","LN_num_c_f")])

dat$Ours <- dat$sum_score
save(dat,file="10_3_ALL_data_OS_valid.Rdata")

dat1 <-dat
load("/home/huanhuan/project/prog/data/sex.Rdata")
# dat1 <-left_join(dat1,sex,by="No")
dat1$class <- NA
dat1[test_set_number,"class"] ="test"
dat1[train_set_number,"class"] ="train"
save(dat1,file="final_3a_0_data20230116.Rdata")
write.table(dat1,"final_3a_data20230116.txt",row.names = F, col.names = T,quote =F,sep="\t")
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.background = element_rect(color="black",size=1.2),
                                                  axis.title.y = element_text(size = 15),
                                                  axis.title.x = element_text(size = 15),
                                                  axis.text.y = element_text(size = 12,colour = "black"),
                                                  axis.text.x = element_text(size = 12,colour = "black"),
                                                  plot.title = element_text(hjust = 0.5),legend.title=element_blank(),legend.text=element_text(size=9))

mycolor <-c("#827717","#1B5E20","#006064","#E53935","#01579B")
#--------------------------------------------------------------------------8
library(plotROC)
longtest <-melt_roc(dat,"dead",c("FLIPI1_count","FLIPI2_count","primapi_re_n","b2mg_ldh1.5s","Ours"))
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
              names(auc)[5]," AUC: ",round(auc[5],3),"\n"),
              size=3.5)
            
ggsave("./figure/10_3_OS_dead_roc5_whole_data.png",p1,dpi=300,width=7,height=5.8)


mycolor <-c("#827717","#1B5E20","#006064","#E53935","#01579B")
# test$Ours <-test$adjust5

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
# p1 <-p1 +annotate("text",x=0.35,y=TPR,label= "5",size=5,color="#E53935" )+geom_point(x=FPR,y=TPR ,size=4,color="#E53935" )

  p2 <-p1 +scale_x_continuous(limits= c(0, 1), breaks= seq(0,1,0.2))+ #expand = c(0, 0)
  scale_y_continuous(limits= c(0, 1), breaks= seq(0,1,0.2)) #+theme(panel.background = element_rect(color="black",size=1.2))
ggsave("./figure/10_3_OS_dead_roc5_whole_data_new.png",p2,dpi=300,width=7,height=5.8)
pdf("bbb.pdf")
print(p2)
dev.off()

