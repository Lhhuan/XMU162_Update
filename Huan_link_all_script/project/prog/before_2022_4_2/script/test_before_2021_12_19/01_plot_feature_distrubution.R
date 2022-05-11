
setwd("D:\\Huan_R\\prog\\data")
library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)
library(cowplot)
library(ggpubr)
library(Hmisc)
data <- read.csv("data.csv",na.strings = "")
colnames(data)[78:79] <- c("pro_time","pro_status")
# data$followup_status[grep("^£¿|",data$followup_status)] <- NA
# data[,c("diagnosis","dead_time","pro_time")]<-apply(data[,c("diagnosis","dead_time","pro_time")],2,as.Date)
data$diagnosis<-as.Date(data$diagnosis)
data$dead_time<-as.Date(data$dead_time)
data$pro_time<-as.Date(data$pro_time)
data$followup<-as.Date(data$followup)
data$pfs_month_new <-difftime(data$pro_time,data$diagnosis,units = c("days"))%>%as.numeric()
data$os_month_new <-difftime(data$dead_time,data$diagnosis,units = c("days"))%>%as.numeric()
data$os_month_new[grep("0",data$dead)] <-difftime(data$followup[grep("0",data$dead)],data$diagnosis[grep("0",data$dead)],units = c("days"))%>%as.numeric()


data$pfs_month_new <-data$pfs_month_new/30
data$os_month_new <-data$os_month_new/30
org_sur <-data[,c("No","dead_time","pfs_month_new","os_month_new","pro_status","pro_time","followup")]
#-----------
org_age<-data%>%select("No","birth","diagnosis")
org_age$diagnosis[grep(" ",org_age$diagnosis)] <-NA
org_age$birth <-as.Date(org_age$birth)
org_age$diagnosis <-as.Date(org_age$diagnosis)
org_age$age <-difftime(org_age$diagnosis,org_age$birth,units = c("days"))
org_age$age <-as.numeric(org_age$age)
org_age$age <-org_age$age/365
org_age <-org_age[,-c(2,3)]

load("step2_filter_sample_and_features.Rdata")
org <-dat_filter
org<-org[,-c(2:5)]
org <-left_join(org,org_age,by="No")
org$grade[grep("3$",org$grade)] <-"3a"
save(org,file = "01_refine_grade_add_age.Rdata")
org1<-merge(org,org_sur,by="No")
org_q <-left_join(org,org_sur,by="No")
save(org,file = "01_refine_grade_add_man_age_pfs_os.Rdata")
org_grade <-org[!is.na(org$grade),]

not_used_var <-c("No","FLIPI1","FLIPI2","primapi",
       "B2mg0","LDH0","HGB0","trans_time","diagnosis","trans_time","relapse_time","followup","dead")

potential_feature <-c("X150b2mg_ldh","b2mg_LDH","CR","Rmaintain","progress","trans","relapse_res")
org_plot_feature <-org_grade%>%select(-c("No","FLIPI1","FLIPI2","primapi","B2mg0","LDH0","HGB0",
                                         "trans_time","diagnosis","trans_time","relapse_time","followup","dead",
                                         "pod_total"))
org_plot_feature <-org_plot_feature%>%select(-c("X150b2mg_ldh","b2mg_LDH","CR","Rmaintain","progress",
                                                "trans","relapse_res","site0","BM","spleen","BM_extend",
                                                "LN6","relapes","gender","Bsym"))

plot_feature_range <-c(1,4:23)
continuous_values<-c(1,23,3,5,7:9,11:17,21,22)

class_values<-setdiff(plot_feature_range,continuous_values)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.background = element_blank(), 
                                                  # axis.title.y = element_text(size = 8),
                                                  # axis.title.x = element_text(size = 10),
                                                  axis.line = element_line(colour = "black"),legend.position ="none",plot.title = element_text(hjust = 0.5))
org_plot_feature1 <-org_plot_feature
org_plot_feature1[,c(3:ncol(org_plot_feature1))]<-apply(org_plot_feature1[,c(3:ncol(org_plot_feature1))],2,as.numeric)
my_comparisons <- list(c("0","3a"),c("0","3b"),c("3a","3b"))

# for(i in plot_feature_range){
plot_f <-function(i){
  colnames(org_plot_feature1)[i]<-"aaaa"
  title_name=colnames(org_plot_feature)[i]
  title_name <-str_replace(title_name,"_"," ")
  title_name <-capitalize(title_name)
  p <-ggplot(org_plot_feature1,aes(x=grade,y=aaaa))+geom_boxplot(aes(fill=grade),width=0.3,outlier.colour = NA)+
    # p <-ggplot(org_plot_feature1,aes(x=grade,y=aaaa))+geom_boxplot(aes(fill=grade),width=0.3)+
    theme_bw()+p_theme+ylab("")+xlab("Grade")+ggtitle(title_name) 
  p1 <- p + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif", method.args = list(alternative = "two.sided"))
  return(p1)
}

plist <-lapply(continuous_values,plot_f)
pdf("./figure/01_continuous_variable1_boxplot.pdf",width = 10,height = 11)
CombinePlots(plist,ncol=3,nrow=3)
# plot_grid(plist=plist,nrow=3,ncol=3,byrow=TRUE) 
dev.off()
pdf("./figure/01_continuous_variable2_boxplot.pdf",width = 10,height = 11)
CombinePlots(plist[9:length(plist)],ncol=3,nrow=3)
# plot_grid(plist=plist,nrow=3,ncol=3,byrow=TRUE) 
dev.off()
plist <-lapply(class_values,plot_f)
pdf("./figure/classified_variable1.pdf",width = 10,height = 11)
CombinePlots(plist,ncol=3,nrow=3)
dev.off()
#--------------------
plist <-lapply(continuous_values,plot_f)
pdf("./figure/continuous_variable1_outlier.pdf",width = 10,height = 7)
CombinePlots(plist,ncol=3,nrow=2)
# plot_grid(plist=plist,nrow=3,ncol=3,byrow=TRUE) 
dev.off()
pdf("./figure/continuous_variable2_outlier.pdf",width = 10,height = 7)
CombinePlots(plist[6:length(plist)],ncol=3,nrow=2)
# plot_grid(plist=plist,nrow=3,ncol=3,byrow=TRUE) 
dev.off()

plist <-lapply(class_values,plot_f)
pdf("./figure/classified_variable1_outlier.pdf",width = 10,height = 11)
CombinePlots(plist,ncol=3,nrow=3)
dev.off()
#-----------------------------------------------violin 
plot_f_v <-function(i){
  colnames(org_plot_feature1)[i]<-"aaaa"
  title_name=colnames(org_plot_feature)[i]
  title_name <-str_replace(title_name,"_"," ")
  title_name <-capitalize(title_name)
  # p <-ggplot(org_plot_feature1,aes(x=grade,y=aaaa))+geom_boxplot(aes(fill=grade),width=0.3,outlier.colour = NA)+
  p <-ggplot(org_plot_feature1,aes(x=grade,y=aaaa,fill=grade))+geom_violin(width=1) +geom_boxplot(width=0.05)+
    theme_bw()+p_theme+ylab("")+xlab("Grade")+ggtitle(title_name) 
  p1 <- p + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif", method.args = list(alternative = "two.sided"))
  return(p1)
}
plist <-lapply(continuous_values,plot_f_v)
pdf("./figure/01_continuous_variable1_outlier_violin.pdf",width = 10,height = 11)
CombinePlots(plist,ncol=3,nrow=3)
# plot_grid(plist=plist,nrow=3,ncol=3,byrow=TRUE) 
dev.off()
pdf("./figure/01_continuous_variable2_outlier_violin.pdf",width = 10,height = 11)
CombinePlots(plist[9:length(plist)],ncol=3,nrow=3)
# plot_grid(plist=plist,nrow=3,ncol=3,byrow=TRUE) 
dev.off()


plist <-lapply(class_values,plot_f_v)
pdf("./figure/01_classified_variable1_outlier_violin.pdf",width = 10,height = 7)
CombinePlots(plist,ncol=3,nrow=2)
dev.off()
#--------------------------------------------------------------------violin_without_outlier
plot_f_v_o <-function(i){
  colnames(org_plot_feature1)[i]<-"aaaa"
  title_name=colnames(org_plot_feature)[i]
  title_name <-str_replace(title_name,"_"," ")
  title_name <-capitalize(title_name)
  # p <-ggplot(org_plot_feature1,aes(x=grade,y=aaaa))+geom_boxplot(aes(fill=grade),width=0.3,outlier.colour = NA)+
  p <-ggplot(org_plot_feature1,aes(x=grade,y=aaaa,fill=grade))+geom_violin(width=1) +geom_boxplot(width=0.05,outlier.colour = NA)+
    theme_bw()+p_theme+ylab("")+xlab("Grade")+ggtitle(title_name) 
  p1 <- p + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif", method.args = list(alternative = "two.sided"))
  return(p1)
}
#-----------------------

plist <-lapply(continuous_values,plot_f_v_o)
pdf("./figure/01_continuous_variable1_violin.pdf",width = 10,height = 11)
CombinePlots(plist,ncol=3,nrow=3)
# plot_grid(plist=plist,nrow=3,ncol=3,byrow=TRUE) 
dev.off()
pdf("./figure/01_continuous_variable2_violin.pdf",width = 10,height = 11)
CombinePlots(plist[9:length(plist)],ncol=3,nrow=3)
# plot_grid(plist=plist,nrow=3,ncol=3,byrow=TRUE) 
dev.off()


plist <-lapply(class_values,plot_f_v)
pdf("./figure/01_classified_variable1_violin.pdf",width = 10,height = 7)
CombinePlots(plist,ncol=3,nrow=2)
dev.off()
#------------------------------------------------------------pie

org_plot_feature1<-org_plot_feature[,class_values]
mycolor<-c( "#70A1D7","#C86B85","#FFD2A5", "#C79ECF", "#C1C0B9", "#A1DE93","#A2D2FF")
for(j in c(1,4:6)){
  org_plot_feature1$grade <-org_plot_feature$grade
  cc<-org_plot_feature1[!is.na(org_plot_feature1[j]),]
  colnames(cc)[j]="aaa" #medi
  aa <-group_by(cc,grade,aaa)%>%summarise(count=n())%>%as.data.frame()
  bb <-group_by(cc,grade)%>%summarise(all_count=n())%>%as.data.frame()
  final_u_p <-left_join(aa,bb,by="grade")
  final_u_p$ratio <-final_u_p$count/final_u_p$all_count*100
  final_u_p$ratio <-round(final_u_p$ratio,2)
  final_u_p$ratio <-paste0(final_u_p$ratio,"%")
  # final_u_p$ratio <-paste(final_u_p$stage,final_u_p$ratio,sep=": ")
  #---------------
  colnames(cc)[j]=colnames(org_plot_feature1)[j]
  fea=colnames(cc)[j]
  fea <-capitalize(fea)
  fea1 <-str_replace(fea,"_"," ")
  p_stage <-function(i){
    final_u_p1 <-filter(final_u_p,grade==i)
    m=paste0(i,": ",fea1)
    f_name <-paste0("./figure/01_classified_var_pie/",fea,"_",i,".pdf")
    pdf(f_name,width=7,height = 7)
    pie(final_u_p1$count,final_u_p1$ratio,main=m,col = mycolor[1:length(final_u_p1$count)])
    legend("topright",c("1","2","3","4"),fill = mycolor[1:length(final_u_p1$count)])
    dev.off()
  }
  lapply(c(0,"3a","3b"),p_stage)
}
#------------------------------------extend
j=2
  org_plot_feature1$grade <-org_plot_feature$grade
  cc<-org_plot_feature1[!is.na(org_plot_feature1[j]),]
  colnames(cc)[j]="aaa" #medi
  aa <-group_by(cc,grade,aaa)%>%summarise(count=n())%>%as.data.frame()
  bb <-group_by(cc,grade)%>%summarise(all_count=n())%>%as.data.frame()
  final_u_p <-left_join(aa,bb,by="grade")
  final_u_p$ratio <-final_u_p$count/final_u_p$all_count*100
  final_u_p$ratio <-round(final_u_p$ratio,2)
  final_u_p$ratio <-paste0(final_u_p$ratio,"%")
  # final_u_p$ratio <-paste(final_u_p$stage,final_u_p$ratio,sep=": ")
  #---------------
  colnames(cc)[j]=colnames(org_plot_feature1)[j]
  fea=colnames(cc)[j]
  fea <-capitalize(fea)
  mycolor<-c( "#70A1D7","#C86B85","#FFD2A5", "#C79ECF", "#C1C0B9", "#A1DE93")
  p_stage <-function(i){
    final_u_p1 <-filter(final_u_p,grade==i)
    m=paste0(i,": ",fea)
    f_name <-paste0("./figure/01_classified_var_pie/",fea,"_",i,".pdf")
    pdf(f_name,width=7,height = 7)
    le = c("0","1","2","3","4","5")
    pie(final_u_p1$count,final_u_p1$ratio,main=m,col = mycolor[1:length(final_u_p1$count)])
    legend("topright",le,fill = mycolor[1:length(le)])
    dev.off()
  }
  lapply(c(0,"3a"),p_stage)
  
i="3b"

final_u_p1 <-filter(final_u_p,grade==i)
m=paste0(i,": ",fea)
f_name <-paste0("./figure/01_classified_var_pie/",fea,"_",i,".pdf")
pdf(f_name,width=7,height = 7)
le = c("0","1","2","3","4","5")
pie(final_u_p1$count,final_u_p1$ratio,main=m,col = mycolor[c(1:4,6)])
# my_c2 = c("#70A1D7","#C86B85","#FFD2A5","#C1C0B9")
# pie(final_u_p1$count,final_u_p1$ratio,main=m,col =my_c2 )
legend("topright",le,fill = mycolor[1:length(le)])
dev.off()

#--------------------------
j=3
org_plot_feature1$grade <-org_plot_feature$grade
cc<-org_plot_feature1[!is.na(org_plot_feature1[j]),]
colnames(cc)[j]="aaa" #medi
aa <-group_by(cc,grade,aaa)%>%summarise(count=n())%>%as.data.frame()
bb <-group_by(cc,grade)%>%summarise(all_count=n())%>%as.data.frame()
final_u_p <-left_join(aa,bb,by="grade")
final_u_p$ratio <-final_u_p$count/final_u_p$all_count*100
final_u_p$ratio <-round(final_u_p$ratio,2)
final_u_p$ratio <-paste0(final_u_p$ratio,"%")
# final_u_p$ratio <-paste(final_u_p$stage,final_u_p$ratio,sep=": ")
#---------------
colnames(cc)[j]=colnames(org_plot_feature1)[j]
fea=colnames(cc)[j]
fea <-capitalize(fea)
mycolor<-c( "#70A1D7","#C86B85","#FFD2A5", "#C79ECF", "#C1C0B9", "#A1DE93")
p_stage <-function(i){
  final_u_p1 <-filter(final_u_p,grade==i)
  m=paste0(i,": ",fea)
  f_name <-paste0("./figure/01_classified_var_pie/",fea,"_",i,".pdf")
  pdf(f_name,width=7,height = 7)
  le = c("0","1","2","3","4")
  pie(final_u_p1$count,final_u_p1$ratio,main=m,col = mycolor[1:length(final_u_p1$count)])
  legend("topright",le,fill = mycolor[1:length(le)])
  dev.off()
}
lapply(c(0,"3a","3b"),p_stage)
#----------------------------------------------------------
#                  survival
library("survival")
library("survminer")
#----------------------------------------------------------
# org_grade$pfs_month_new <-as.numeric(org_grade$pfs_month)
# org_grade$OS_month <-as.numeric(org_grade$OS_month)
fit <- survfit(Surv(pfs_month_new, pro_status) ~ grade, data=org_grade)
pdf("./figure/01_survival/01_survival_pfs.pdf")
p1 <- ggsurvplot(fit,
                 legend.title="Grade",
                 title="Pfs",
                 xlab = " Time (Months)",
                 xlim = c(0, 100))
print(p1) 
dev.off()
org_grade$pfs_month <-as.numeric(org_grade$pfs_month)
fit <- survfit(Surv(pfs_month, pro_status) ~ grade, data=org_grade)
pdf("./figure/01_survival/01_survival_pfs_old.pdf")
p1 <- ggsurvplot(fit,
                 legend.title="Grade",
                 title="Pfs",
                 xlab = " Time (Months)")
print(p1) 
dev.off()
#-----------------------
fit <- survfit(Surv(os_month_new, dead) ~ grade, data=org_grade)
pdf("./figure/01_survival/01_survival_os.pdf")
p1 <- ggsurvplot(fit,
                 legend.title="Grade",
                 title="Os",
                 xlab = " Time (Months)")

# pdf("./figure/01_survival/01_survival_os.pdf")
print(p1) 
dev.off()
#--------------
org_grade$OS_month <-as.numeric(org_grade$OS_month)
fit <- survfit(Surv(OS_month, dead) ~ grade, data=org_grade)
pdf("./figure/01_survival/01_survival_os_old.pdf")
p1 <- ggsurvplot(fit,
                 legend.title="Grade",
                 title="Os",
                 xlab = " Time (Months)")

# pdf("./figure/01_survival/01_survival_os.pdf")
print(p1) 
dev.off()

#----------

aaaaa <-org_grade[,c("No","diagnosis","dead_time","dead","pfs_month","OS_month","pro_time","pro_status","pfs_month_new","os_month_new","followup")]

aaaaa$pfs_diff <-as.numeric(aaaaa$pfs_month)-aaaaa$pfs_month_new
aaaaa$os_diff <-as.numeric(aaaaa$os_month)-aaaaa$os_month_new
