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
# dat <- read.csv("01_add_age_raw_pfs_os_filter_grade_mergrPod_3A.txt",sep="\t")
load("06_3a_0_file.Rdata")

dat$age_s=NA
dat$age_s[dat$age_raw<=60]=0
dat$age_s[dat$age_raw>60]=1

dat$stage_s=NA
dat$stage_s[dat$stage<3]=0
dat$stage_s[dat$stage>=3]=1
dat$stage_s[is.na(dat$stage_s) &dat$new_pod_total==0] =0
dat$stage_s[is.na(dat$stage_s) &dat$new_pod_total==1] =1


dat$HGB_s =NA
dat$HGB_s[dat$HGB <120]=1
dat$HGB_s[dat$HGB >=120]=0
# dat$HGB_s[is.na(dat$HGB_s)]=dat$HGB0[is.na(dat$HGB_s)]
dat$HGB_s[is.na(dat$HGB_s) &dat$new_pod_total==0] =0
dat$HGB_s[is.na(dat$HGB_s) &dat$new_pod_total==1] =1

dat$LDH_s =dat$LDH_re0

dat$LN_num_s =NA
dat$LN_num_s[dat$LN_num <=4]= 0
dat$LN_num_s[dat$LN_num >4]= 1
dat$LN_num_s[is.na(dat$LN_num_s)&dat$stage_s==1]=1
dat$LN_num_s[is.na(dat$LN_num_s)&dat$stage_s==0]=0

dat$FLIPI1_count <- base::rowSums(dat[,c("age_s","stage_s","HGB_s","LDH_s","LN_num_s")], na.rm = TRUE)

# dat$diff_f1 <-dat$FLIPI1_count -dat$FLIPI1
# dat$diff_f2 <-dat$FLIPI1_count -dat$FLIPI2
#-------------------------------FLIPI1
dat$FLIPI1_count_re =dat$FLIPI1_count
dat$FLIPI1_count_re[grep("0|1",dat$FLIPI1_count)] <- "Low"
dat$FLIPI1_count_re[grep("2",dat$FLIPI1_count)] <- "Intermediate"
dat$FLIPI1_count_re[grep("3|4|5",dat$FLIPI1_count)] <- "High"

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}

fit <- survfit(Surv(pfs_month_new, pro_status) ~ FLIPI1_count_re, data=dat)
p1 <- ggsurvplot(fit,
                palette=c("#B71C1C", "#38aa34","#21618C"),
                legend.labs=c("High","Intermediate","Low"), #标签
                pval = TRUE,
                legend.title="FLIPI1",
                title="PFS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
# print(p1$plot) 

ggsave("./figure/070_pfs_FLIPI1_ori_3a_0.png",p1$plot,dpi=300,height=5.2,width=5)

fit <- survfit(Surv(os_month_new, dead) ~ FLIPI1_count_re, data=dat)
p1 <- ggsurvplot(fit,
                palette=c("#B71C1C", "#38aa34","#21618C"),
                legend.labs=c("High","Intermediate","Low"), #标签
                pval = TRUE,
                legend.title="FLIPI1",
                title="OS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
)
ggsave("./figure/070_os_FLIPI1_ori_3a_0.png",p1$plot,dpi=300,height=5.2,width=5)
#-----------------------------------

#---------------------------------------------flipi2
dat$B2mg_s <- dat$B2MG_re0
dat$LN6_s <- dat$LN6
dat$LN6_s[is.na(dat$LN6_s)&dat$new_pod_total==1]=1
dat$LN6_s[is.na(dat$LN6_s)&dat$new_pod_total==0]=0

dat$BM_s <- dat$BM
dat$BM_s[is.na(dat$BM_s)&dat$new_pod_total==1]=1
dat$BM_s[is.na(dat$BM_s)&dat$new_pod_total==0]=0
#HGB_s, age_s
dat$FLIPI2_count <- base::rowSums(dat[,c("B2mg_s","LN6_s","BM_s","HGB_s","age_s")], na.rm = TRUE)

dat$FLIPI2_count_re =dat$FLIPI2_count
dat$FLIPI2_count_re[grep("0",dat$FLIPI2_count)] <- "Low"
dat$FLIPI2_count_re[grep("1|2",dat$FLIPI2_count)] <- "Intermediate"
dat$FLIPI2_count_re[grep("3|4|5",dat$FLIPI2_count)] <- "High"
#--------------------------------------------------------

fit <- survfit(Surv(pfs_month_new, pro_status) ~ FLIPI2_count_re, data=dat)
p1 <- ggsurvplot(fit,
                palette=c("#B71C1C", "#38aa34","#21618C"),
                legend.labs=c("High","Intermediate","Low"), #标签
                pval = TRUE,
                legend.title="FLIPI2",
                title="PFS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
# print(p1$plot) 

ggsave("./figure/070_pfs_FLIPI2_ori_3a_0.png",p1$plot,dpi=300,height=5.2,width=5)

fit <- survfit(Surv(os_month_new, dead) ~ FLIPI2_count_re, data=dat)
p1 <- ggsurvplot(fit,
                palette=c("#B71C1C", "#38aa34","#21618C"),
                legend.labs=c("High","Intermediate","Low"), #标签
                pval = TRUE,
                legend.title="FLIPI2",
                title="OS",
                xlab = "Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
)
ggsave("./figure/070_os_FLIPI2_ori_3a_0.png",p1$plot,dpi=300,height=5.2,width=5)
#------------------------------------------primapi

dat$primapi_re  <- NA
dat$primapi_re[dat$B2mg>3]= "High"
dat$primapi_re[dat$B2mg<=3 &dat$BM_s==1 ]="Intermediate"
dat$primapi_re[dat$B2mg<=3 &dat$BM_s==0 ]="Low"
dat$primapi_re[is.na(dat$primapi_re) & dat$pod_total==0 ]="Low"
dat$primapi_re[is.na(dat$primapi_re) & dat$pod_total>=3 ]="High"
dat$primapi_re[is.na(dat$primapi_re)]="Intermediate"



# save(dat,file ="07_add_model_index_age_raw_pfs_os_filter_grade_mergrPod_3A_0.Rdata")

fit <- survfit(Surv(pfs_month_new, pro_status) ~ primapi_re, data=dat)
p1 <- ggsurvplot(fit,
                palette=c("#B71C1C", "#38aa34","#21618C"),
                legend.labs=c("High","Intermediate","Low"), #标签
                pval = TRUE,
                legend.title="PRIMA-PI",
                title="PFS",
                xlab = "Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
# print(p1$plot) 

ggsave("./figure/070_Primapi_pfs_ori_3a_0.png",p1$plot,dpi=300,height=5.2,width=5)
#--------------os
fit <- survfit(Surv(os_month_new, dead) ~ primapi_re, data=dat)
p1 <- ggsurvplot(fit,
                palette=c("#B71C1C", "#38aa34","#21618C"),
                legend.labs=c("High","Intermediate","Low"), #标签
                pval = TRUE,
                legend.title="PRIMA-PI",
                title="OS",
                xlab = "Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
)
ggsave("./figure/070_Primapi_os_ori_3a_0.png",p1$plot,dpi=300,height=5.2,width=5)
#----------------------------------b2mg
dat$ldh_1.5 <-NA
dat$ldh_1.5[dat$LDH >1.5*dat$LDH_upper]=1
dat$ldh_1.5[dat$LDH <=1.5*dat$LDH_upper]=0
dat$ldh_1.5[is.na(dat$ldh_1.5)&dat$new_pod_total==1]=1
dat$ldh_1.5[is.na(dat$ldh_1.5)&dat$new_pod_total==0]=0

dat$b2m_1.5 <-NA
dat$b2m_1.5[dat$B2mg >1.5*dat$B2M_upper]=1
dat$b2m_1.5[dat$B2mg <=1.5*dat$B2M_upper]=0
dat$b2m_1.5[is.na(dat$b2m_1.5)&dat$new_pod_total==1]=1
dat$b2m_1.5[is.na(dat$b2m_1.5)&dat$new_pod_total==0]=0

dat$b2mg_ldh1.5s <- base::rowSums(dat[,c("ldh_1.5","b2m_1.5")], na.rm = TRUE)

fit <- survfit(Surv(pfs_month_new, pro_status) ~ b2mg_ldh1.5s, data=dat)
p1 <- ggsurvplot(fit,
                palette=c("#B71C1C", "#38aa34","#21618C"),
                legend.labs=c("High","Intermediate","Low"), #标签
                pval = TRUE,
                legend.title="B2M_LDH",
                title="PFS",
                xlab = "Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
ggsave("./figure/070_B2M_ori_3a_0_PFS.png",p1$plot,dpi=300,height=5.2,width=5)
#------------------------
fit <- survfit(Surv(os_month_new, dead) ~ b2mg_ldh1.5s, data=dat)
# pdf("./figure/070_3a_b2mg_ldh1.5s_count_os.pdf")
p1 <- ggsurvplot(fit,
                palette=c("#B71C1C", "#38aa34","#21618C"),
                legend.labs=c("High","Intermediate","Low"), #标签
                pval = TRUE,
                legend.title="B2M_LDH",
                title="OS",
                xlab = "Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
ggsave("./figure/070_B2M_ori_3a_0_OS.png",p1$plot,dpi=300,height=5.2,width=5)


dat$pod_total_merge12_34 <-  dat$pod_total
dat$pod_total_merge12_34[grep("1|2",dat$pod_total_merge12_34)] <- "1"
dat$pod_total_merge12_34[grep("3|4",dat$pod_total_merge12_34)] <- "2"


fit <- survfit(Surv(pfs_month_new, pro_status) ~pod_total_merge12_34, data=dat)
# pdf("./figure/08_pod24_pfs_cutoff1_adjust3.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C", "#38aa34","#B71C1C"),
                legend.labs=c("Low","Intermediate","High"), #标签
                pval = TRUE,
                legend.title="POD",
                title="PFS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
ggsave("./figure/07_pod_merge12_34_pfs_ori.png",p1$plot,dpi=300,height=5.2,width=5)

fit <- survfit(Surv(os_month_new, dead) ~pod_total_merge12_34, data=dat)
# pdf("./figure/08_pod24_Os_cutoff1_adjust3.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C", "#38aa34","#B71C1C"),
                legend.labs=c("Low","Intermediate","High"), #标签
                pval = TRUE,
                legend.title="POD",
                title="OS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
ggsave("./figure/07_pod_merge12_34_os_ori.png",p1$plot,dpi=300,height=5.2,width=5)

save(dat,file ="07_add_model_index_age_raw_pfs_os_filter_grade_mergrPod_3A_0.Rdata")


#-----------------------------------------------

#-----------------------------------
