library(ggplot2)
library(dplyr)
library(stringr)
library(regplot)
library("survival")
rm(list=ls())

all <-read.csv("all_data.2022-3-21-new.csv",na.strings = "")

dat <-all[,c("No","grade","转化与否","复发与否")]
colnames(dat)[3:4] <-c("transform","relapse")
dat$No[grep("1620",dat$No)] <- 1620

dat$grade[grep("低级别",dat$grade)] <- 0
dat$grade[grep("转化",dat$grade)] <- "3b"
dat$grade[grep("3B",dat$grade)] <- "3b"
dat$grade[grep("高级别",dat$grade)] <- "3b"
dat$grade[grep("滤泡弥漫混合型",dat$grade)] <- "3b"
dat$grade[grep("滤泡弥漫型",dat$grade)] <- "3b"
dat$grade[grep("滤泡中心母细胞型",dat$grade)] <- "3b"
dat$grade[grep("弥漫大B",dat$grade)] <- "3b"
dat$grade[grep("3A",dat$grade)] <- "3a"
dat$grade[grep("3$",dat$grade)] <-"3a"
dat$grade[grep("儿童",dat$grade)] <- NA
dat$grade[grep("NA",dat$grade)] <- NA
# dat$grade[dat$grade == 1] <- 0
# dat$grade[dat$grade == 2] <- 0
dat$relapse[grep("NA",dat$relapse)] <- NA
dat$re_transform <-dat$transform
dat$re_transform[grep("0",dat$re_transform)] <- 2
dat$re_transform[grep("NA",dat$re_transform)] <- 0

dat$re_relapse <-dat$relapse
dat$re_relapse[grep("0",dat$re_relapse)] <- 2
dat$re_relapse[is.na(dat$re_relapse)] <-0
dat$re_relapse[grep("NA",dat$re_relapse)] <- 0
colnames(dat)[2]<-"ori_grade"


load("./04_add_age_raw_pfs_os_filter_grade_refine_ldh_b2m.Rdata")
dat1 <-filter(dat1,grade!="3b")

dat_filter <-filter(dat1,!is.na(pfs_month_new))

# features = c('stage','Bsym','LN6','BM','spleen','extend_num','ECOG','B2MG_re0_train','LDH_re0_train',"pfs_month_new","pro_status","grade","No","age_raw")
#,"LDH","B2mg"
features = c('stage','Bsym','LN6','BM','spleen','extend_num','ECOG','B2mg','LDH',"pfs_month_new","pro_status","grade","No","age_raw")
dat_used <-left_join(dat_filter[,features],dat,by="No")
colnames(dat_used)[c(1:10,14,15)] <-c('Stage','Bsym','LoDLIN_6cm','BM','Spleen','Extra_sites','ECOG','B2MG','LDH','Time','Age','Grade')

fdat<-filter(dat_used,Time!=0)
LDH_error <-filter(fdat,LDH >1000)
B2M_error <-filter(fdat,B2MG >100)
fdat <-anti_join(fdat,LDH_error,by="No")
fdat <-anti_join(fdat,B2M_error,by="No")
# for (i in c(1:9)){
# for (i in c(1:7)){
#   fdat[,i] <- as.factor(fdat[,i])
#   print(i)
# }
# 
# 


save(fdat,file = "crprep_nomogram_f.Rdata")


#------------------------------------------factor

dat3a <- filter(fdat,grade=="3a")
for (i in c(1:7)){
  dat3a[,i] <- as.factor(dat3a[,i])
  print(i)
}


library(mstate)
df.w <- crprep('Time', 're_relapse',
                data=dat3a, trans=c(1,2),
                cens=0, id='No',
                keep=c('Stage','Bsym','LoDLIN_6cm','BM','Spleen','Extra_sites','ECOG','B2MG','LDH','Age'))

df.w$T<- df.w$Tstop - df.w$Tstart


m.crr<- coxph(Surv(T,status==1)~Stage+Bsym+LoDLIN_6cm+BM+Spleen+Extra_sites+ECOG+B2MG+LDH+Age,
              data=df.w[df.w$failcode==1,],
              weight=weight.cens,
              subset=failcode==1)
summary(m.crr)


library(regplot)
regplot(m.crr,plots=c("density","spikes"), observation=FALSE,points=TRUE,dencol="#84B1ED",spkcol="#84B1ED",
        failtime = c(12, 24),title="Relapse: 3A")

#---------------------------------------------
dat0 <- filter(fdat,grade=="0")
# dat0 <-filter(dat0,!is.na(LDH))


for (i in c(1:7)){
  dat0[,i] <- as.factor(dat0[,i])
  print(i)
}

library(mstate)
df.w <- crprep('Time', 're_relapse',
               data=dat0, trans=c(1,2),
               cens=0, id='No',
               keep=c('Stage','Bsym','LoDLIN_6cm','BM','Spleen','Extra_sites','ECOG','B2MG','LDH','Age'))

df.w$T<- df.w$Tstop - df.w$Tstart


m.crr<- coxph(Surv(T,status==1)~Stage+Bsym+LoDLIN_6cm+BM+Spleen+Extra_sites+ECOG+B2MG+LDH+Age,
              data=df.w[df.w$failcode==1,],
              weight=weight.cens,
              subset=failcode==1)
summary(m.crr)


library(regplot)
regplot(m.crr,plots=c("density","spikes"), observation=FALSE,points=TRUE,dencol="#84B1ED",spkcol="#84B1ED",
        failtime = c(12, 24),title="Relapse: grade0-2")


# 
# 
# library(mstate)
# df.w <- crprep('Time', 're_relapse',
#                data=fdat, trans=c(1,2),
#                cens=0, id='No',
#                keep=c('Stage','Bsym','LoDLIN_6cm','BM','Spleen','Extra_sites','ECOG','B2MG','LDH','Grade'))
# 
# df.w$T<- df.w$Tstop - df.w$Tstart
# 
# 
# m.crr<- coxph(Surv(T,status==1)~Stage+Bsym+LoDLIN_6cm+BM+Spleen+Extra_sites+ECOG+B2MG+LDH+Grade,
#               data=df.w[df.w$failcode==1,],
#               weight=weight.cens,
#               subset=failcode==1)
# summary(m.crr)
# 
# 
# library(regplot)
# regplot(m.crr,plots=c("density","spikes"), observation=FALSE,points=TRUE,dencol="#84B1ED",spkcol="#84B1ED",
#         failtime = c(12, 24),title="Relapse: continuous variables")
# 
