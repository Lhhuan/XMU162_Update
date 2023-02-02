library(ggplot2)
library(dplyr)
library(stringr)
library(regplot)
library("survival")
rm(list=ls())

load("crprep_nomogram_f.Rdata")

#------------------------------------------factor

dat3a <- filter(fdat,grade=="3a")
Extra_sites_error <-filter(dat3a,Extra_sites==5|Extra_sites==6)
dat3a <-anti_join(dat3a,Extra_sites_error,by="No")

for (i in c(1:7)){
  dat3a[,i] <- as.factor(dat3a[,i])
  print(i)
}
# dat3a <-dat3a[!is.na(dat3a$B2MG),]

library(mstate)
df.w <- crprep('Time', 're_transform',
               data=dat3a, trans=c(1,2),
               cens=0, id='No',
               keep=c('Stage','Bsym','LoDLIN_6cm','BM','Spleen','Extra_sites','ECOG','B2MG','LDH','Age')) #,'B2MG'

df.w$T<- df.w$Tstop - df.w$Tstart


m.crr<- coxph(Surv(T,status==1)~Stage+Bsym+LoDLIN_6cm+BM+Spleen+Extra_sites+ECOG+B2MG+LDH+Age,
              data=df.w[df.w$failcode==1,],
              weight=weight.cens,
              subset=failcode==1)
summary(m.crr)
a <-m.crr

library(regplot)
regplot(m.crr,plots=c("density","spikes"), observation=FALSE,points=TRUE,dencol="#84B1ED",spkcol="#84B1ED",
        failtime = c(12, 24),title="Transform: 3A")

#---------------------------------------------
dat0 <- filter(fdat,grade=="0")
for (i in c(1:7)){
  dat0[,i] <- as.factor(dat0[,i])
  print(i)
}

library(mstate)
df.w <- crprep('Time', 're_transform',
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
        failtime = c(12, 24),title="Transform: grade0-2")
