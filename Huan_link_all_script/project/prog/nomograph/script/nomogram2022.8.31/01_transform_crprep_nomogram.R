library(ggplot2)
library(dplyr)
library(stringr)
library(regplot)
library("survival")
rm(list=ls())

load("crprep_nomogram_f.Rdata")
library(mstate)
df.w <- crprep('Time', 're_transform',
               data=fdat, trans=c(1,2),
               cens=0, id='No',
               keep=c('Stage','Bsym','LoDLIN_6cm','BM','Spleen','Extra_sites','ECOG','B2MG','LDH','Grade'))

df.w$T<- df.w$Tstop - df.w$Tstart


m.crr<- coxph(Surv(T,status==1)~Stage+Bsym+LoDLIN_6cm+BM+Spleen+Extra_sites+ECOG+B2MG+LDH+Grade,
              data=df.w[df.w$failcode==1,],
              weight=weight.cens,
              subset=failcode==1)
summary(m.crr)


library(regplot)
regplot(m.crr,plots=c("density","spikes"), observation=FALSE,points=TRUE,dencol="#84B1ED",spkcol="#84B1ED",
        failtime = c(12, 24),title="Transform: continuous variables")
#------------------------------------------facctor

for (i in c(1:9,14)){
  fdat[,i] <- as.factor(fdat[,i])
  print(i)
}

library(mstate)
df.w <- crprep('Time', 're_transform',
                data=fdat, trans=c(1,2),
                cens=0, id='No',
                keep=c('Stage','Bsym','LoDLIN_6cm','BM','Spleen','Extra_sites','ECOG','B2MG','LDH','Grade'))

df.w$T<- df.w$Tstop - df.w$Tstart


m.crr<- coxph(Surv(T,status==1)~Stage+Bsym+LoDLIN_6cm+BM+Spleen+Extra_sites+ECOG+B2MG+LDH+Grade,
              data=df.w[df.w$failcode==1,],
              weight=weight.cens,
              subset=failcode==1)
summary(m.crr)


library(regplot)
regplot(m.crr,plots=c("density","spikes"), observation=FALSE,points=TRUE,dencol="#84B1ED",spkcol="#84B1ED",
        failtime = c(12, 24),title="Transform: classification")



