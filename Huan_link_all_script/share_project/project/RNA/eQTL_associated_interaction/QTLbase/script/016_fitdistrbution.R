library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(data.table)
# library(parallel)
library(fitdistrplus)
setwd("~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL")


org <- read.table("NHPoisson_emplambda_interval_16cutoff_7.3_all_eQTL.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org2<-filter(org,emplambda != "NA")
fitdist(org2$emplambda,"gamma")
aaa<-head(org2,10000)
bbb <-aaa$emplambda
ccc <-org2$emplambda
descdist(ccc)
o <-fitdist(bbb,"gamma",lower=c(0,0), start=list(scale=1,shape=1))
descdist(bbb)
o<-
fitdist(bbb,"gamma",method= "mse")
fitdist(ccc,"exp")


descdistr(bbb,"gamma")


plotdist(bbb,histo=TRUE,demp=TRUE)

library(MASS)

fitdistr(bbb,"gamma")
#---------------------------

library(car)

qqPlot(bbb,dist="gamma",estimate.params = TRUE, add.line = TRUE)
car::qqPlot(x, dist = "gamma", shape="mle")

qqPlot(x, dist="chisq", df=2)

library(dgof)

#--------------------------------
EnvStats::qqPlot(bbb,dist = "gamma", estimate.params = TRUE, digits = 2, points.col = "blue",   add.line = TRUE)
library(qualityTools)
qualityTools::qqPlot(bbb,"gamma")

x<-bbb

 x[which(x==0)]=x[which(x==0)]+1e-12


 qualityTools::qqPlot(x,"gamma3",main="huan")
 bbb<-bbb[which(bbb!=0)]


 EnvStats::qqPlot(bbb,dist = "gamma", estimate.params = TRUE)



result <-ks.test(x,"pgamma", 3, 2)




tmp<-data.frame(Pvalue=result$p.value,D = result$statistic,j=1)

write.table(tmp,"1234.txt",row.names = F, col.names = T,quote =F,sep="\t")

o <-fitdist(x,"gamma",lower=c(0,0), start=list(scale=1,shape=1))



qqplot(x = qgamma(ppoints(bbb), 
                 shape = 1, scale = 1),
       y = bbb,
       
       main = "QQ-Plot: NZ River Length, Gamma Distribution",
       xlab = "Theoretical Quantiles, Gamma Distribution",
       ylab = "Sample Quantiles, River Length")
abline(a = 0, b = 1, col = "dodgerblue", lwd = 2)
grid()


 xlim = c(0, 1), ylim = c(0, 1),