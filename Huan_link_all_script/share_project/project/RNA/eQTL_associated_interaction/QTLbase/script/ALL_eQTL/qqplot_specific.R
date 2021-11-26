library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(data.table)
library(parallel)
# library(fitdistrplus)
library(qualityTools)
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/ALL_eQTL/")

ProcessBedGz <- function(j = NULL){
     file_name<-paste("../../output/ALL_eQTL/NHPoisson_emplambda_interval_",j,"cutoff_7.3_all_eQTL.txt.gz",sep="")
     org <- read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
     org2<-filter(org,emplambda != "NA")

     emplambda<-org2$emplambda

     emplambda[which(emplambda==0)]=emplambda[which(emplambda==0)]+1e-12
#----------------------------
      for (alpha in c(1:3)){
            for(beta in c(1:3)){
                  title =paste("Q-Q Plot for gamma distribution: interval",j,alpha,beta,sep =" " )
                  figure_name1<-paste("./add_qqplot/alpha",alpha,"/Q_Q_Plot_for_gamma_distribution_interval_",j,"_",alpha,"_",beta,"_add0",".png",sep ="")
                  png(figure_name1)
                  qqplot(x = qgamma(ppoints(org2$emplambda), 
                              shape = alpha, scale = beta),
                        y = org2$emplambda,
                        xlim = c(0, 1), ylim = c(0, 1),
                        main = title,
                        xlab = "emplambda, Gamma Distribution",
                        ylab = "emplambda")
                  abline(a = 0, b = 1, col = "dodgerblue", lwd = 2)
                  grid()
                  dev.off()
            }
      }
}



# # file2<-files[1:2]
# interval <-c(10:17)
# # interval <-c(18,20,24,28,32)
# test_result <-mclapply(interval, ProcessBedGz, mc.cores = 10)

interval <-c(12,18,20,24,28,32)
# interval <-c(10:17)
mclapply(interval, ProcessBedGz, mc.cores = 20)




#----------------------------
j=6
file_name<-paste("../../output/ALL_eQTL/NHPoisson_emplambda_interval_",j,"cutoff_7.3_all_eQTL.txt.gz",sep="")
org <- read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
org2<-filter(org,emplambda != "NA")


# calculating sample moments
empl_samp_moment_1 = mean(org2$emplambda)
empl_samp_moment_2 = mean(org2$emplambda ^ 2)

# method of moments estimators
len_alpha_mom = empl_samp_moment_1 ^ 2 / (empl_samp_moment_2 - empl_samp_moment_1 ^ 2)
len_beta_mom  = empl_samp_moment_1 / len_alpha_mom

# estimates for this dataset
c(len_alpha_mom, len_beta_mom)


hist(org2$emplambda, probability = TRUE, 
     main = "NZ River Lengths", xlab = "River Length (km)", 
     breaks = 12)
box()
grid()
curve(dgamma(x, shape = len_alpha_mom, scale = len_beta_mom), 
      add = TRUE, col = "darkorange", lwd = 2)

png("1234.png")
qqplot(x = qgamma(ppoints(org2$emplambda), 
                 shape = len_alpha_mom, scale = len_beta_mom),
       y = org2$emplambda,
       xlim = c(0, 1), ylim = c(0, 1),
       main = "QQ-Plot: emplambda, Gamma Distribution",
       xlab = "emplambda, Gamma Distribution",
       ylab = "emplambda")
abline(a = 0, b = 1, col = "dodgerblue", lwd = 2)
grid()
dev.off()
