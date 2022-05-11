library("survival")
library("survminer")
library(randomForest)
library(stringr)
library(dplyr)
library(Hmisc)

setwd("/home/huanhuan/project/prog/output/")
load("../data/01_add_age_raw_pfs_os_filter_grade.Rdata")
dat$new_pod_total <-dat$pod_total
dat$new_pod_total[grep("1|2",dat$new_pod_total)] <- 0
dat$new_pod_total[grep("3|4",dat$new_pod_total)] <- 1
dat$SPD <-as.numeric(dat$SPD)
p1 <-lapply(1:nrow(dat),function(i){
  if(dat[i,"age_raw"]>60){Var=1}else{Var=0}
  return(Var)
})
dat$age_raw_60 <-unlist(p1)

p1 <-lapply(1:nrow(dat),function(i){
  x=dat[i,"Ki.67"]
  if(is.na(x)){Var=NA}else if(x>20){Var=1}else{Var=0}
  return(Var)
})
dat$ki67_20 <-unlist(p1)
#-------------------------
p1 <-lapply(1:nrow(dat),function(i){
  x=dat[i,"LN_num"]
  if(is.na(x)){Var=NA}else if(x>=6){Var=1}else{Var=0}
  return(Var)
})
dat$LN_num_6 <-unlist(p1)
#--------------------
p1 <-lapply(1:nrow(dat),function(i){
  x=dat[i,"extend_num"]
  if(is.na(x)){Var=NA}else if(x>0){Var=1}else{Var=0}
  return(Var)
})
dat$extend_num_0 <-unlist(p1)
#--------------------------
p1 <-lapply(1:nrow(dat),function(i){
  x=dat[i,"stage"]
  if(is.na(x)){Var=NA}else if(x>2){Var=1}else{Var=0}
  return(Var)
})
dat$stage_III_IV <-unlist(p1)
#-------------------------
p1 <-lapply(1:nrow(dat),function(i){
  x=dat[i,"stage"]
  if(is.na(x)){Var=NA}else if(x<3){Var=1}else{Var=0}
  return(Var)
})
dat$stage_I_II <-unlist(p1)
#---------------------------------
p1 <-lapply(1:nrow(dat),function(i){
  x=dat[i,"SUVmax"]
  if(is.na(x)){Var=NA}else if(x>2){Var=1}else{Var=0}
  return(Var)
})
dat$SUVmax_2 <-unlist(p1)
#-----------------------------------
p1 <-lapply(1:nrow(dat),function(i){
  x=dat[i,"LDH"]
  if(is.na(x)){Var=NA}else if(x>300){Var=1}else{Var=0}
  return(Var)
})
dat$LDH_300 <-unlist(p1)
#------------------------------
p1 <-lapply(1:nrow(dat),function(i){
  x=dat[i,"Lym"]
  y=dat[i,"Mono"]
  if(is.na(x)| is.na(y)){Var=NA}else if(x/y>10){Var=1}else{Var=0}
  return(Var)
})
dat$Lym_Mono_10 <-unlist(p1)
#_----------------------------------------
p1 <-lapply(1:nrow(dat),function(i){
  x=dat[i,"B2mg"]
  if(is.na(x)){Var=NA}else if(x>3.4){Var=1}else{Var=0}
  return(Var)
})
dat$B2mg_3.4 <-unlist(p1)
#-----------------------------------------
p1 <-lapply(1:nrow(dat),function(i){
  x=dat[i,"SPD"]
  if(is.na(x)){Var=NA}else if(x>0){Var=1}else{Var=0}
  return(Var)
})
dat$SPD_0 <-unlist(p1)
#---------------------------
save(dat,file ="01_add_age_raw_pfs_os_filter_grade_mergrPod.Rdata")

prog <-dat
# prog <-dat 
#----------------------------------------
# load("/home/huanhuan/project/prog/data/dat_rawfilter.Rdata")

# dat <-dat%>%select("NO","Ki67")
# colnames(dat)[1] <-"No"
# dat$No <-as.character(dat$No)
# dat_m <-inner_join(prog,dat,by="No")

# datk <-dat_m%>% select("No","Ki.67","Ki67")
# datk$diff <-abs(datk$Ki.67 -as.numeric(datk$Ki67))


fill <- function(x,y){ #x: ref, y:na
  re <- lapply(which(is.na(y)),function(i){
    print(i)
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    yy <- y[!is.na(y)]
    a <- quantile(yy,range01(rank(x,na.last = TRUE))[i])
    return(a)
  })
}
prog$SUVmax[is.na(prog$SUVmax)] <- unlist(fill(prog$Ki.67,prog$SUVmax))
prog$Ki.67[is.na(prog$Ki.67)] <- unlist(fill(prog$SUVmax,prog$Ki.67))
prog$SPD[is.na(prog$SPD)] <- unlist(fill(prog$SUVmax,prog$SPD))
#---------
prog$B2mg[is.na(prog$B2mg)] <- unlist(fill(prog$SUVmax,prog$B2mg))
prog$LDH[is.na(prog$LDH)] <- unlist(fill(prog$SUVmax,prog$LDH))

prog[,c("gender","stage","Bsym","LN_num","site0","extend","BM","spleen","extend_num","BM_extend","LN6","X150b2mg_ldh","b2mg_LDH","ECOG","LDH0","WBC","HGB","HGB0","PLT","Mono","Lym","interm_res","end_res","CR","Rmaintain","progress","trans","relapes","relapse_res","age_raw")] <-na.roughfix(prog[,c("gender","stage","Bsym","LN_num","site0","extend","BM","spleen","extend_num","BM_extend","LN6","X150b2mg_ldh","b2mg_LDH","ECOG","LDH0","WBC","HGB","HGB0","PLT","Mono","Lym","interm_res","end_res","CR","Rmaintain","progress","trans","relapes","relapse_res","age_raw")])
#-------------------
# prog[,c("gender","stage","Bsym","LN_num","site0","extend","BM","spleen","extend_num","BM_extend","LN6","X150b2mg_ldh","b2mg_LDH","ECOG","LDH0","WBC","HGB","HGB0","PLT","Mono","Lym","interm_res","end_res","CR","Rmaintain","progress","trans","relapes","relapse_res","age_raw","age_raw_60","ki67_20","LN_num_6","extend_num_0","stage_III_IV","stage_I_II","SUVmax_2","LDH_300","Lym_Mono_10","B2mg_3.4","SPD_0")] <-na.roughfix(prog[,c("gender","stage","Bsym","LN_num","site0","extend","BM","spleen","extend_num","BM_extend","LN6","X150b2mg_ldh","b2mg_LDH","ECOG","LDH0","WBC","HGB","HGB0","PLT","Mono","Lym","interm_res","end_res","CR","Rmaintain","progress","trans","relapes","relapse_res","age_raw","age_raw_60","ki67_20","LN_num_6","extend_num_0","stage_III_IV","stage_I_II","SUVmax_2","LDH_300","Lym_Mono_10","B2mg_3.4","SPD_0")])


# covariates <- c("gender","Ki.67","stage","Bsym","LN_num","site0","extend","BM","spleen","extend_num","BM_extend","LN6","SUVmax","SPD","X150b2mg_ldh","b2mg_LDH","ECOG","B2mg","LDH","LDH0","WBC","HGB","HGB0","PLT","Mono","Lym","interm_res","end_res","CR","Rmaintain","progress","trans","relapes","relapse_res","age_raw","age_raw_60","ki67_20","LN_num_6","extend_num_0","stage_III_IV","stage_I_II","SUVmax_2","LDH_300","Lym_Mono_10","B2mg_3.4","SPD_0")
covariates <- c("gender","Ki.67","stage","Bsym","LN_num","site0","extend","BM","spleen","extend_num","BM_extend","LN6","SUVmax","SPD","X150b2mg_ldh","b2mg_LDH","ECOG","B2mg","LDH","LDH0","WBC","HGB","HGB0","PLT","Mono","Lym","interm_res","end_res","CR","Rmaintain","progress","trans","relapes","relapse_res","age_raw")



univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(os_month_new, dead)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = prog)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         HR <-signif(x$coef[2], digits=2);
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR_re <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR_re,HR,HR.confint.lower,HR.confint.upper)
                         names(res)<-c("Pvalue","HR (95% CI for HR)","HR","HR_confint_lower","HR_confint_upper")
                         return(res)
                       })

res <- t(as.data.frame(univ_results, check.names = FALSE))
result1 <- as.data.frame(res)
result1[,c(1,3:5)] <-apply(result1[,c(1,3:5)],2,as.numeric)
# result1$var <- rownames(result1)
aa <-filter(result1,Pvalue <0.05)
aa$var <-rownames(aa)

write.table(aa, "03_fill_univariate_cox_result_fill_all.txt",quote=F,sep="\t",row.names=FALSE)

# my_comparisons
# prog$pod_total <-as.factor(prog$pod_total)
# p1 <-ggplot(prog,aes(x=pod_total,y=Ki.67))+geom_boxplot(aes(fill=pod_total),width=0.3,outlier.colour = NA)
# print(p1)
# dev.off()
#-----------------------------------------------plot
pstar <-function(i){
  Pvalue <-i
  if(Pvalue <= 0.0001){
    pval_star = "****"
  }else if(Pvalue <= 0.001){
    pval_star = "***"
  }else if(Pvalue <= 0.01){
    pval_star = "**"
  }else if(Pvalue <= 0.05){
    pval_star = "*"
  }else{
    pval_star = "ns"
  }
  return(pval_star)
}

pp <-lapply(aa[c(1:nrow(aa)),1],pstar)
aa$significant <-unlist(pp)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
aa <-aa%>%arrange(HR) 
aa$var <-as.factor(aa$var)
aa$var <-gsub("_"," ",aa$var)
aa$var <-capitalize(aa$var)
pdf("./figure/03_univariate_cox_fill_all.pdf",width=5,height=5)
p1<-ggplot(data=aa, aes(x=reorder(var,X=HR),y=HR, ymin=HR_confint_lower , ymax=HR_confint_upper))
  p1<-p1+geom_pointrange(size = 0.3) 
  p1<-p1+coord_flip()# +  # flip coordinates (puts labels on y axis)
  p1<-p1+labs(y="Hazard ratio",x="",title="Univariable analysis") 
  p1 <-p1+p_theme+theme(axis.text.x = element_text(color = "black"),axis.text.y = element_text(color = "black"),plot.title = element_text(hjust = 0.5))
  p1 <-p1+geom_text(aes(x=var,y=7.2,label = significant)) +
      geom_hline(yintercept=1,colour="#c00024",linetype=3)
print(p1)
dev.off()

#------------------  
# res.cox <- coxph(Surv(pod.total) ~ gender + ECOG + stage +Bsymptom + Lymph.node+ spleen +BM +transform+ PRIME.PI+regimen.R. + interim.time+ interim.evaluation +
# end.evaluation + R.course+ progress + HGB.g.l.+ MONO + SPD, data =  prog_roughfix)
# res.cox <- coxph(Surv(os_month_new, dead) ~ gender+Ki.67+stage+Bsym+LN_num+site0+extend+BM+spleen+extend_num+BM_extend+LN6+SUVmax+SPD+X150b2mg_ldh+b2mg_LDH+ECOG+B2mg+LDH+LDH0+WBC+HGB+HGB0+PLT+Mono+Lym+interm_res+end_res+CR+Rmaintain+progress+trans+relapes+relapse_res+age_raw, data =  prog)
res.cox <- coxph(Surv(os_month_new, dead) ~ Rmaintain+CR+HGB+LDH+age_raw+extend_num+stage+site0+HGB0+end_res+Bsym+ECOG+interm_res+LN6+X150b2mg_ldh+b2mg_LDH+progress+trans+LDH0, data =  prog)
x <- summary(res.cox)
pvalue=signif(as.matrix(x$coefficients)[,5],2)
HR=signif(as.matrix(x$coefficients)[,2],2)
low=signif(x$conf.int[,3],2)
high=signif(x$conf.int[,4],2)
multi_res=data.frame(p.value=pvalue,
                     HR_re=paste(HR," (",low,"-",high,")",sep=""),
                     HR =HR,
                     low=low,
                     high=high,
                     stringsAsFactors = F
)
multi_res$value <-rownames(multi_res)
aa <-filter(multi_res,p.value<0.05)
aa$var <-rownames(aa)
pp <-lapply(aa[c(1:nrow(aa)),1],pstar)
aa$significant <-unlist(pp)
aa <-aa%>%arrange(HR)
aa$var <-as.factor(aa$var)
aa$var <-gsub("_"," ",aa$var)
aa$var <-capitalize(aa$var)
pdf("./figure/03_Multivariable_cox_fill_all.pdf",width=5,height=5)
p1<-ggplot(data=aa, aes(x=reorder(var,X=HR),y=HR, ymin=low , ymax=high))
  p1<-p1+geom_pointrange(size = 0.3) 
  p1<-p1+coord_flip()# +  # flip coordinates (puts labels on y axis)
  p1<-p1+labs(y="Hazard ratio",x="",title="Multivariable analysis") 
  p1 <-p1+p_theme+theme(axis.text.x = element_text(color = "black"),axis.text.y = element_text(color = "black"),plot.title = element_text(hjust = 0.5))
  p1 <-p1+geom_text(aes(x=var,y=5,label = significant)) +
      geom_hline(yintercept=1,colour="#c00024",linetype=3)
print(p1)
dev.off()


write.table(file="03_multivariate_cox_result_fill_all.txt",multi_res,quote=F,sep="\t")

# write.table(file="multivariate_cox_result_sig.txt",multi_res_sig,quote=F,sep="\t")