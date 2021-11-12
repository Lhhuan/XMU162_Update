library("survival")
library("survminer")
library(randomForest)
library(stringr)
library(dplyr)

setwd("/home/huanhuan/project/prog/output/")
load("../data/01_add_age_raw_pfs_os_filter_grade.Rdata")




############################
# load("step1.Rdata")
prog <- dat
prog[,c("ECOG","regimen.R.","interim.time","interim.evaluation","end.evaluation","R.course","progress")] <-apply(prog[,c("ECOG","regimen.R.","interim.time","interim.evaluation","end.evaluation","R.course","progress")],2,as.integer)%>%as.data.frame()


prog$Ki67 <-as.numeric(prog$Ki67)

fill <- function(x,y){ #x: ref, y:na
  re <- lapply(which(is.na(y)),function(i){
    print(i)
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    yy <- y[!is.na(y)]
    a <- quantile(yy,range01(rank(x,na.last = TRUE))[i])
    return(a)
  })
}

# dat$SUVmax[is.na(dat$SUVmax)] <- unlist(fill(dat$Ki67,dat$SUVmax))  
prog$SUVmax[is.na(prog$SUVmax)] <- unlist(fill(prog$Ki67,prog$SUVmax))
prog$Ki67[is.na(prog$Ki67)] <- unlist(fill(prog$SUVmax,prog$Ki67))
prog$SPD[is.na(prog$SPD)] <- unlist(fill(prog$SUVmax,prog$SPD))
#---------
prog$B2.MG.ng.L.[is.na(prog$B2.MG.ng.L.)] <- unlist(fill(prog$SUVmax,prog$B2.MG.ng.L.))
prog$LDH.U.L.[is.na(prog$LDH.U.L.)] <- unlist(fill(prog$SUVmax,prog$LDH.U.L.))

fff_c <-c("gender","age","ECOG","stage","Bsymptom","Lymph.node","extranodal","spleen","BM","transform","PRIME.PI","Radiotherpy","regimen.R.","interim.time","interim.evaluation","end.evaluation","R.course","progress","grade.classify")

prog[,fff_c] <-apply(prog[,fff_c],2,as.factor)%>%as.data.frame()

prog_roughfix <-na.roughfix(prog)

prog_roughfix[,fff_c] <-apply(prog_roughfix[,fff_c],2,as.character)%>%as.data.frame()
prog_roughfix[,fff_c] <-apply(prog_roughfix[,fff_c],2,as.integer)%>%as.data.frame()


save(prog_roughfix,file="01_relate_fill.Rdata")

# res.cox <- coxph(Surv(pod.total) ~ gender, data = prog_roughfix)
# colnames()
# res.cox <- coxph(Surv(time, status) ~ sex, data = lung)
# res.cox

# covariates <- c("gender","age","ECOG","stage","Bsymptom","Lymph.node","extranodal","spleen","BM","transform","PRIME.PI","Radiotherpy","regimen.R.","interim.time","interim.evaluation","end.evaluation","R.course","progress","grade.classify","WBC","HGB.g.l.","PLT","MONO","LYMPH","M.L","LDH.U.L.","B2.MG.ng.L.","Ki67","SUVmax","SPD")
covariates <- c("gender","age","ECOG","stage","Bsymptom","Lymph.node","extranodal","spleen","BM","transform","PRIME.PI","Radiotherpy","regimen.R.","interim.time","interim.evaluation","end.evaluation","progress","grade.classify","WBC","HGB.g.l.","PLT","MONO","LYMPH","M.L","LDH.U.L.","B2.MG.ng.L.","Ki67","SUVmax","SPD")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(pod.total)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = prog_roughfix)})


univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #��ȡpֵ
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         #��ȡHR
                         HR <-signif(x$coef[2], digits=2);
                         #��ȡ95%��������
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })

res <- t(as.data.frame(univ_results, check.names = FALSE))
result1 <- as.data.frame(res)
result1$var <- rownames(result1)
write.table(file="univariate_cox_result.txt",result1 ,quote=F,sep="\t",row.names=FALSE)




#------------------  
# res.cox <- coxph(Surv(pod.total) ~ gender + ECOG + stage +Bsymptom + Lymph.node+ spleen +BM +transform+ PRIME.PI+regimen.R. + interim.time+ interim.evaluation +
# end.evaluation + R.course+ progress + HGB.g.l.+ MONO + SPD, data =  prog_roughfix)
res.cox <- coxph(Surv(pod.total) ~ gender + ECOG + stage +Bsymptom + Lymph.node+ spleen +BM +transform+ PRIME.PI+regimen.R. + interim.time+ interim.evaluation +
                   end.evaluation + progress + HGB.g.l.+ MONO + SPD, data =  prog_roughfix)

x <- summary(res.cox)
pvalue=signif(as.matrix(x$coefficients)[,5],2)
HR=signif(as.matrix(x$coefficients)[,2],2)
low=signif(x$conf.int[,3],2)
high=signif(x$conf.int[,4],2)
multi_res=data.frame(p.value=pvalue,
                     HR=paste(HR," (",low,"-",high,")",sep=""),
                     stringsAsFactors = F
)
multi_res$value <-rownames(multi_res)
multi_res_sig <-filter(multi_res,p.value<0.05)
write.table(file="multivariate_cox_result.txt",multi_res,quote=F,sep="\t")

write.table(file="multivariate_cox_result_sig.txt",multi_res_sig,quote=F,sep="\t")