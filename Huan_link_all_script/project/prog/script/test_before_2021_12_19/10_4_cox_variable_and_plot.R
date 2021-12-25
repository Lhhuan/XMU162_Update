library("survival")
library("survminer")
library(randomForest)
library(stringr)
library(dplyr)
library(Hmisc)

setwd("/home/huanhuan/project/prog/output/")
load("07_add_model_index_age_raw_pfs_os_filter_grade_mergrPod_3A.Rdata")
covariates <- c('Ki.67','stage','Bsym','LN_num','site0','BM','spleen','extend_num','BM_extend','SUVmax','SPD','ECOG','B2mg','LDH','HGB','age_raw','Lym_Mono')
# X = dataset.loc[:, ['B2mg','LN_num','LDH','age_raw','Lym_Mono','HGB','Ki.67','SPD','SUVmax','Bsym','BM']]
# c('gender','Ki.67','stage','Bsym','LN_num','site0','extend','BM','spleen','extend_num','BM_extend','LN6','SUVmax','SPD','X150b2mg_ldh','b2mg_LDH','ECOG','B2mg','LDH','LDH0','HGB','HGB0','Mono','Lym','age_raw','age_raw_60','ki67_20','LN_num_6','extend_num_0','SUVmax_2','LDH_300','Lym_Mono','B2mg_3.4','SPD_0')

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(os_month_new, dead)~', x)))
                        # function(x) as.formula(paste('Surv(pfs_month_new, pro_status)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = dat)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         HR <-signif(x$coef[2], digits=4);
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 4)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],4)
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
# aa <-filter(result1,Pvalue <0.05)
aa<-result1
aa$var <-rownames(aa)
write.table(aa,"10_4_fill_univariate_cox_result_not_fill_3a_os.txt" ,quote=F,sep="\t",row.names=FALSE)
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

 # ,axis.line = element_line(colour = "black")

aa$var <-capitalize(aa$var)
# aa$var <-gsub("Ki.67","Ki.67 >20%",aa$var)
aa$var <-gsub("LN_num","Lymph nodes",aa$var)
aa$var <-gsub("Extend_num","Extra sites",aa$var)
aa$var <-gsub("BM_extend","Extra BM",aa$var)
aa$var <-gsub("Lym_Mono","Lym/Mono",aa$var)
aa$var <-gsub("B2mg","B2M",aa$var)
aa$var <-gsub("Age_raw","Age",aa$var)

aa <-aa%>%arrange(HR) 
aa$var <-factor(aa$var,levels=aa$var)
#----------
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), 
                                                panel.grid.minor = element_blank(), 
                                                panel.background = element_rect(color="black",size=1.2), 
                                                axis.title.x = element_text(size =10),
                                                axis.text.x = element_text(color="black",size=9),
                                                axis.text.y = element_text(color="black",size = 10),
                                                axis.ticks.y=element_blank())

p1<-ggplot(data=aa, aes(x=HR,y=var))+
  geom_vline(xintercept=1,colour="#12497f",linetype=2,size=0.6)+
  geom_errorbarh(aes(xmax =HR_confint_upper, xmin = HR_confint_lower), height = 0.3)+
  geom_point(size=2,shape=22,fill="#c00024",colour="#c00024",)+
  theme_bw()+
  p_theme +
  labs(x="Hazard ratio",y="")+
  scale_x_continuous(limits= c(0, 6), breaks= seq(0,6,1))

p2 <-ggplot(data=aa, aes(x=-log10(Pvalue),y=var))+
    geom_bar(stat = "identity", width = 0.3,fill="#12497f") +
    geom_vline(xintercept=-log10(0.05),colour="#c00024",linetype=2,size=0.6)+
    theme_bw()+
    p_theme+
    labs(x="-log10(P-value)",y="")+
    scale_x_continuous(limits= c(0, 7.2), breaks= seq(0,7,1))+
    theme(axis.text.y=element_blank(),
    axis.title.y = element_blank()) 
    # +margin(2, 2, 2, 2, "cm")

library(patchwork)
p3<- (p1+theme(plot.margin=unit(c(0,0,0,0),"cm")))+(p2+theme(plot.margin=unit(c(0,0,0,0),"cm"))) + plot_layout(nrow = 1, width = c(100,30))

# p3 <-p1+plot_spacer()+ p2 +plot_layout(nrow = 1, width = c(100,0,40))
ggsave("./figure/10_4_univariate_cox_not_fill_na_3a_os_forest.png",p3,dpi=300,width=6.5,height=5)


#-----------------------------------------------pfs


univ_formulas <- sapply(covariates,
                        # function(x) as.formula(paste('Surv(os_month_new, dead)~', x)))
                        function(x) as.formula(paste('Surv(pfs_month_new, pro_status)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = dat)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         HR <-signif(x$coef[2], digits=4);
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 4)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],4)
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
# aa <-filter(result1,Pvalue <0.05)
aa<-result1
aa$var <-rownames(aa)
write.table(aa,"10_4_fill_univariate_cox_result_not_fill_3a_pfs.txt" ,quote=F,sep="\t",row.names=FALSE)
#-----------------------------------------------plot

pp <-lapply(aa[c(1:nrow(aa)),1],pstar)
aa$significant <-unlist(pp)

 # ,axis.line = element_line(colour = "black")

aa$var <-capitalize(aa$var)
# aa$var <-gsub("Ki.67","Ki.67 >20%",aa$var)
aa$var <-gsub("LN_num","Lymph nodes",aa$var)
aa$var <-gsub("Extend_num","Extra sites",aa$var)
aa$var <-gsub("BM_extend","Extra BM",aa$var)
aa$var <-gsub("Lym_Mono","Lym/Mono",aa$var)
aa$var <-gsub("B2mg","B2M",aa$var)
aa$var <-gsub("Age_raw","Age",aa$var)

aa <-aa%>%arrange(HR) 
aa$var <-factor(aa$var,levels=aa$var)
#----------

p1<-ggplot(data=aa, aes(x=HR,y=var))+
  geom_vline(xintercept=1,colour="#12497f",linetype=2,size=0.6)+
  geom_errorbarh(aes(xmax =HR_confint_upper, xmin = HR_confint_lower), height = 0.3)+
  geom_point(size=2,shape=22,fill="#c00024",colour="#c00024",)+
  theme_bw()+
  p_theme +
  labs(x="Hazard ratio",y="")+
  scale_x_continuous(limits= c(0, 6), breaks= seq(0,6,1))

p2 <-ggplot(data=aa, aes(x=-log10(Pvalue),y=var))+
    geom_bar(stat = "identity", width = 0.3,fill="#12497f") +
    geom_vline(xintercept=-log10(0.05),colour="#c00024",linetype=2,size=0.6)+
    theme_bw()+
    p_theme+
    labs(x="-log10(P-value)",y="")+
    scale_x_continuous(limits= c(0, 7.2), breaks= seq(0,7,1))+
    theme(axis.text.y=element_blank(),
    axis.title.y = element_blank()) 
    # +margin(2, 2, 2, 2, "cm")

library(patchwork)
p3<- (p1+theme(plot.margin=unit(c(0,0,0,0),"cm")))+(p2+theme(plot.margin=unit(c(0,0,0,0),"cm"))) + plot_layout(nrow = 1, width = c(100,30))
# p3 <-p1+plot_spacer()+ p2 +plot_layout(nrow = 1, width = c(100,0,40))
ggsave("./figure/10_4_univariate_cox_not_fill_na_3a_pfs_forest.png",p3,dpi=300,width=6.5,height=5)