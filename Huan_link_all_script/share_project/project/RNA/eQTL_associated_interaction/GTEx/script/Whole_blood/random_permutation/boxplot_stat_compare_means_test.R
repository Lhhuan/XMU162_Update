library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/enrichment/interval_18/ALL/")
hotspot<-read.table("hotspot_15_state_count_Whole_Blood_cutoff_0.176.txt",header = T,sep = "\t") %>% as.data.frame()
random<-read.table("original_random_15_state_count_Whole_Blood_cutoff_0.176.txt",header = T,sep = "\t") %>% as.data.frame()
random_permutation<-read.table("random_permutation_15_state_count_Whole_Blood_cutoff_0.176.txt",header = T,sep = "\t") %>% as.data.frame()

hotspot$random_number <- 1
hotspot$Class <- "hotspot"
random_permutation$Class <- "random permutation"
random$Class <- "random"



rs<- bind_rows(hotspot,random_permutation,random)


unique_chromatin_state <- unique(rs$chromatin_state)
figure_list <-list()
my_comparisons <- list(c("hotspot","random"), c("hotspot","random permutation"),c("random permutation","random"))
my_comparisons2 <- list(c("hotspot","random"))


setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/random_permutation/test/")
for(state in unique_chromatin_state){
    aa <-filter(rs,chromatin_state==state)
    title_name<-state
    state <-str_replace(state,"/","_")
    # figure_list[[i]] <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ theme
    pdf(paste0("15_state_",state,".pdf")) 
    p <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ theme(legend.position ="none")+ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5))+
    # p <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ 
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "greater"))
    # i=i+1
    print(state)
    print(p)
    dev.off()
}


state = "15_Quies"

    aa <-filter(rs,chromatin_state==state)
    title_name<-state
    state <-str_replace(state,"/","_")
    # figure_list[[i]] <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ theme
    pdf(paste0("15_state_less_",state,".pdf")) 
    p <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ theme(legend.position ="none")+ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5))+
    # p <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ 
    stat_compare_means(comparisons = my_comparisons2, method = "wilcox.test", method.args = list(alternative = "less"))
    # i=i+1
    print(state)
    print(p)
    dev.off()


figure_list <-list()
i=1

for(state in unique_chromatin_state){
    aa <-filter(rs,chromatin_state==state)
    title_name<-state
    state <-str_replace(state,"/","_")
    figure_list[[i]] <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ 
    theme(legend.position ="none")+ggtitle(title_name) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "greater"))
    i=i+1
    print(state)
}

pdf("boxplot_compare_background.pdf") 
p1<-marrangeGrob(figure_list,nrow=2,ncol=2) 
print(p1)
dev.off()

#-------------
figure_list <-list()
i=1
rs1 <-filter(rs,random_number <501)
for(state in unique_chromatin_state){
    aa <-filter(rs1,chromatin_state==state)
    title_name<-state
    state <-str_replace(state,"/","_")
    figure_list[[i]] <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ 
    theme(legend.position ="none")+ggtitle(title_name) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "greater"))
    i=i+1
    print(state)
}

pdf("boxplot_compare_background_random500.pdf") 
p1<-marrangeGrob(figure_list,nrow=2,ncol=2) 
print(p1)
dev.off()


setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/enrichment/interval_18/ALL/")
hotspot<-read.table("hotspot_15_state_count_Whole_Blood_cutoff_0.176.txt",header = T,sep = "\t") %>% as.data.frame()

random_permutation<-read.table("random_permutation_15_state_count_Whole_Blood_cutoff_0.176.txt",header = T,sep = "\t") %>% as.data.frame()


rs<- bind_rows(hotspot,random_permutation)
