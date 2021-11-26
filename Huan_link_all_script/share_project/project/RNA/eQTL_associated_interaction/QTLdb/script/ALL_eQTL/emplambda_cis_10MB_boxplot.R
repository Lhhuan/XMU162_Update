library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gcookbook)
library(data.table)

cutoff =7.3
intervals =1000
j=intervals 

xQTL = "eQTL"


current_directory <- paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/ALL_",xQTL,"/", sep = "")
setwd(current_directory)
input_file <-paste("../../output/ALL_", xQTL, "/NHPoisson_emplambda_interval_",j,"_cutoff_",cutoff,"_all_",xQTL,"_cis_trans.txt.gz", sep = "")
org<-fread(input_file,header = T,sep = "\t") %>% as.data.frame()
#-------------------------------------
figure_name <-paste("boxpolt_emplambda_distance_p_cutoff_",cutoff,"_interval_",j,"_all_",xQTL,".pdf", sep = "")
pdf(figure_name,height = 5,width = 5) 
p <- ggplot(org,aes(x=cis_trans_10MB,y=emplambda))+geom_violin()+geom_boxplot(aes(fill=cis_trans_10MB),width=0.3,outlier.colour = NA)+ #plot by legend
  scale_y_continuous(name = "Emplambda") +scale_x_discrete(name = "Cis or trans")+guides(fill=FALSE)
my_comparisons <- list(c("cis","tran") ) 
# p+coord_cartesian(ylim = c(-0.5,0.85))+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", method.args = list(alternative = "less"),
p+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
##组内进行比较，指定p bar的长度并指定坐标轴p值标的位置,#+guides(fill=FALSE)是移除图例
print(p)
dev.off()





#---------------------------------------------------------
QTLs = c("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL","QTL")
n<-length(QTLs)

cutoff =7.3
intervals =1000
j=intervals 

fun=function(i){
    xQTL =QTLs[i]
    print(i)
    print(xQTL)
    current_directory <- paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/ALL_",xQTL,"/", sep = "")
    setwd(current_directory)
    input_file <-paste("../../output/ALL_", xQTL, "/NHPoisson_emplambda_interval_",j,"_cutoff_",cutoff,"_all_",xQTL,"_cis_trans.txt.gz", sep = "")
    org<-fread(input_file,header = T,sep = "\t") %>% as.data.frame()
    #-------------------------------------
    figure_name <-paste("boxpolt_emplambda_distance_cutoff_",cutoff,"_interval_",j,"_all_",xQTL,".pdf", sep = "")
    pdf(figure_name,height = 5,width = 5) 
    p <- ggplot(org,aes(x=cis_trans_10MB,y=emplambda))+geom_boxplot(aes(fill=cis_trans_10MB),width=0.3,outlier.colour = NA)+ #plot by legend
      scale_y_continuous(name = "Emplambda") +scale_x_discrete(name = "Cis or trans")+guides(fill=FALSE)
    my_comparisons <- list(c("cis","tran") ) 
    # p+coord_cartesian(ylim = c(-0.5,0.85))+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", method.args = list(alternative = "less"),
    p+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")+guides(fill=FALSE)
    ##组内进行比较，指定p bar的长度并指定坐标轴p值标的位置,#+guides(fill=FALSE)是移除图例
    dev.off()
}








setwd("/f/mulinlab/huan/ALL_result_ICGC_ALL_drug/gene_network_merge_repurposing_model/V2/validation/TCGA_28847918/")
org1<-read.table("./output/35_pandrug_top_0_0.35_final.txt",header = T,sep = "\t") %>% as.data.frame()
org2<-read.table("./output/35_pandrug_top_0.35_1_final.txt",header = T,sep = "\t") %>% as.data.frame()



rs<- bind_rows(org1,org2)
setwd("/f/mulinlab/huan/ALL_result_ICGC_ALL_drug/gene_network_merge_repurposing_model/V2/validation/TCGA_28847918/figure/")

pdf("36_pandrug_top_0.35_1_group_comparison_p.pdf",height = 4,width = 3.8) #把图片存下来
p <- ggplot(rs,aes(x=proportion,y=value_in_paper))+geom_boxplot(aes(fill=proportion),width=0.3,outlier.colour = NA)+ #按照标签画图
  scale_y_continuous(name = "Drug sensitivity in paper") +
  scale_x_discrete(name = "Top proportion")
my_comparisons <- list(c("0_0.35","0.35_1") ) 
p+coord_cartesian(ylim = c(-0.5,0.85))+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", method.args = list(alternative = "less"),
                                                           tip.length = 0.01,label.y=c(0.63,0.73,0.85))+guides(fill=FALSE)
##组内进行比较，指定p bar的长度并指定坐标轴p值标的位置,#+guides(fill=FALSE)是移除图例
dev.off()

#----------------violin
pdf("36_pandrug_top_0.35_1_group_comparison_violin.pdf",height = 4,width = 3.8) #把图片存下来
p1 <- ggplot(rs,aes(x=proportion,y=value_in_paper))+geom_violin()+geom_boxplot(aes(fill=proportion),width=0.3,outlier.colour = NA)+ #按照标签画图
  scale_y_continuous(name = "Drug sensitivity in paper") +
  scale_x_discrete(name = "Top proportion")
p1
dev.off()