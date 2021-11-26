library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/enrichment/interval_18/ALL/")
hotspot<-read.table("hotspot_cutoff_0.176_histone_marker_jaccard_index.txt.gz",header = T,sep = "\t") %>% as.data.frame()
# random<-fread("0_0.176_jaccard_index_histone_marker_1000.txt.gz",header = T,sep = "\t")
random<-read.table("0_0.176_jaccard_index_histone_marker_1000.txt.gz",header = T,sep = "\t") %>% as.data.frame()


active_mark <-c("H3K27ac","H3K9ac","H3K36me3","H3K4me1","H3K4me3")
repressed_mark <-c("H3K27me3","H3K9me3")

hotspot$Random_number <- 1
hotspot$Class <- "hotspot"
random$Class <- "random"

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 8),
                                                # axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))

rs<- bind_rows(hotspot,random)


# unique_chromatin_state <- unique(rs$chromatin_state)
figure_list <-list()
my_comparisons <- list(c("hotspot","random"))

Marker ="H3K27ac"
aa <-filter(rs,Marker==Marker)

title_name <-Marker
pdf("./figure/0_0.176/0_0.176_hotspot_jaccard_H3K27ac_no_limit.pdf")
p <-ggplot(aa,aes(x=Class,y=jaacard_index))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ theme(legend.position ="none")+ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "greater"))
print(p)
dev.off()


pdf("./figure/0_0.176/0_0.176_hotspot_jaccard_H3K27ac_no_limit_org.pdf")
p <-ggplot(aa,aes(x=Class,y=jaacard_index))+geom_boxplot(aes(fill=Class),width=0.3)+ theme(legend.position ="none")+ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "greater"))
print(p)
dev.off()


png("./figure/0_0.176/0_0.176_hotspot_jaccard_H3K27ac_no_limit_org.png")
p <-ggplot(aa,aes(x=Class,y=jaacard_index))+geom_boxplot(aes(fill=Class),width=0.3)+ theme(legend.position ="none") +ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "greater"))
print(p)
dev.off()



png("./figure/0_0.176/0_0.176_hotspot_jaccard_H3K27ac_violin_boxplot_org.png")
p <-ggplot(aa,aes(x=Class,y=jaacard_index,fill = Class))+geom_boxplot() + geom_violin()+ theme(legend.position ="none") +ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "greater"))
print(p)
dev.off()

png("./figure/0_0.176/0_0.176_hotspot_jaccard_H3K27ac_violin_boxplot_without_outlier.png")
# bbb <-head(aa,400000)
p <-ggplot(aa,aes(x=Class,y=jaacard_index,fill = Class))+geom_boxplot(outlier.colour = NA) + geom_violin()+ theme(legend.position ="none") +ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "greater"))+
    p_theme
print(p)
dev.off()


pdf("./figure/0_0.176_hotspot_jaccard_0_0.25.pdf")
p <-ggplot(aa,aes(x=Class,y=jaacard_index))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ theme(legend.position ="none")+ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "greater"))+
    scale_y_continuous(limits=c(0,0.25))

print(p)
dev.off()


for(Marker in active_mark ){
    aa <-filter(rs,Marker==Marker)
    title_name<-state
    state <-str_replace(state,"/","_")
    # figure_list[[i]] <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ theme
    pdf(paste0("15_state_",state,".pdf")) 
    p <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ theme(legend.position ="none")+ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5))+
    # p <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ 
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "greater"))
    i=i+1
    print(state)
    print(p)
    dev.off()
}

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