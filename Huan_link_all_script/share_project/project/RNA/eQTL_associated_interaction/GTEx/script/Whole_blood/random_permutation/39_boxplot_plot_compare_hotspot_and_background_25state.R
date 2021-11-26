library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/enrichment/interval_18/ALL/")
hotspot<-read.csv("hotspot_25_state_count_Whole_Blood_cutoff_0.176.csv",header = T) %>% as.data.frame()
random<-read.csv("original_random_25_state_count_Whole_Blood_cutoff_0.176.csv",header = T) %>% as.data.frame()

hotspot$random_number <- 1
hotspot$Class <- "Hotspot"
random$Class <- "Random expectaion"
rs<- bind_rows(hotspot,random)

unique_chromatin_state <- unique(rs$chromatin_state)



p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 13),
                                                # axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))


#---------------
#-------------"1_TssA"
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

active_chromatin <-c("1_TssA","2_PromU","3_PromD1","4_PromD2","5_Tx5'","6_Tx","7_Tx3'","8_TxWk","9_TxReg","10_TxEnh5'","11_TxEnh3'","12_TxEnhW","13_EnhA1","14_EnhA2","15_EnhAF","16_EnhW1","17_EnhW2","18_EnhAc","19_DNase","20_ZNF_Rpts")
bivalent_chromatin<-c("22_PromP","23_PromBiv")
silent_chromatin <-c("21_Het","24_ReprPC","25_Quies")

# state = "1_TssA"
# state = "15_Quies"
# for(state in unique_chromatin_state){
#     aa <-filter(rs,chromatin_state==state)
#     Hotspot <-filter(aa,Class=="Hotspot")
#     Random_expectaion <-filter(aa,Class!="Hotspot")
#     observed <-Hotspot$number
#     used_count_p <- nrow(filter(Random_expectaion,number >observed))
#     p_value = used_count_p/nrow(Random_expectaion)


#     title_name<-state
#     state <-str_replace(state,"/","_")
#     pdf(paste0("./figure/15_state_",state,".pdf"),width=4, height=4) 
#     p <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ theme(legend.position ="none")+ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),axis.text.x = element_text(size = 12, color = "black"),axis.text.y = element_text(color = "black")) +ylab("Number of factor")+p_theme
#     # add_pval(p,annotation = p_value)
#     add_pval(p,annotation = p_value)
#     print(c(state,p_value))
#     print(p)
#     dev.off()
#     rm(p)
# }
#----active_chromatin

add_pppp="****"
active_chromatin <-c("1_TssA","2_PromU","3_PromD1","4_PromD2","5_Tx5'","6_Tx","7_Tx3'","8_TxWk","9_TxReg","10_TxEnh5'","11_TxEnh3'","12_TxEnhW","13_EnhA1","14_EnhA2","15_EnhAF","16_EnhW1","17_EnhW2","18_EnhAc","19_DNase","20_ZNF/Rpts")
state = active_chromatin[19]
# for(state in active_chromatin){
    aa <-filter(rs,chromatin_state==state)
    Hotspot <-filter(aa,Class=="Hotspot")
    Random_expectaion <-filter(aa,Class!="Hotspot")
    observed <-Hotspot$number
    used_count_p <- nrow(filter(Random_expectaion,number >observed))
    p_value = used_count_p/nrow(Random_expectaion)


    title_name<-state
    state <-str_replace(state,"/","_")
    state <-str_replace(state,"'","")
    pdf(paste0("./figure/25_state_",state,".pdf"),width=4, height=4) 
    p <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ theme(legend.position ="none")+ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),axis.text.x = element_text(size = 12, color = "black"),axis.text.y = element_text(color = "black")) +ylab("Number of factor")+p_theme
    # add_pval(p,annotation = p_value)
    add_pppp = "*"
    add_pval(p,annotation = add_pppp,pval_star = T)
    print(c(state,p_value))
    print(p)
    dev.off()
    # # rm(p)
# }

bivalent_chromatin<-c("22_PromP","23_PromBiv")

state = bivalent_chromatin[2]
# for(state in bivalent_chromatin){
    aa <-filter(rs,chromatin_state==state)
    Hotspot <-filter(aa,Class=="Hotspot")
    Random_expectaion <-filter(aa,Class!="Hotspot")
    observed <-Hotspot$number
    used_count_p <- nrow(filter(Random_expectaion,number >observed))
    p_value = used_count_p/nrow(Random_expectaion)


    title_name<-state
    state <-str_replace(state,"/","_")
    state <-str_replace(state,"'","")
    pdf(paste0("./figure/25_state_",state,".pdf"),width=4, height=4) 
    p <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ theme(legend.position ="none")+ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),axis.text.x = element_text(size = 12, color = "black"),axis.text.y = element_text(color = "black")) +ylab("Number of factor")+p_theme
    # add_pval(p,annotation = p_value)
    # add_pppp = "****"
    # add_pval(p,annotation = "****",pval_star = T)
    add_pval(p,annotation = "*",pval_star = T)
    print(c(state,p_value))
    print(p)
    dev.off()
    # rm(p)
# }
silent_chromatin <-c("21_Het","24_ReprPC","25_Quies")

state = silent_chromatin[2]
# for(state in silent_chromatin){
    aa <-filter(rs,chromatin_state==state)
    Hotspot <-filter(aa,Class=="Hotspot")
    Random_expectaion <-filter(aa,Class!="Hotspot")
    observed <-Hotspot$number
    used_count_p <- nrow(filter(Random_expectaion,number >observed))
    p_value = used_count_p/nrow(Random_expectaion)


    title_name<-state
    state <-str_replace(state,"/","_")
    state <-str_replace(state,"'","")
    pdf(paste0("./figure/25_state_",state,".pdf"),width=4, height=4) 
    p <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ theme(legend.position ="none")+ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),axis.text.x = element_text(size = 12, color = "black"),axis.text.y = element_text(color = "black")) +ylab("Number of factor")+p_theme
    # add_pval(p,annotation = p_value)
    # add_pppp = "****"
    add_pval(p,annotation = "****",pval_star = T)
    # add_pval(p,annotation = "ns",pval_star = T)
    print(c(state,p_value))
    print(p)
    dev.off()
# }

#----------"15_Quies"
state = silent_chromatin[3]
# for(state in silent_chromatin){
    aa <-filter(rs,chromatin_state==state)
    Hotspot <-filter(aa,Class=="Hotspot")
    Random_expectaion <-filter(aa,Class!="Hotspot")
    observed <-Hotspot$number
    used_count_p <- nrow(filter(Random_expectaion,number <observed))
    p_value = used_count_p/nrow(Random_expectaion)


    title_name<-state
    state <-str_replace(state,"/","_")
    pdf(paste0("./figure/25_state_",state,".pdf"),width=4, height=4) 
    p <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ theme(legend.position ="none")+ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),axis.text.x = element_text(size = 12, color = "black"),axis.text.y = element_text(color = "black")) +ylab("Number of factor")+p_theme
    # add_pval(p,annotation = p_value)
    # add_pppp = "****"
    add_pval(p,annotation = "****",pval_star = T)
    # add_pval(p,annotation = "ns",pval_star = T)
    print(c(state,p_value))
    print(p)
    dev.off()
#----




# #-------------
# for(state in unique_chromatin_state){
#     aa <-filter(rs,chromatin_state==state)
#     title_name<-state
#     state <-str_replace(state,"/","_")
#     # figure_list[[i]] <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ theme
#     pdf(paste0("./figure/15_state_",state,".pdf")) 
#     p <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ theme(legend.position ="none")+ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank()) +ylab("Number of factor")
#     print(state)
#     print(p)
#     dev.off()
# }
# #---------------





# #------------
# figure_list <-list()
# # my_comparisons <- list(c("hotspot","random"))

figure_list <-list()
i=1

for(state in unique_chromatin_state){
    aa <-filter(rs,chromatin_state==state)
    title_name<-state
    state <-str_replace(state,"/","_")
    figure_list[[i]] <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ 
    theme(legend.position ="none")+ggtitle(title_name)+theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank()) +ylab("Number of factor")
    i=i+1
    print(state)
}

pdf("./figure/boxplot_compare_background_state25.pdf") 
p1<-marrangeGrob(figure_list,nrow=2,ncol=2) 
print(p1)
dev.off()



# #----------------------
# figure_list <-list()
# i=1

# for(state in unique_chromatin_state){
#     #---
#     aa <-filter(rs,chromatin_state==state)
#     Hotspot <-filter(aa,Class=="Hotspot")
#     Random_expectaion <-filter(aa,Class!="Hotspot")
#     observed <-Hotspot$number
#     used_count_p <- nrow(filter(Random_expectaion,number >observed))
#     p_value = used_count_p/nrow(Random_expectaion)
#     # #---
#     # aa <-filter(rs,chromatin_state==state)
#     # title_name<-state
#     # state <-str_replace(state,"/","_")
#     # figure_list[[i]] <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ 
#     # theme(legend.position ="none")+ggtitle(title_name)+theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank()) +ylab("Number of factor")+theme_bw()
#     # i=i+1
#     print(c(state,p_value))
#     # print(p_value)
# }

# pdf("./figure/boxplot_compare_background_random_10000_p.pdf") 
# p1<-marrangeGrob(figure_list,nrow=2,ncol=2) 
# print(p1)
# dev.off()
