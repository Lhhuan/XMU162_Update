library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/enrichment/interval_18/ALL/")
hotspot<-read.table("hotspot_15_state_count_Whole_Blood_cutoff_0.176.txt",header = T,sep = "\t") %>% as.data.frame()
random<-read.table("original_random_15_state_count_Whole_Blood_cutoff_0.176.txt",header = T,sep = "\t") %>% as.data.frame()

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

active_chromatin <-c("1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","8_ZNF/Rpts")
bivalent_chromatin<-c("10_TssBiv","11_BivFlnk","12_EnhBiv")
silent_chromatin <-c("9_Het","13_ReprPC","14_ReprPCWk","15_Quies")

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
#----1-8

# add_pppp=symnum.args[[2]][1]
add_pppp="****"
active_chromatin <-c("1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","8_ZNF/Rpts")
state = active_chromatin[8]
# for(state in active_chromatin){
    aa <-filter(rs,chromatin_state==state)
    Hotspot <-filter(aa,Class=="Hotspot")
    Random_expectaion <-filter(aa,Class!="Hotspot")
    observed <-Hotspot$number
    used_count_p <- nrow(filter(Random_expectaion,number >observed))
    p_value = used_count_p/nrow(Random_expectaion)


    title_name<-state
    state <-str_replace(state,"/","_")
    pdf(paste0("./figure/15_state_",state,".pdf"),width=4, height=4) 
    p <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ theme(legend.position ="none")+ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),axis.text.x = element_text(size = 12, color = "black"),axis.text.y = element_text(color = "black")) +ylab("Number of factor")+p_theme
    # add_pval(p,annotation = p_value)
    # add_pppp = "****"
    add_pval(p,annotation = add_pppp,pval_star = T)
    print(c(state,p_value))
    print(p)
    dev.off()
    # rm(p)
# }

bivalent_chromatin<-c("10_TssBiv","11_BivFlnk","12_EnhBiv")

state = bivalent_chromatin[3]
# for(state in bivalent_chromatin){
    aa <-filter(rs,chromatin_state==state)
    Hotspot <-filter(aa,Class=="Hotspot")
    Random_expectaion <-filter(aa,Class!="Hotspot")
    observed <-Hotspot$number
    used_count_p <- nrow(filter(Random_expectaion,number >observed))
    p_value = used_count_p/nrow(Random_expectaion)


    title_name<-state
    state <-str_replace(state,"/","_")
    pdf(paste0("./figure/15_state_",state,".pdf"),width=4, height=4) 
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

#----
silent_chromatin <-c("9_Het","13_ReprPC","14_ReprPCWk","15_Quies")

state = silent_chromatin[3]
# for(state in silent_chromatin){
    aa <-filter(rs,chromatin_state==state)
    Hotspot <-filter(aa,Class=="Hotspot")
    Random_expectaion <-filter(aa,Class!="Hotspot")
    observed <-Hotspot$number
    used_count_p <- nrow(filter(Random_expectaion,number >observed))
    p_value = used_count_p/nrow(Random_expectaion)


    title_name<-state
    state <-str_replace(state,"/","_")
    pdf(paste0("./figure/15_state_",state,".pdf"),width=4, height=4) 
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
state = silent_chromatin[4]
# for(state in silent_chromatin){
    aa <-filter(rs,chromatin_state==state)
    Hotspot <-filter(aa,Class=="Hotspot")
    Random_expectaion <-filter(aa,Class!="Hotspot")
    observed <-Hotspot$number
    used_count_p <- nrow(filter(Random_expectaion,number <observed))
    p_value = used_count_p/nrow(Random_expectaion)


    title_name<-state
    state <-str_replace(state,"/","_")
    pdf(paste0("./figure/15_state_",state,".pdf"),width=4, height=4) 
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

# # figure_list <-list()
# i=1

# for(state in unique_chromatin_state){
#     aa <-filter(rs,chromatin_state==state)
#     title_name<-state
#     state <-str_replace(state,"/","_")
#     figure_list[[i]] <-ggplot(aa,aes(x=Class,y=number))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ 
#     theme(legend.position ="none")+ggtitle(title_name)+theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank()) +ylab("Number of factor")+theme_bw()
#     i=i+1
#     print(state)
# }

# pdf("./figure/boxplot_compare_background_random_10000.pdf") 
# p1<-marrangeGrob(figure_list,nrow=2,ncol=2) 
# print(p1)
# dev.off()



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
