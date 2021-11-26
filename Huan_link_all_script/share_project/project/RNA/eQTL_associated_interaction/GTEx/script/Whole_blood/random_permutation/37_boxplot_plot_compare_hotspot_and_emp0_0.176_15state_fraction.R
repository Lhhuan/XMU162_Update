library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(Seurat)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/enrichment/interval_18/ALL/")
hotspot<-read.table("hotspot_15_state_count_Whole_Blood_cutoff_0.176.txt",header = T,sep = "\t") %>% as.data.frame()
random<-read.table("0_0.176_15_state_count_Whole_Blood_cutoff_0.176.txt",header = T,sep = "\t") %>% as.data.frame()

hotspot$random_number <- 1
hotspot$Class <- "Hotspot"
random$Class <- "Random expectaion"
hotspot$Fraction <- hotspot$number/46369
random$Fraction <-random$number/46369
rs<- bind_rows(hotspot,random)

unique_chromatin_state <- unique(rs$chromatin_state)



p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 8),
                                                # axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))


#---------------
#-------------"1_TssA"
symnum.args <- list(cutpoints = c(0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

active_chromatin <-c("1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","8_ZNF/Rpts")
bivalent_chromatin<-c("10_TssBiv","11_BivFlnk","12_EnhBiv")
silent_chromatin <-c("9_Het","13_ReprPC","14_ReprPCWk","15_Quies")
silent_chromatin1 <-c("9_Het","13_ReprPC","14_ReprPCWk")

# state = "1_TssA"
# state = "15_Quies"
# for(state in unique_chromatin_state){
#     aa <-filter(rs,chromatin_state==state)
#     Hotspot <-filter(aa,Class=="Hotspot")
#     Random_expectaion <-filter(aa,Class!="Hotspot")
#     observed <-Hotspot$number
#     used_count_p <- nrow(filter(Random_expectaion,number >observed))
#     p_value = used_count_p/nrow(Random_expectaion)
#     print(p_value)

#     title_name<-state
#     state <-str_replace(state,"/","_")
#     #pdf(paste0("./figure/Fraction/15_state_",state,".pdf"),width=4, height=4) 
#     p <-ggplot(aa,aes(x=Class,y=Fraction))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ 
#     #scale_y_continuous(limits=c(0,5)) + 
#     theme(legend.position ="none")+ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),axis.text.x = element_text(size = 10, color = "black"),axis.text.y = element_text(color = "black")) +ylab("Fraction of segments")+p_theme
#     # add_pval(p,annotation = p_value)
#     pdf(paste0("./figure/random_permutation/15_state_",state,".pdf"),width=4, height=4)
#     print(add_pval(p,annotation = p_value))
#     print(c(state,p_value))
#     #print(p)
#     dev.off()
#     # rm(p)
# }


large <-function(state){
    aa <-filter(rs,chromatin_state==state)
    Hotspot <-filter(aa,Class=="Hotspot")
    Random_expectaion <-filter(aa,Class!="Hotspot")
    observed <-Hotspot$number
    used_count_p <- nrow(filter(Random_expectaion,number >observed))
    p_value = used_count_p/nrow(Random_expectaion)
    print(p_value)

    title_name<-state
    state <-str_replace(state,"/","_")
    #pdf(paste0("./figure/Fraction/15_state_",state,".pdf"),width=4, height=4) 
    p <-ggplot(aa,aes(x=Class,y=Fraction))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ 
    #scale_y_continuous(limits=c(0,5)) + 
    theme(legend.position ="none")+ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),axis.text.x = element_text(size = 7, color = "black"),axis.text.y = element_text(color = "black")) +ylab("Fraction of segments")+p_theme

    #add_pval(p,annotation = p_value)

    if(p_value <= 0.0001){
        add_pval(p,annotation = "****",pval_star = T)
    }else if(p_value <= 0.001){
        add_pval(p,annotation = "***",pval_star = T)
    }else if(p_value <= 0.01){
        add_pval(p,annotation = "**",pval_star = T)
    }else if(p_value <= 0.05){
        add_pval(p,annotation = "*",pval_star = T)
    }else{
        add_pval(p,annotation = "ns",pval_star = T)
    }
    
}
 
plist = lapply(active_chromatin,large)

pdf("./figure/percentage/emplambda0_0.176_15_state_active_chromatin_Fraction.pdf",width=8.8, height=4)
CombinePlots(plist,ncol=4,nrow=2)
dev.off()

plist = lapply(bivalent_chromatin,large)

pdf("./figure/percentage/emplambda0_0.176_15_state_bivalent_chromatin_Fraction.pdf",width=6.6, height=2)
CombinePlots(plist,ncol=3)
dev.off()

plist = lapply(silent_chromatin1,large)

pdf("./figure/percentage/emplambda0_0.176_15_state_silent_chromatin1_Fraction.pdf",width=6.6, height=2)
CombinePlots(plist,ncol=3)
dev.off()

small <-function(state){
    aa <-filter(rs,chromatin_state==state)
    Hotspot <-filter(aa,Class=="Hotspot")
    Random_expectaion <-filter(aa,Class!="Hotspot")
    observed <-Hotspot$number
    used_count_p <- nrow(filter(Random_expectaion,number <observed))
    p_value = used_count_p/nrow(Random_expectaion)
    print(p_value)

    title_name<-state
    state <-str_replace(state,"/","_")
    #pdf(paste0("./figure/Fraction/15_state_",state,".pdf"),width=4, height=4) 
    p <-ggplot(aa,aes(x=Class,y=Fraction))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ 
    #scale_y_continuous(limits=c(0,5)) + 
    theme(legend.position ="none")+ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),axis.text.x = element_text(size = 7, color = "black"),axis.text.y = element_text(color = "black")) +ylab("Fraction of segments")+p_theme

    #add_pval(p,annotation = p_value)

    if(p_value <= 0.0001){
        add_pval(p,annotation = "****",pval_star = T)
    }else if(p_value <= 0.001){
        add_pval(p,annotation = "***",pval_star = T)
    }else if(p_value <= 0.01){
        add_pval(p,annotation = "**",pval_star = T)
    }else if(p_value <= 0.05){
        add_pval(p,annotation = "*",pval_star = T)
    }else{
        add_pval(p,annotation = "ns",pval_star = T)
    }
    
}
plist = lapply("15_Quies",small)
pdf("./figure/percentage/emplambda0_0.176_15_state_Quies_Fraction.pdf",width=2.2, height=2)
# CombinePlots(plist,ncol=3)
plist
dev.off()

large_fc <-function(state){
    aa <-filter(rs,chromatin_state==state)
    Hotspot <-filter(aa,Class=="Hotspot")
    Random_expectaion <-filter(aa,Class!="Hotspot")
    observed <-Hotspot$number
    used_count_p <- nrow(filter(Random_expectaion,number >observed))
    p_value = used_count_p/nrow(Random_expectaion)
    print(p_value)

    title_name<-state
    state <-str_replace(state,"/","_")
    fold_change=Hotspot$Fraction/(mean(Random_expectaion$Fraction))
    tmp <-data.frame(Marker=state,fold_change=fold_change)
    print(c(state,fold_change))
    # print("aaa")
    return(tmp)
    
}
small_fc <-function(state){
    aa <-filter(rs,chromatin_state==state)
    Hotspot <-filter(aa,Class=="Hotspot")
    Random_expectaion <-filter(aa,Class!="Hotspot")
    observed <-Hotspot$number
    used_count_p <- nrow(filter(Random_expectaion,number >observed))
    p_value = used_count_p/nrow(Random_expectaion)
    print(p_value)

    title_name<-state
    state <-str_replace(state,"/","_")
    fold_change=(mean(Random_expectaion$Fraction))/Hotspot$Fraction
    tmp <-data.frame(Marker=state,fold_change=fold_change)
    print(c(state,fold_change))
    # print("aaa")
    return(tmp)
    
}

large_state=c(active_chromatin,bivalent_chromatin,silent_chromatin1)
plist = lapply(large_state,large_fc)
large_result <-do.call(rbind,plist)

plist = lapply("15_Quies",small_fc)
small_result <-do.call(rbind,plist)
final_result <-bind_rows(large_result,small_result)
write.table(final_result,"./figure/percentage/compare_emplambda0_0.176_15_state_fraction_fold_change.txt",row.names = F, col.names = T,quote =F,sep="\t")