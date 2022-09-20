library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(Seurat)
library(reshape2)
library(R.utils)
# library(mclust)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/Knowledge_Graph/script/prediction/script/chr1/output/")
gidx <-read.table("01_gene_idx.txt.gz",header = T,sep = "\t")%>% as.data.frame()
hidx <-read.table("01_hotspot_idx.txt.gz",header = T,sep = "\t")%>% as.data.frame()
hotspot_egene <-read.table("02_hotspot_egene_idx.txt.gz",header = T,sep = "\t")%>% as.data.frame()
interaction <-read.table("02_egene_interaction_idx.txt.gz",header = T,sep = "\t")%>% as.data.frame()
co_expression <-read.table("02_egene_co-expression_idx.txt.gz",header = T,sep = "\t")%>% as.data.frame()

hotspot_egene_new <-hotspot_egene[sample(1:nrow(hotspot_egene), nrow(hotspot_egene), replace = FALSE),]
line=round(nrow(hotspot_egene_new)*0.1)
train_h <- hotspot_egene_new[1:(line*8),]
test_p <- hotspot_egene_new[(line*8+1):(line*9),]
val_p <- hotspot_egene_new[(line*9+1):nrow(hotspot_egene_new),]

set.seed(1)
random_h <- hidx[sample(1:nrow(hidx), nrow(train_h), replace = TRUE),"hotspot_idx"]
random_e <- gidx[sample(1:nrow(gidx), nrow(train_h), replace = TRUE),"Gene_idx"]

rs <-data.frame(Hotspot_idx=random_h,Egene_idx=random_e)
potential_negative <- anti_join(rs,hotspot_egene,by=c("Hotspot_idx","Egene_idx"))

test_n <- potential_negative[1:nrow(test_p),]
val_n  <- potential_negative[(nrow(test_p)+1):(nrow(test_p)+nrow(val_p)),]

test_n$label<-0
test_p$label<-1
test <-bind_rows(test_p,test_n)
test$edge_type <-1
test <-test[,c(4,1:3)]

val_n$label<-0
val_p$label<-1
val <-bind_rows(val_p,val_n)
val$edge_type <-1
val <-val[,c(4,1:3)]

write.table(test,"./huan_GATNE/test.txt",row.names = F, col.names = F,quote =F,sep=" ")
write.table(val,"./huan_GATNE/valid.txt",row.names = F, col.names = F,quote =F,sep=" ")
#---------

train_h$edge_type <-1


# for(i in nrow(interaction)){
inter <-function(i=NULL){
    vs=unlist(str_split(interaction[i,"ReactomeFI_idx"], ";"))
    rs <-data.frame()
    for(j in 1:length(vs)){
        tmp =data.frame(Gene1=interaction[i,"Egene_idx"],Gene2=as.numeric(vs[j]))
        rs <-bind_rows(rs,tmp)
    }
    return(rs)
}
interaction_re <-lapply(1:nrow(interaction),inter)
f_inter <-do.call(rbind,interaction_re)
f_inter$edge_type <-2

co <-function(i=NULL){
    vs=unlist(str_split(co_expression[i,"co_expression_idx"], ";"))
    rs <-data.frame()
    for(j in 1:length(vs)){
        tmp =data.frame(Gene1=co_expression[i,"Egene_idx"],Gene2=as.numeric(vs[j]))
        rs <-bind_rows(rs,tmp)
    }
    return(rs)
}
co_expression_re <-lapply(1:nrow(co_expression),co)
f_coexpre <-do.call(rbind,co_expression_re)
f_coexpre$edge_type <-3

colnames(train_h)[1:2] <-c("Gene1","Gene2")
train <-bind_rows(train_h,f_inter,f_coexpre)
train <-train[,c(3,1,2)]
write.table(train,"./huan_GATNE/train.txt",row.names = F, col.names = F,quote =F,sep=" ")

train_hf <-train_h[,c(3,1,2)]
write.table(train_hf,"./huan_GATNE/train.txt",row.names = F, col.names = F,quote =F,sep=" ")
