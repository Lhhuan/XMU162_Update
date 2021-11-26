library(readxl)
library(dplyr)
library(stringr)
library(igraph)

setwd("/home/huanhuan/project/kinase/kinase_ppi/script/")
di<-"~/project/kinase/kinase_ppi"

tfs<-read.table(file.path(di, "output/08_fisher_need_data.txt"),header = T,sep = "\t") %>% as.data.frame()
N<-nrow(tfs)
rs <- data.frame(symbol =character(),p_value = numeric()) #all result   #构造一个表
for(i in 1:N ){#对每一行进行循环
    # data<-c(tfs$kinase_in_the_top_final,tfs$kinase_out_of_top,tfs$gene_not_kinase_in_the_top,tfs$gene_not_kinase_out_of_the_top)
    symbol= tfs[i,1]
    # a= tfs[i,1]
    kinase_in_the_top_final= tfs[i,2] #提第i行的第二列
    gene_not_kinase_in_the_top = tfs[i,3]
    kinase_out_of_top = tfs[i,4]
    gene_not_kinase_out_of_the_top = tfs[i,5]
    data<-c(kinase_in_the_top_final,kinase_out_of_top,gene_not_kinase_in_the_top,gene_not_kinase_out_of_the_top)
    data2<-matrix(data,nrow=2)
    test_result <-fisher.test(data2)  #%>%as.data.frame()
    p1 <-test_result$p.value
    temp <- data.frame(symbol = symbol,p_value=p1)
    rs <- bind_rows(rs,temp) ##把数据集temp作为新的行添加到rs中
}

write.table(rs,file.path(di,"output/09_fisher_test_result.txt"),quote = F,sep = "\t",row.names = F)