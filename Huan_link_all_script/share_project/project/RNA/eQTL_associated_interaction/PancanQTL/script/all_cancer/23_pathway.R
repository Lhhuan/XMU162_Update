library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(patchwork)
library(Seurat)
# library(clusterProfiler)
# library(org.Hs.eg.db)
#-------------------------------------------------------
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/cancer_total/")
#-------------normal
org <- read.table("15_cancer_specfic_hotspot_gene.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org_sg <-filter(org,pvalue < 5e-8)
org_sg$cancer_tissue <-paste(org_sg$cancer,org_sg$tissue,sep=":")
org_sg_gene <-org_sg %>% dplyr::select(cancer_tissue,egene)%>%unique
unique_cancer_tissue <-unique(org_sg$cancer_tissue)

# gene <-data.frame(unique(filter(org_sg_gene,cancer_tissue==unique_cancer_tissue[6])%>%dplyr::select(egene)))
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/cancer_total/specific/")
# a <-(2:6,)
for(i in c(1:length(unique_cancer_tissue))){
    gene <-data.frame(unique(filter(org_sg_gene,cancer_tissue==unique_cancer_tissue[i])%>%dplyr::select(egene)))
    if(nrow(gene)>100){
        print(c(unique_cancer_tissue[i],nrow(gene)))
        cancer2 <- str_replace(unique_cancer_tissue[i],":",".")
        file = paste0(cancer2,".txt")
        write.table(gene,file,row.names = F, col.names = F,quote =F,sep="\t")
    }
   
    print(c(i,i))
}


all_gene <-org_sg %>% dplyr::select(egene)%>%unique()
write.table(all_gene,"all_gene.txt",row.names = F, col.names = F,quote =F,sep="\t")

#----------------plot cancer specific
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/cancer_total/specific/david_analyse/")

all <-function(class=NULL){
    a <- filter(h_sig,Category==class)   
    class <-gsub("_"," ",class)
    ggplot(a,aes(x=gene_ratio,y=reorder(Term,X=gene_ratio)))+
    geom_point(aes(size=Count*2,color=-1*log10(p.adjust)))+scale_color_gradient(low="blue",high ="red")+
    labs(color=expression(-log[10](p.adjust)),size="Gene",x="gene ratio",y=class)+theme_bw()+
    theme(axis.text.y=element_text(size=10.5),plot.title=element_text(hjust = 0.5))
}



for(tissue in unique_cancer_tissue){
# fig <-function(tissue){
    gene <-data.frame(unique(filter(org_sg_gene,cancer_tissue==tissue)%>%dplyr::select(egene)))
    if(nrow(gene)>100){
        # print(c(tissue,nrow(gene)))
        cancer2 <- str_replace(tissue,":",".")
        file = paste0(cancer2,".txt")
        david <- read.table(file,header = T,sep = "\t") %>% as.data.frame()
        colnames(david)[4]="gene_ratio"
        david$p.adjust <-p.adjust(david$PValue,method="fdr")
        h_sig <-filter(david,p.adjust<=0.05)%>%select("Category","Term","gene_ratio","Count","PValue","FDR","p.adjust")
        if(nrow(h_sig)>0){
            print(c(tissue,nrow(h_sig)))
            h_sig$Term <-gsub("^.*:","",h_sig$Term)
            h_sig$Term <-gsub("^.*~","",h_sig$Term)
            unique(h_sig$Term)
            # unique_Category <-unique(h_sig$Category)
            # unique_Category = factor(unique_Category, levels=c("KEGG_PATHWAY","GOTERM_CC_DIRECT", "GOTERM_BP_DIRECT","GOTERM_MF_DIRECT"))
            unique_Category <- c("KEGG_PATHWAY","GOTERM_CC_DIRECT", "GOTERM_BP_DIRECT","GOTERM_MF_DIRECT")
            # print(length(unique_Category))
            plist = lapply(unique_Category,all)
            fig_name <-paste0("./sig/",cancer2,".pdf")
            pdf(fig_name,width=15, height=10)
            p1 <-CombinePlots(plist,ncol=2,nrow=2)
            print(p1)
            dev.off()
            print(c(tissue,"fig"))
        }
    }
}

# lapply(unique_cancer_tissue,fig)


# enrichKK
# dotplot(enrichKK,title= "1111")
# dev.off()


# p1 <- barplot(enrichKK,showCategory = 20)
# print(p1)
# dev.off()


# barplot(enrichKK,showCategory = 20)
# dev.off()

# heatplot(enrichKK)
# dev.off()

# write.table(rs3,output_file3,row.names = F, col.names = T,quote =F,sep="\t")