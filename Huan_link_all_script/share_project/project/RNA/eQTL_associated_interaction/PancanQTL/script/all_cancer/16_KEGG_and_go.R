library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
#-------------------------------------------------------
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/cancer_total/")
#-------------normal
org <- read.table("15_cancer_specfic_hotspot_gene.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org_sg <-filter(org,pvalue < 5e-8)
org_sg$cancer_tissue <-paste(org_sg$cancer,org_sg$tissue,sep=":")
org_sg_gene <-org_sg %>% dplyr::select(cancer_tissue,egene)%>%unique
unique_cancer_tissue <-unique(org_sg$cancer_tissue)

gene <-data.frame(unique(filter(org_sg_gene,cancer_tissue==unique_cancer_tissue[6])%>%dplyr::select(egene)))

# gene <-data.frame(unique(filter(org_sg_gene,cancer_tissue=="KIRC:Kidney_Cortex")%>%dplyr::select(egene)))
gene <-data.frame(unique(filter(org_sg_gene,cancer_tissue=="UCS:Uterus")%>%dplyr::select(egene)))

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/cancer_total/KEGG/")

ssss <- c(2:length(unique_cancer_tissue))
ssss <-ssss[-c(7,23)]

a <-(2:6,)
for(i in ssss){
    gene <-data.frame(unique(filter(org_sg_gene,cancer_tissue==unique_cancer_tissue[i])%>%dplyr::select(egene)))

    gene_up_id <- bitr(gene$egene,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = "org.Hs.eg.db"
                    )

    enrichKK <- enrichKEGG(gene = gene_up_id$ENTREZID,
                        organism = "hsa",
                        pvalueCutoff = .05,
                        qvalueCutoff = .05)

    enrichKK = DOSE::setReadable(enrichKK,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')
    pathway_number = length(enrichKK$Description)
    if(pathway_number >0){
        cancer2 <- str_replace(unique_cancer_tissue[i],":",".")
        fig <- paste0(cancer2,".pdf")
        pdf(fig)
        p1 <-dotplot(enrichKK,title= unique_cancer_tissue[i])
        print(p1)
        dev.off()
        print(i)
    }
    print(c(i,i))
}


all_gene <-org_sg %>% dplyr::select(egene)%>%unique
write.table(all_gene,"all_gene.txt",row.names = F, col.names = T,quote =F,sep="\t")


enrichKK
dotplot(enrichKK,title= "1111")
dev.off()


# p1 <- barplot(enrichKK,showCategory = 20)
# print(p1)
# dev.off()


# barplot(enrichKK,showCategory = 20)
# dev.off()

# heatplot(enrichKK)
# dev.off()

# write.table(rs3,output_file3,row.names = F, col.names = T,quote =F,sep="\t")