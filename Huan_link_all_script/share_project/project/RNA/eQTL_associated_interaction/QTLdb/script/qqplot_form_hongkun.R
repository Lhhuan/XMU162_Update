all_mRNA <- read.table(paste0("/state/partition1/hongkun/Genotype_TCGA/Pancancer_exon_1M_pvOutputThreshold_1/",cancer,"/all_mRNA_",cancer,"_tra_greater_than1M.eqtl.gz"),sep="\t",header=T)
 png(paste0("/state/partition1/hongkun/Genotype_TCGA/Pancancer_exon_1M_pvOutputThreshold_1/",cancer,"/",cancer,"_qqplot_mRNA_tra_all.png"))
 par(las = 1, pty="s")
 plot(c(0, 7.5), c(0, 7.5), col = 8, lwd = 3, type= "l", xlab="Expected (-log10P)", ylab="Observed (-log10P)", xlim = c(0, 7.5), ylim = c(0, 7.5), xaxs="i", yaxs="i", bty="l", main = paste(" raw Pvalue qqplot of mRNA in PRAD"))
 pvalue <- all_mRNA$p.value
 o <- -log10(sort(na.omit(pvalue), decreasing=F))
 e <- -log10(1:length(o) / (length(o) + 1))
 points(e, o, pch= 18, cex= .4, col = "blue")
 lab <- quantile(o, 0.475) / quantile(e, 0.475)
 legend("bottomright", legend = paste(expression("lambda = "), signif(lab, 3), sep = ""), bty = "n")
 dev.off()


x <-org2$emplambda
x[which(x==0)]=x[which(x==0)]+1e-12


png("1234.png")

par(las = 1, pty="s")
plot(c(0, 7.5), c(0, 7.5), col = 8, lwd = 3, type= "l", xlab="Expected (-log10P)", ylab="Observed (-log10P)", xlim = c(0, 7.5), ylim = c(0, 7.5), xaxs="i", yaxs="i", bty="l", main = paste(" raw Pvalue qqplot of mRNA in PRAD"))
 pvalue <- x
 o <- -log10(sort(na.omit(pvalue), decreasing=F))
 e <- -log10(1:length(o) / (length(o) + 1))
 points(e, o, pch= 18, cex= .4, col = "blue")
 lab <- quantile(o, 0.475) / quantile(e, 0.475)
 legend("bottomright", legend = paste(expression("lambda = "), signif(lab, 3), sep = ""), bty = "n")
 dev.off()