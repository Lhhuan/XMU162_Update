ARGS <- read.delim("filter.args", stringsAsFactor = F, sep = "=", quote = "\"", head = F, row.names = 1)
CIRCOS <- ARGS["circos", 1]
#plot circos using OmicCircos
options(stringsAsFactors = FALSE)
set.seed(999)
library("OmicCircos")
data("UCSC.hg19.chr")
#data("TCGA.BC.gene.exp.2k.60")
dat <- UCSC.hg19.chr
#dat$chrom <- gsub("chr", "",dat$chrom)
indx = which(dat$chrom != "chrX" & dat$chrom != "chrY")
seg.d = dat[indx, ]
seg.name = paste("chr", seq(1, 22, by = 1), sep = "")
seg.c = segAnglePo(seg.d, seg = seg.name)
colors <- rainbow(22, alpha = 0.8)

####################AI_focal circos
ai.clust = read.table(paste(CIRCOS, "AI_focal_clust_inputMC1.4.txt", sep = "/"), header = F, sep = " ", stringsAsFactors = F)
ai.clust = read.table("AI_focal_clust_inputMC1.4.txt", sep = "/"), header = F, sep = " ", stringsAsFactors = F)
indx = which(ai.clust[, 6] < 10)
ai.clustSL = ai.clust[indx, 1:6]
ai.clustSM = ai.clust[-indx, 1:6]
ai.clust = ai.clust[, 1:4]
# out.file <- paste(CIRCOS, "AI_focal_clust_circosMC1.4.pdf", sep = "/")
# pdf(file = out.file, height = 8, width = 8, compress = TRUE)
pdf(file = "test_circos.pdf", height = 8, width = 8, compress = TRUE)
par(las = 1, mar = c(0.5, 0.5, 1, 0.5))
plot(c(1,800), c(1,800), type= "n", axes = FALSE, xlab = "", ylab = "")
circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
circos(R = 10, cir = seg.c, type= "arc", W = 380, mapping = ai.clust, B = TRUE, col = colors[c(1)], lwd = 5, col.v = 4)

plot(c(1, 800), c(1, 800), type= "n", axes = FALSE, xlab = "", ylab = "", main = "Cutoff Seg_length = 10MB", cex.main = 0.8)
circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
circos(R = 300, cir = seg.c, type = "arc", W = 90, mapping = ai.clustSM, B = TRUE, col = "purple", lwd = 5, col.v = 6)
circos(R = 20, cir = seg.c, type = "arc", W = 275, mapping = ai.clustSL, B = T, col = "orange", lwd = 5, col.v = 6)
dev.off()


#---------------------------------------------
pdf(file = "test_circos1.pdf", height = 8, width = 8, compress = TRUE)
par(las = 1, mar = c(0.5, 0.5, 1, 0.5))
plot(c(1,800), c(1,800), type= "n", axes = FALSE, xlab = "", ylab = "")
circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
circos(R = 10, cir = seg.c, type= "heatmap2", W = 380, mapping = ai.clust, B = TRUE, col = colors[c(1)], lwd = 5, col.v = 4)

plot(c(1, 800), c(1, 800), type= "n", axes = FALSE, xlab = "", ylab = "", main = "Cutoff Seg_length = 10MB", cex.main = 0.8)
circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
circos(R = 300, cir = seg.c, type = "arc", W = 90, mapping = ai.clustSM, B = TRUE, col = "purple", lwd = 5, col.v = 6)
circos(R = 20, cir = seg.c, type = "arc", W = 275, mapping = ai.clustSL, B = T, col = "orange", lwd = 5, col.v = 6)
dev.off()