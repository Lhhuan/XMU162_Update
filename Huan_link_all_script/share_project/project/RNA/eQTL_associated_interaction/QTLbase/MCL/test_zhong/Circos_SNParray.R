ARGS <- read.delim("filter.args", stringsAsFactor = F, sep = "=", quote = "\"", head = F, row.names = 1)
PROJECT <- ARGS["project", 1]
WORKDIR <- ARGS["workdir", 1]
PROCDIR <- ARGS["procdir", 1]
FILTDIR  <- ARGS["filtdir", 1]
CDF <- ARGS["cdf", 1]
GENOME <- ARGS["genome", 1]
ANNOTATION <- ARGS["probe.info", 1] 
ANCESTRY <- ARGS["population.info", 1]
CLINICAL <- ARGS["clin.info", 1]
CYTOBAND <- ARGS["cytoband", 1]
PURITY <- ARGS["purity", 1]
CIRCOS <- ARGS["circos", 1]

setwd(WORKDIR)

library(rARPACK)
library(parallel)
#library(fastcluster)

#load(ANNOTATION)
load(CYTOBAND)

source("genomefilter.5.2.R")

files = list.files(path = FILTDIR, pattern = "*.MCclust1.4", full = T)
ai.focal = grep("AI_focal", files, value = T)
del.focal <- grep("loss_focal", files, value = T)
amp.focal <- grep("gain_focal", files, value = T)

sapply(del.focal, function(f) {
	cat(f, "\n")
	seg = try(read.table(f, stringsAsFactor = F, head = F))
	if (class(seg) == "try-error") return(NA)
	#seg[, 1] = sub("chr", "hs", seg[, 1])
	cl <- sort(unique(seg[, 6]))
	posit = t(sapply(cl, function(i) {
		indx = which(seg[, 6] == i & !is.na(seg[, 1]))
		seg_start <- min(seg[indx, 2])
		seg_end <- max(seg[indx, 3])
		seg_length = (seg_end - seg_start) / 1E6
		write.table(matrix(c(seg[1, 1], seg_start, seg_end, i, seg[indx, 10][1], seg_length), nr = 1), quote = F, col.names = F, row.names = F, file = paste(CIRCOS, "/loss_focal_clust_inputMC1.4.txt", sep = ""), append = T)
}))
	posit = t(sapply(cl, function(i) {
		indx = which(seg[, 6] == i & !is.na(seg[, 1]) & seg[, 10] < 2.5)
		if(length(indx) == 0) return(NULL)
		seg_start <- min(seg[indx, 2])
		seg_end <- max(seg[indx, 3])
		seg_length = (seg_end - seg_start) / 1E6
		write.table(matrix(c(seg[1, 1], seg_start, seg_end, i, seg[indx, 10][1], seg_length), nr = 1), quote = F, col.names = F, row.names = F, file = paste(CIRCOS, "/loss_focal_clustL2.5_inputMC1.4.txt", sep = ""), append = T)
	}))
	posit = t(sapply(cl, function(i) {
		indx = which(seg[, 6] == i & !is.na(seg[, 1]) & seg[, 10] >= 2.5)
		if(length(indx) == 0) return(NULL)
		seg_start <- min(seg[indx, 2])
		seg_end <- max(seg[indx, 3])
		seg_length = (seg_end - seg_start) / 1E6
		write.table(matrix(c(seg[1, 1], seg_start, seg_end, i, seg[indx, 10][1], seg_length), nr = 1), quote = F, col.names = F, row.names = F, file = paste(CIRCOS, "/loss_focal_clustU2.5_inputMC1.4.txt", sep = ""), append = T)
	}))

	ans = seg[, c(1:3, 6, 10)]
	write.table(ans, quote = F, row.names = F, col.names = F, file = paste(CIRCOS, "/loss_focal_inputMC1.4.txt", sep = ""), append = T)
	indx = which(seg[, 10] < 2.5)
	if(length(indx) == 0) {
		write.table(ans, quote = F, row.names = F, col.names = F, file = paste(CIRCOS, "/loss_focal_more2.5_inputMC1.4.txt", sep = ""), append = T)
	} else {
		write.table(ans[indx, ], quote = F, row.names = F, col.names = F, file = paste(CIRCOS, "/loss_focal_less2.5_inputMC1.4.txt", sep = ""), append = T)
		write.table(ans[-indx, ], quote = F, row.names = F, col.names = F, file = paste(CIRCOS, "/loss_focal_more2.5_inputMC1.4.txt", sep = ""), append = T)
	}
})

sapply(amp.focal, function(f) {
	cat(f, "\n")
	seg = try(read.table(f, stringsAsFactor = F, head = F))
	if (class(seg) == "try-error") return(NA)
	#seg[, 1] = sub("chr", "hs", seg[, 1])
	cl <- sort(unique(seg[, 6]))
	posit = t(sapply(cl, function(i) {
		indx = which(seg[, 6] == i & !is.na(seg[, 1]))
		seg_start <- min(seg[indx, 2])
		seg_end <- max(seg[indx, 3])
			#seg_start = round(quantile(seg[indx, 2], 0.25))
			#seg_end = round(quantile(seg[indx, 3], 0.75))
		seg_length = (seg_end - seg_start) / 1E6
		write.table(matrix(c(seg[1, 1], seg_start, seg_end, i, seg[indx, 10][1], seg_length), nr = 1), quote = F, col.names = F, row.names = F, file = paste(CIRCOS, "/gain_focal_clust_inputMC1.4.txt", sep = ""), append = T)
	}))
	posit = t(sapply(cl, function(i) {
		indx = which(seg[, 6] == i & !is.na(seg[, 1]) & seg[, 10] < 2.5)
		if(length(indx) == 0) return(NULL)
		seg_start <- min(seg[indx, 2])
		seg_end <- max(seg[indx, 3])
			#seg_start = round(quantile(seg[indx, 2], 0.25) / 1E5) / 10
			#seg_end = round(quantile(seg[indx, 3], 0.75) / 1E5) / 10
		seg_length = (seg_end - seg_start) / 1E6
		write.table(matrix(c(seg[1, 1], seg_start, seg_end, i, seg[indx, 10][1], seg_length), nr = 1), quote = F, col.names = F, row.names = F, file = paste(CIRCOS, "/gain_focal_clustL2.5_inputMC1.4.txt", sep = ""), append = T)
	}))
	posit = t(sapply(cl, function(i) {
		indx = which(seg[, 6] == i & !is.na(seg[, 1]) & seg[, 10] >= 2.5)
		if(length(indx) == 0) return(NULL)
		seg_start <- min(seg[indx, 2])
		seg_end <- max(seg[indx, 3])
		#seg_start = round(quantile(seg[indx, 2], 0.25))
		#seg_end = round(quantile(seg[indx, 3], 0.75))
		seg_length = (seg_end - seg_start) / 1E6
		write.table(matrix(c(seg[1, 1], seg_start, seg_end, i, seg[indx, 10][1], seg_length), nr = 1), quote = F, col.names = F, row.names = F, file = paste(CIRCOS, "/gain_focal_clustU2.5_inputMC1.4.txt", sep = ""), append = T)
	}))

	ans = seg[, c(1:3, 6, 10)]
	write.table(ans, quote = F, row.names = F, col.names = F, file = paste(CIRCOS, "/gain_focal_inputMC1.4.txt", sep = ""), append = T)
	indx = which(seg[, 10] < 2.5)
	if(length(indx) == 0) {
		write.table(ans, quote = F, row.names = F, col.names = F, file = paste(CIRCOS, "/gain_focal_more2.5_inputMC1.4.txt", sep = ""), append = T)
	} else {
		write.table(ans[indx, ], quote = F, row.names = F, col.names = F, file = paste(CIRCOS, "/gain_focal_less2.5_inputMC1.4.txt", sep = ""), append = T)
		write.table(ans[-indx, ], quote = F, row.names = F, col.names = F, file = paste(CIRCOS, "/gain_focal_more2.5_inputMC1.4.txt", sep = ""), append = T)
	}
})

sapply(ai.files, function(f) {
	cat(f, "\n")
	seg = try(read.table(f, stringsAsFactor = F, head = F))
	if (class(seg) == "try-error") return(NULL)
	seg[, 1] = sub("chr", "hs", seg[, 1])
	cl <- sort(unique(seg[, 6]))
	posit = t(sapply(cl, function(i) {
		indx = which(seg[, 6] == i & !is.na(seg[, 1]))
		seg_start = round(quantile(seg[indx, 2], 0.05))
		seg_end = round(quantile(seg[indx, 3], 0.95))
		#write.table(matrix(c(seg[1, 1], seg_start, seg_end, i), nr = 1), quote = F, col.names = F, row.names = F, file = paste(CIRCOS, "/AI_clust_input.txt", sep = ""), append = T)
	}))
	posit = t(sapply(cl, function(i) {
		indx = which(seg[, 6] == i & !is.na(seg[, 1]) & seg[, 10] < 2.3)
		if(length(indx) == 0) return(NULL)
		seg_start = round(quantile(seg[indx, 2], 0.05))
		seg_end = round(quantile(seg[indx, 3], 0.95))
		write.table(matrix(c(seg[1, 1], seg_start, seg_end, i), nr = 1), quote = F, col.names = F, row.names = F, file = paste(CIRCOS, "/AI_clustL2.3_input.txt", sep = ""), append = T)
	}))
	posit = t(sapply(cl, function(i) {
		indx = which(seg[, 6] == i & !is.na(seg[, 1]) & seg[, 10] >= 2.3)
		if(length(indx) == 0) return(NULL)
		seg_start = round(quantile(seg[indx, 2], 0.05))
		seg_end = round(quantile(seg[indx, 3], 0.95))
		write.table(matrix(c(seg[1, 1], seg_start, seg_end, i), nr = 1), quote = F, col.names = F, row.names = F, file = paste(CIRCOS, "/AI_clustU2.3_input.txt", sep = ""), append = T)
	}))

	ans = seg[, c(1:3, 6)]
	#write.table(ans, quote = F, row.names = F, col.names = F, file = paste(CIRCOS, "/AI_input.txt", sep = ""), append = T)
	indx = which(seg[, 10] < 2.3)
	if(length(indx) == 0) {
		write.table(ans, quote = F, row.names = F, col.names = F, file = paste(CIRCOS, "/AI_more2.3_input.txt", sep = ""), append = T)
	} else {
		write.table(ans[indx, ], quote = F, row.names = F, col.names = F, file = paste(CIRCOS, "/AI_less2.3_input.txt", sep = ""), append = T)
		write.table(ans[-indx, ], quote = F, row.names = F, col.names = F, file = paste(CIRCOS, "/AI_more2.3_input.txt", sep = ""), append = T)
	}
})

sapply(ai.focal, function(f) {
	cat(f, "\n")
	seg = try(read.table(f, stringsAsFactor = F, head = F))
	if (class(seg) == "try-error") return(NULL)
	#seg[, 1] = sub("chr", "hs", seg[, 1])
	cl <- sort(unique(seg[, 6]))
	posit = t(sapply(cl, function(i) {
		indx = which(seg[, 6] == i & !is.na(seg[, 1]))
		seg_start <- min(seg[indx, 2])
		seg_end <- max(seg[indx, 3])
			#seg_start = round(quantile(seg[indx, 2], 0.25))
			#seg_end = round(quantile(seg[indx, 3], 0.75))
		seg_length = (seg_end - seg_start) / 1E6
		write.table(matrix(c(seg[1, 1], seg_start, seg_end, i, seg[indx, 10][1], seg_length), nr = 1), quote = F, col.names = F, row.names = F, file = paste(CIRCOS, "/AI_focal_clust_inputMC1.4.txt", sep = ""), append = T)
	}))
	posit = t(sapply(cl, function(i) {
		indx = which(seg[, 6] == i & !is.na(seg[, 1]) & seg[, 10] < 2.5)
		if(length(indx) == 0) return(NULL)
		seg_start <- min(seg[indx, 2])
		seg_end <- max(seg[indx, 3])
			#seg_start = round(quantile(seg[indx, 2], 0.25) / 1E5) / 10
			#seg_end = round(quantile(seg[indx, 3], 0.75) / 1E5) / 10
		seg_length = (seg_end - seg_start) / 1E6
		write.table(matrix(c(seg[1, 1], seg_start, seg_end, i, seg[indx, 10][1], seg_length), nr = 1), quote = F, col.names = F, row.names = F, file = paste(CIRCOS, "/AI_focal_clustL2.5_inputMC1.4.txt", sep = ""), append = T)
	}))
	posit = t(sapply(cl, function(i) {
		indx = which(seg[, 6] == i & !is.na(seg[, 1]) & seg[, 10] >= 2.5)
		if(length(indx) == 0) return(NULL)
		seg_start <- min(seg[indx, 2])
		seg_end <- max(seg[indx, 3])
		#seg_start = round(quantile(seg[indx, 2], 0.25))
		#seg_end = round(quantile(seg[indx, 3], 0.75))
		seg_length = (seg_end - seg_start) / 1E6
		write.table(matrix(c(seg[1, 1], seg_start, seg_end, i, seg[indx, 10][1], seg_length), nr = 1), quote = F, col.names = F, row.names = F, file = paste(CIRCOS, "/AI_focal_clustU2.5_inputMC1.4.txt", sep = ""), append = T)
	}))

	ans = seg[, c(1:3, 6, 10)]
	write.table(ans, quote = F, row.names = F, col.names = F, file = paste(CIRCOS, "/AI_focal_inputMC1.4.txt", sep = ""), append = T)
	indx = which(seg[, 10] < 2.5)
	if(length(indx) == 0) {
		write.table(ans, quote = F, row.names = F, col.names = F, file = paste(CIRCOS, "/AI_focal_more2.5_inputMC1.4.txt", sep = ""), append = T)
	} else {
		write.table(ans[indx, ], quote = F, row.names = F, col.names = F, file = paste(CIRCOS, "/AI_focal_less2.5_inputMC1.4.txt", sep = ""), append = T)
		write.table(ans[-indx, ], quote = F, row.names = F, col.names = F, file = paste(CIRCOS, "/AI_focal_more2.5_inputMC1.4.txt", sep = ""), append = T)
	}
})


##############################plot circos using OmicCircos
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
#ai.clust[, 1] = gsub("chr", "", ai.clust[, 1])
pdf(paste(CIRCOS, "AI_focal_clust_seglength1.4.pdf", sep = "/"))
par(las = 1, mfrow = c(3, 1))
segleng = ai.clust[, 6]
hist(segleng,  freq = F, probability = T, right = F,  xlab = "Seg_length", mlossn = "All seglength")
hist(segleng[which(segleng < 10)], xlim = c(0, 10), freq = F, probability = T, right = F, xlab = "Seg_length", mlossn = "Seglength < 10MB")
hist(segleng[which(segleng >= 10)], xlim = c(10, max(segleng)), freq = F, probability = T, right = F, xlab = "Seg_length", main = "Seglength >= 10MB")
dev.off()

indx = which(ai.clust[, 6] < 10)
ai.clustSL = ai.clust[indx, 1:6]
ai.clustSM = ai.clust[-indx, 1:6]
ai.clust = ai.clust[, 1:4]

ai.clustL = read.table(paste(CIRCOS, "AI_focal_clustL2.5_inputMC1.4.txt", sep = "/"), header = F, sep = " ", stringsAsFactors = F)
#ai.clustL[, 1] = gsub("chr", "", ai.clustL[, 1])
indx = which(ai.clustL[, 6] < 10)
ai.clustLL = ai.clustL[indx, ]
ai.clustLM = ai.clustL[-indx, ]
ai.clustL = ai.clustL[, 1:4]

ai.clustU = read.table(paste(CIRCOS, "AI_focal_clustU2.5_inputMC1.4.txt", sep = "/"), header = F, sep = " ", stringsAsFactors = F)
#ai.clustU[, 1] = gsub("chr", "", ai.clustU[, 1])
indx = which(ai.clustU[, 6] < 10)
ai.clustUL = ai.clustU[indx, ]
ai.clustUM = ai.clustU[-indx, ]
ai.clustU = ai.clustU[, 1:4]

out.file <- paste(CIRCOS, "AI_focal_clust_circosMC1.4.pdf", sep = "/")
pdf(file = out.file, height = 8, width = 8, compress = TRUE)
par(las = 1, mar = c(0.5, 0.5, 1, 0.5))
plot(c(1,800), c(1,800), type= "n", axes = FALSE, xlab = "", ylab = "")
circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
circos(R = 10, cir = seg.c, type= "arc", W = 380, mapping = ai.clust, B = TRUE, col = colors[c(1)], lwd = 5, col.v = 4)

plot(c(1, 800), c(1, 800), type= "n", axes = FALSE, xlab = "", ylab = "", main = "Cutoff Seg_length = 10MB", cex.main = 0.8)
circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
circos(R = 300, cir = seg.c, type = "arc", W = 90, mapping = ai.clustSM, B = TRUE, col = "purple", lwd = 5, col.v = 6)
circos(R = 20, cir = seg.c, type = "arc", W = 275, mapping = ai.clustSL, B = T, col = "orange", lwd = 5, col.v = 6)

plot(c(1, 800), c(1, 800), type= "n", axes = FALSE, xlab = "", ylab = "", main = "Cutoff DI = 2.5", cex.main = 0.8)
circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
circos(R = 190, cir = seg.c, type = "arc", W = 200, mapping = ai.clustU, B = TRUE, col = "purple", lwd = 5, col.v = 4)
circos(R = 10, cir = seg.c, type = "arc", W = 180, mapping = ai.clustL, B = T, col = "orange", lwd = 5, col.v = 4)
dev.off()

################loss_focal and gain_focal circos
loss.clust = read.table(paste(CIRCOS, "loss_focal_clust_inputMC1.4.txt", sep = "/"), header = F, sep = " ", stringsAsFactors = F)
#loss.clust[, 1] = gsub("chr", "", loss.clust[, 1])
pdf(paste(CIRCOS, "loss_focal_clust_seglength1.4.pdf", sep = "/"))
par(las = 1, mfrow = c(3, 1))
segleng = loss.clust[, 6]
hist(segleng,  freq = F, probability = T, right = F,  xlab = "Seg_length", main = "All seglength")
hist(segleng[which(segleng < 10)], xlim = c(0, 10), freq = F, probability = T, right = F, xlab = "Seg_length", main = "Seglength < 10MB")
hist(segleng[which(segleng >= 10)], xlim = c(10, max(segleng)), freq = F, probability = T, right = F, xlab = "Seg_length", main = "Seglength >= 10MB")
dev.off()

indx = which(loss.clust[, 6] < 10)
loss.clustSL = loss.clust[indx, 1:6]
loss.clustSM = loss.clust[-indx, 1:6]
loss.clust = loss.clust[, 1:4]

loss.clustL = read.table(paste(CIRCOS, "loss_focal_clustL2.5_inputMC1.4.txt", sep = "/"), header = F, sep = " ", stringsAsFactors = F)
#loss.clustL[, 1] = gsub("chr", "", loss.clustL[, 1])
indx = which(loss.clustL[, 6] < 10)
loss.clustLL = loss.clustL[indx, ]
loss.clustLM = loss.clustL[-indx, ]
loss.clustL = loss.clustL[, 1:4]

loss.clustU = read.table(paste(CIRCOS, "loss_focal_clustU2.5_inputMC1.4.txt", sep = "/"), header = F, sep = " ", stringsAsFactors = F)
#loss.clustU[, 1] = gsub("chr", "", loss.clustU[, 1])
indx = which(loss.clustU[, 6] < 10)
loss.clustUL = loss.clustU[indx, ]
loss.clustUM = loss.clustU[-indx, ]
loss.clustU = loss.clustU[, 1:4]

out.file <- paste(CIRCOS, "loss_focal_clust_circosMC1.4.pdf", sep = "/")
pdf(file = out.file, height = 8, width = 8, compress = TRUE)
par(las = 1, mar = c(0.5, 0.5, 1, 0.5))
plot(c(1,800), c(1,800), type= "n", axes = FALSE, xlab = "", ylab = "")
circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
circos(R = 10, cir = seg.c, type= "arc", W = 380, mapping = loss.clust, B = TRUE, col = colors[c(1)], lwd = 5, col.v = 4)

plot(c(1, 800), c(1, 800), type= "n", axes = FALSE, xlab = "", ylab = "", main = "Cutoff Seg_length = 10MB", cex.main = 0.8)
circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
circos(R = 300, cir = seg.c, type = "arc", W = 90, mapping = loss.clustSM, B = TRUE, col = "purple", lwd = 5, col.v = 6)
circos(R = 20, cir = seg.c, type = "arc", W = 275, mapping = loss.clustSL, B = T, col = "orange", lwd = 5, col.v = 6)

plot(c(1, 800), c(1, 800), type= "n", axes = FALSE, xlab = "", ylab = "", main = "Cutoff DI = 2.5", cex.main = 0.8)
circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
circos(R = 190, cir = seg.c, type = "arc", W = 200, mapping = loss.clustU, B = TRUE, col = "purple", lwd = 5, col.v = 4)
circos(R = 10, cir = seg.c, type = "arc", W = 180, mapping = loss.clustL, B = T, col = "orange", lwd = 5, col.v = 4)
dev.off()

###########gene in each cluster
gene.dat = read.table("/share/data0/reference/Genome/hg19.RefSeqGene.uniq.bed", header = F, sep = "\t", stringsAsFactors = F)
ai.clust = read.table(paste(CIRCOS, "AI_focal_clust_inputMC1.4.txt", sep = "/"), header = F, sep = " ", stringsAsFactors = F)
sapply(1:nrow(ai.clust), function(x) {
	indx = which(gene.dat[, 1] == ai.clust[x, 1] & gene.dat[, 2] >= ai.clust[x, 2] & gene.dat[, 3] <= ai.clust[x, 3])
	gene = paste(gene.dat[indx, 4], collapse = ",")
	write.table(matrix(c(ai.clust[x, 1:4], gene), nr = 1), quote = F, col.names = F, row.names = F, sep = "\t", append = T,
		file = paste(CIRCOS, "AI_focal_clust_gene.bed", sep = "/"))
})

#################snp in each cluster
load("map_SNP6_1KG_dbSNP138.RData")
load(paste(WORKDIR, PROCDIR, "allele.RData", sep = "/"))
load(paste(WORKDIR, PROCDIR, "genotype_TR.RData", sep = "/"))
colnames(genotype) <- sub("-\\w+-\\w+-\\w+-\\w+$", "", colnames(genotype))
sample_id <- sub("-\\w+-\\w+-\\w+-\\w+$", "", colnames(a))
SNP = SNPs[grep("^SNP", SNPs$Affy_Probe_ID, value = F), ]
ai.clust = read.table(paste(CIRCOS, "AI_focal_clust_inputMC1.4.txt", sep = "/"), header = F, sep = " ", stringsAsFactors = F)

bands = unique(ai.clust[, 1])
aiclust_snp = sapply(bands, function(j) {
	clust = ai.clust[which(ai.clust[, 1] == j), ]
	
	aisnp = NULL
	for(x in 1:nrow(clust)) {
		indx = which(SNP$Chromsome_hg19 == clust[x, 1] & SNP$Physical_Position_hg19 >= clust[x, 2] & SNP$Physical_Position_hg19 <= clust[x, 3])
		snp = SNP$Affy_Probe_ID[indx]
		rs = SNP$dbSNP_RS_ID[indx]
		chr = SNP$Chromsome_hg19[indx]
		position = unlist(SNP$Physical_Position_hg19[indx])
		strand = SNP$Strand[indx]
		ref = unlist(SNP$Allele_A_1KG[indx])
		alt = unlist(SNP$Allele_B_1KG[indx])
		
		abdat = NULL
		for(i in unique(sample_id)) {
			cat(paste(j, x, i, sep = "_"), "\n")
			if (!i %in% colnames(genotype)) next
			indx.n <- which(sample_id == i & (barc2type(colnames(a)) == "NB" | barc2type(colnames(a)) == "NT"))[1]
			indx.t <- which(sample_id == i & barc2type(colnames(a)) == "TP")[1]
			if (is.na(indx.n)) next
			if (is.na(indx.t)) next
			pair_i <- c(indx.n, indx.t) 
			snpdat <- array(cbind(a[, pair_i], b[, pair_i]), dim = c(nrow(a), 2, 2), dimnames = list(rownames(a), colnames(a)[pair_i], c("A", "B")))
			snpid = match(snp, dimnames(snpdat)[[1]])
			ab.n = paste(snpdat[snpid, 1, 1], snpdat[snpid, 1, 2], sep = "/")
			ab.t = paste(snpdat[snpid, 2, 1], snpdat[snpid, 2, 2], sep = "/")
			ans = cbind(ab.n, ab.t)
			colnames(ans) = paste(gsub("\\w-\\w+-\\w+-\\w+$", "", dimnames(snpdat)[[2]]), "A/B", sep = "_")
			abdat = cbind(abdat, ans)	
		}
		bns = cbind(chr = chr, clust_start = clust[x, 2], clust_end = clust[x, 3], clust = clust[x, 4], rs = rs, position = position, strand = strand, ref = ref, alt = alt, abdat)
		aisnp = rbind(aisnp, bns)
	}
	save(aisnp, file = paste(CIRCOS, paste(j, "_AI_clust_snp.RData", sep = ""), sep = "/"))
	gc()		
	aisnp
})
save(aiclust_snp, file = paste(CIRCOS, "AI_clust_snp.RData", sep = "/"))

#########################genotype in each clust, not used
load("/share/data4/TCGA/SNP6_Genotype/00_gtypRData/gtyp_LUAD.RData")
gtyp = gtyp.normal.luad
SNPchg = rownames(SNPs)[which(SNPs$flipped == "TRUE")]
odx = which(rownames(gtyp) %in% SNPchg)
gtyp[odx, ] = 2 - gtyp[odx, ]
colnames(gtyp) =  sub("-\\w+-\\w+-\\w+-\\w+$", "", colnames(gtyp))
SNP = SNPs[grep("^SNP", SNPs$Affy_Probe_ID, value = F), ]

sample_id <- read.table(paste(FILTDIR, "/", "sample_id.txt", sep = ""), stringsAsFactor = F)[, 1]

loh.focal = read.csv(paste0(FILTDIR, "/loh_focalMC.csv"), header = T, stringsAsFactor = F)
colnames(loh.focal) = gsub("\\.", "-", colnames(loh.focal))
ai.clust = read.table(paste(CIRCOS, "AI_focal_clust_inputMC.txt", sep = "/"), header = F, sep = " ", stringsAsFactors = F)
loh.sample = cbind(ai.clust[, 1:4], loh.focal[, -c(1:6)])

bands = unique(loh.sample[, 1])
sapply(bands[1:2], function(i) {
	lohs = loh.sample[which(loh.sample[, 1] == i), ]
	write.table(matrix(c("Chr", "Start", "End", "Cluster", "Affy_Probe_ID", "RS_ID", "Physical_Position_hg19", sample_id), nr = 1), 
		quote = F, col.names = F, row.names = F, 
		file = paste(CIRCOS, paste(i, "SNPgtype_loh_focalMC.bed", sep = "_"), sep = "/"), sep = "\t")
	
	for(j in 1: nrow(lohs)) {
		x = lohs[j, ]
		x = unlist(x)
		print(x[1:4])
		indx = which(SNP$Chromsome_hg19 == x[1] & SNP$Physical_Position_hg19 >= x[2] & SNP$Physical_Position_hg19 <= x[3])
		snp = SNP$Affy_Probe_ID[indx]
		rs = SNP$dbSNP_RS_ID[indx]
		position = SNP$Physical_Position_hg19[indx]
		gtype = matrix(NA, length(snp), length(sample_id))
		rownames(gtype) = snp
		colnames(gtype) = sample_id
		id = sample_id[which(x[-c(1:4)] == 1)]
		coldx = match(id, colnames(gtyp))
		rowdx = match(snp, rownames(gtyp))
		gtype[, id] = gtyp[rowdx, coldx]
		ans = cbind(matrix(x[1:4], length(snp), 4, byrow = T), snp, rs, position, gtype)
		write.table(ans, quote = F, col.names = F, row.names = F, 
			file = paste(CIRCOS, paste(i, "SNPgtype_loh_focalMC.bed", sep = "_"), sep = "/"), sep = "\t", append = T)
	}
})
############################
