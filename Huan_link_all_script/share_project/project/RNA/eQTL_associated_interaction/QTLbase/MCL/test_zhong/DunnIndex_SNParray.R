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
DUNNINDEX <- ARGS["dunnindex", 1]

setwd(WORKDIR)

library(rARPACK)
library(cluster)
library(parallel)
library(fastcluster)
library(clv)
#library(RWeka)
library(matrixLaplacian)
#library(sparcl)
#library(SNFtool)
library(fpc)
library(MCL)
library(igraph)

#source("xmeans_function.R")

#load(ANNOTATION)
load(CYTOBAND)

ancestry <- read.csv(ANCESTRY, head = T, stringsAsFactor = F, row.names = 1)
rownames(ancestry) = sub("-\\w+-\\w+-\\w+-\\w+$", "", rownames(ancestry))
#sample_id <- rownames(ancestry)
sample_id <- read.table(paste(FILTDIR, "/", "sample_id.txt", sep = ""), stringsAsFactor = F)[, 1]

readEmat <- function(file) {
	x <- readLines(file)
	n <- length(x)
	ans <- matrix(NA, n, n)
	for(i in 1:n) {
		e <- as.numeric(strsplit(x[i], split = ",")[[1]])
		ans[i, 1:i] <- e[1:i]
		ans[1:i, i] <- e[1:i]
	}
	ans	
}
	
spec.clust <- function(Emat, method = "Laplacian", p = 2, clust = "hierarchical", ...) {
	Dmat <- diag(apply(Emat, 1, sum))
	if (method == "Laplacian") {
		Lmat <- Dmat - Emat
	} else if (method == "normalised.Laplacian") {
		D_minus_sqrt <- diag(1 / sqrt(diag(Dmat)))
		Lmat <- diag(nrow(Emat)) - D_minus_sqrt %*% Emat %*% D_minus_sqrt
	}
	eig <- eigs(Lmat, which = "SM", k = p, opts = list(maxitr = 50000))
	if (clust == "hierarchical") {
		cl <- cutree(hclust(dist(eig$vectors), method = "complete"), ...)
	}
	if (clust == "kmeans") {
		cl <- kmeans(eig$vectors, ...)$cluster
	}
	if (clust == "skmeans") {
		cl <- skmeans(eig$vectors, ...)$cluster
	}
	if (clust == "pam") {
		cl <- pam(eig$vectors, ...)$cluster
	}
	cl
}

files <- list.files(path = FILTDIR, pattern = "*.bed", full = T)
loss.focal = grep("loss_focal", files, value = T)
gain.focal = grep("gain_focal", files, value = T)
ai.focal = grep("AI_focal", files, value = T)

seg_out = sapply(c(gain.focal), function(f) {
	cat(f, "\n")
	seg <- read.table(f, stringsAsFactor = F, head = T)
	if(nrow(seg) == 0) return(NULL)
	system(paste0("perl ", WORKDIR, "/seg_overlap.pl < ", f, " > ", sub(".bed", ".emat", f)))
	Emat <- readEmat(sub(".bed", ".emat", f))
	if(class(Emat) != "matrix") return(NULL)
	out = rep(NA, nrow(Emat))
	for(i in 1:nrow(Emat)) {
		e <- as.numeric(Emat[i, ])
		out[i] = ifelse(all(e[-i] == 0) & e[i] == 1, i, NA)
		#rate = length(which(e[-i] != 0)) / nrow(Emat)
		#out[i] = ifelse(rate <= 0.1, i, NA)
	}
	outid = which(!is.na(out))
	if(length(outid) != 0) {
		seg = seg[-outid, ]
		Emat = Emat[-outid, -outid]
	}
	#Dmat <- diag(apply(Emat, 1, sum))
	#D_minus_sqrt <- diag(1 / sqrt(diag(Dmat)))
	#Lmat <- diag(nrow(Emat)) - D_minus_sqrt %*% Emat %*% D_minus_sqrt
	#eig <- eigs(Lmat, which = "SM", k = 3, opts = list(maxitr = 500000))
	
		#NL <- matrixLaplacian(Emat, plot2D = F, plot3D = F)
		#Lmat <- NL$LaplacianMatrix
	#eig <- eigs(Lmat, which = "SM", k = 20)
		#eig <- eigen(Lmat)$values
	#eig = matrixLaplacian(Dmat)
	#cl <- kmeans(eig$vectors[, which(eig < 1E-10)], centers = sqrt(nrow(seg) / 2))$cluster
	#cl <- xmeans(eig$vectors[, which(eig$values < 1E-10)], ik = 2)$cluster
	#cl <- xmeans(NL$eigenvector[, which(eig < 1E-10)], ik = sqrt(nrow(seg) / 2), iter.max = 40)$cluster
	#km.perm <- KMeansSparseCluster.permute(eig$vectors, K = sqrt(nrow(seg) / 2))
	#km.out <- KMeansSparseCluster(eig$vectors, K = sqrt(nrow(seg) / 2), wbounds = km.perm$bestw,  maxiter = 100)
	#cl <- km.out[[1]]$Cs
	
	#aft = affinityMatrix(Emat, K = 20, alpha = 0.5)
	#cl = spectralClustering(aft, K = 20, type = 3)
		#p = ifelse((nrow(seg) - 1) < 80, nrow(seg) - 1, 80)
		#n = ncol(NL$eigenvector)
	#cl <- pamk(NL$eigenvector[, which(eig < 1E-10)], krange = sqrt(nrow(seg) / 2) : p)[[1]]$clustering
		#nc <- pamk(NL$eigenvector[, (n - 2) : n], krange = sqrt(nrow(seg) / 2) : p)$nc
		#cl <- pamk(NL$eigenvector[, (n - nc + 1) : n], krange = nc)[[1]]$clustering
	mat = Emat
	mat[which(mat != 0)] = 1	
	nc <- try(mcl(mat, addLoops = T, expansion = 2, inflation = 1.4)) # ai_focal:expansion = 2, inflation = 1.4(THCA_chr18q: inflation = 4; chr9q: inflation = 2; chr1p: Emat, inflation = 6; chr3q: Emat, i = 1.4);
	if(class(nc) == "character") {
		mat = Emat
		mat[which(mat < 0.3)] = 0
		mat[which(mat >= 0.3)] = 1
		nc <- try(mcl(mat, addLoops = T, expansion = 2, inflation = 1.4))
	}
	
	cls = nc$Cluster
	clu = sort(unique(cls))
	cln = 1 : length(clu)
	names(cln) = clu
	cl = cln[as.character(cls)]
	cls.scatt = cls.scatt.diss.mx(1 - Emat, cl)
	gc()
	
	ans <- t(sapply(1:max(cl), function(x) {
		indx <- which(cl == x)
		N <- length(unique(seg[indx, 4]))
		m <- Emat[indx, indx]
		S <- sum(m[lower.tri(m)], na.rm = T) / sum(as.numeric(lower.tri(m)))
		if(max(cl) == 1) DIC = Inf else {
			interD = min(cls.scatt$intercls.average[x, -x], na.rm = T)
			intraD = cls.scatt$intracls.average[x]
			DIC = interD/intraD	
		}
		c(N, S, length(indx), DIC)
	}))
	ans <- cbind(seg[, 1:5], cl, ans[cl, ])
	indx <- which(ans[, 7] > 5)
	write.table(ans[indx, ], quote = F, row.names = F, col.names = F, file = sub(".bed", ".MCclust1.4", f))
	# system(paste("rm ", sub(".bed", ".emat", f), sep = ""))
	gc()
   #c(list(cl), list(cls.scatt), list(ans))
   outid
})

#save(clust_sc, file = paste(FILTDIR, "/DunnIndex.RData", sep = ""))

####################################

files = list.files(path = FILTDIR, pattern = "*.MCclust1.4", full = T)
ai.focal = grep("AI_focal", files, value = T)
del.focal = grep("loss_focal", files, value = T)
amp.focal = grep("gain_focal", files, value = T)
cn.files = c(del.files, amp.files)


focal.cnDIC = as.vector(unlist(sapply(c(del.focal, amp.focal), function(f) {
	seg = try(read.table(f, stringsAsFactor = F, head = F))
	if (class(seg) == "try-error") return(NA)
	dic = seg[, 10]
	dic[which(dic != Inf)]
})))
focal.aiDIC = as.vector(unlist(sapply(ai.focal, function(f) {
	cat(f, "\n")
	seg = try(read.table(f, stringsAsFactor = F, head = F))
	if (class(seg) == "try-error") return(NA)
	dic = seg[, 10]
	dic[which(dic != Inf)]
})))

pdf(paste(FILTDIR, "/DunnIndexMC1.4.pdf", sep = ""))
par(las = 1, mfrow = c(2, 1))
hist(focal.cnDIC, freq = F, right = F, xlab = "Dunn Index", cex.axis = 0.8, cex.lab = 0.8, main = "Cnv_focal_DI")
hist(log(focal.cnDIC), freq = F, right = F, xlab = "log(Dunn Index)", cex.axis = 0.8, cex.lab = 0.8, main = "log(Cnv_focal_DI)")
hist(focal.cnDIC[which(focal.cnDIC < 2.5)], xlim = c(0, 2.5), breaks = 10, freq = F, right = F, xlab = "Dunn Index", cex.axis = 0.8, cex.lab = 0.8, main = "Cnv_focal_DI<2.5")
hist(focal.cnDIC[which(focal.cnDIC >= 2.5)], freq = F, right = F, xlab = "Dunn Index", cex.axis = 0.8, cex.lab = 0.8, main = "Cnv_focal_DI>=2.5")
hist(focal.aiDIC, xlim = c(0, max(focal.aiDIC, na.rm = T)), freq = F, right = F, xlab = "Dunn Index", cex.axis = 0.8, cex.lab = 0.8, main = "AI_focal_DI")
hist(log(focal.aiDIC), freq = F, right = F, xlab = "log(Dunn Index)", cex.axis = 0.8, cex.lab = 0.8, main = "log(AI_focal_DI)")
hist(focal.aiDIC[which(focal.aiDIC < 2.5)], xlim = c(0, 2.5), breaks = 10, freq = F, right = F, xlab = "Dunn Index", cex.axis = 0.8, cex.lab = 0.8, main = "AI_focal_DI<2.5")
hist(focal.aiDIC[which(focal.aiDIC >= 2.5)], freq = F, right = F, xlab = "Dunn Index", cex.axis = 0.8, cex.lab = 0.8, main = "AI_focal_DI>=2.5")
dev.off()

files <- list.files(path = FILTDIR, pattern = "loss_focal.*.\\.MCclust1.4", full = T)
write.table(matrix(c("chr", "cluster", "start", "end", "DI", "event", sample_id), nr = 1), 
		quote = F, col.names = F, row.names = F, 
		file = paste0(FILTDIR, "/loss_focalMC1.4.csv"), sep = ",")
del_focal_rate = sapply(files, function(f) {
	cat(f, "\n")
	seg <- try(read.table(f, stringsAsFactor = F, head = F))
	if (class(seg) == "try-error") return(NULL)
	cl <- sort(unique(seg[, 6]))
	info <- strsplit(sub(".+\\/", "", sub(".MCclust1.4", "", f)), split = "_")[[1]]
	event <- info[1]
	band <- info[3]
	txt <- read.table(paste0(PROCDIR, "/", band, "_CNstate.txt"), stringsAsFactor = F)
	sna <- rep(NA, length(sample_id))
	names(sna) <- sample_id
	sna[unique(txt[!is.na(txt[, 4]), 1])] <- 0
	ans = sapply(cl, function(i) {
		indx <- which(seg[, 6] == i & !is.na(seg[, 1]))
		#if(seg[indx[1], 8] > 0.1) {
			seg_start <- min(round(seg[indx, 2] / 1E5) / 10)
			seg_end <- max(round(seg[indx, 3] / 1E5) / 10)
			#seg_start = round(quantile(seg[indx, 2], 0.05) / 1E5) / 10
			#seg_end = round(quantile(seg[indx, 3], 0.95) / 1E5) / 10
			DI = seg[indx, 10][1]
			event <- "del"
			id <- seg[indx, 4]
			state <- sna			
			state[id] <- -1
			write.table(matrix(c(band, i, seg_start, seg_end, DI, event, state), nr = 1), 
					quote = F, col.names = F, row.names = F, 
					file = paste0(FILTDIR, "/loss_focalMC1.4.csv"), sep = ",", append = T)
			
		#}
		length(unique(seg[indx, 4])) / length(sample_id)	
	})
	names(ans) = cl
	ans
})

files <- list.files(path = FILTDIR, pattern = "gain_focal.*.\\.MCclust1.4", full = T)
write.table(matrix(c("chr", "cluster", "start", "end", "DI", "event", sample_id), nr = 1), 
		quote = F, col.names = F, row.names = F, 
		file = paste0(FILTDIR, "/gain_focalMC1.4.csv"), sep = ",")
amp_focal_rate = sapply(files, function(f) {
	seg <- try(read.table(f, stringsAsFactor = F, head = F))
	if (class(seg) == "try-error") return(NULL)
	cl <- sort(unique(seg[, 6]))
	info <- strsplit(sub(".+\\/", "", sub(".MCclust1.4", "", f)), split = "_")[[1]]
	event <- info[1]
	band <- info[3]
	txt <- read.table(paste0(PROCDIR, "/", band, "_CNstate.txt"), stringsAsFactor = F)
	sna <- rep(NA, length(sample_id))
	names(sna) <- sample_id
	sna[unique(txt[!is.na(txt[, 4]), 1])] <- 0
	ans = sapply(cl, function(i) {
		indx <- which(seg[, 6] == i & !is.na(seg[, 1]))
		#if (seg[indx[1], 8] > 0.1) {
			seg_start <- min(round(seg[indx, 2] / 1E5) / 10)
			seg_end <- max(round(seg[indx, 3] / 1E5) / 10)
			#seg_start = round(quantile(seg[indx, 2], 0.05) / 1E5) / 10
			#seg_end = round(quantile(seg[indx, 3], 0.95) / 1E5) / 10
			DI = seg[indx, 10][1]
			event <- "amp"
			id <- seg[indx, 4]
			state <- sna			
			state[id] <- 1
			write.table(matrix(c(band, i, seg_start, seg_end, DI, event, state), nr = 1), 
					quote = F, col.names = F, row.names = F, 
					file = paste0(FILTDIR, "/gain_focalMC1.4.csv"), sep = ",", append = T)
		#}
		length(unique(seg[indx, 4])) / length(sample_id)
	})
	names(ans) = cl
	ans
})

files <- list.files(path = FILTDIR, pattern = "AI_focal_chr.*.\\.MCclust1.4", full = T)
write.table(matrix(c("chr", "cluster", "start", "end", "DI", "event", sample_id), nr = 1), 
		quote = F, col.names = F, row.names = F, 
		file = paste0(FILTDIR, "/loh_focalMC1.4.csv"), sep = ",")
ai_focal_rate = sapply(files, function(f) {
	cat(f, "\n")
	seg <- try(read.table(f, stringsAsFactor = F, head = F))
	if (class(seg) == "try-error") return(NULL)
	cl <- sort(unique(seg[, 6]))
	info <- strsplit(sub(".+\\/", "", sub(".MCclust1.4", "", f)), split = "_")[[1]]
	event <- info[1]
	band <- info[3]
	txt <- read.table(paste0(PROCDIR, "/", band, "_allelicratio.txt"), stringsAsFactor = F)
	sna <- rep(NA, length(sample_id))
	names(sna) <- sample_id
	sna[unique(txt[!is.na(txt[, 3]), 1])] <- 0
	ans = sapply(cl, function(i) {
		indx <- which(seg[, 6] == i & !is.na(seg[, 1]))
		#if (seg[indx[1], 8] > 0.01) {                         #>0.3
		seg_start <- round(min(seg[indx, 2]) / 1E6, 2)
		seg_end <- round(max(seg[indx, 3]) / 1E6, 2)
		#seg_start = round(quantile(seg[indx, 2], 0.05) / 1E5) / 10
		#seg_end = round(quantile(seg[indx, 3], 0.95) / 1E5) / 10
		DI = seg[indx, 10][1]
		event <- "loh"
		id <- seg[indx, 4]
		state <- sna			
		state[id] <- 1
		write.table(matrix(c(band, i, seg_start, seg_end, DI, event, state), nr = 1), 
					quote = F, col.names = F, row.names = F, 
					file = paste0(FILTDIR, "/loh_focalMC1.4.csv"), sep = ",", append = T)
		#}
		length(unique(seg[indx, 4])) / length(sample_id)
	})
	names(ans) = cl
	ans
})
save(del_focal_rate, amp_focal_rate, ai_focal_rate, file = paste(FILTDIR, "/cl_rateMC1.4.RData", sep = ""))

pdf(paste(FILTDIR, "/clust_incidenceMC1.4.pdf", sep = ""), width = 8.27, height = 12.7)
par(las = 1, mfrow = c(2, 1))
hist(unlist(del_focal_rate), freq = F, probability = T, right = F, xlab = "Rate", labels = T, cex = 0.7, main = "Delition_focal_incidence")
#hist(unlist(del_focal_rate)[which(unlist(del_focal_rate) > 0.05)], freq = F, probability = T, right = F, xlab = "Rate", labels = T, cex = 0.7, main = "Delition_focal_incidence(>0.05)")
hist(unlist(amp_focal_rate), freq = F, probability = T, right = F, xlab = "Rate", labels = T,  cex = 0.7, main = "Amplification_focal_incidence")
#hist(unlist(amp_focal_rate)[which(unlist(amp_focal_rate) > 0.05)], freq = F, probability = T, right = F, xlab = "rate", labels = T, cex = 0.7, main = "Amplification_focal_incidence(>0.05)")
hist(unlist(ai_focal_rate), xlim = c(0, max(unlist(ai_focal_rate))), freq = F, probability = T, right = F, xlab = "Rate", labels = T, cex = 0.7, main = "Allelicratio_focal_incidence")
#hist(unlist(ai_focal_rate)[which(unlist(ai_focal_rate) > 0.05)], xlim = c(0.05, max(unlist(ai_focal_rate))), freq = F, probability = T, right = F, xlab = "Rate", labels = T, cex = 0.7, main = "Allelicratio_focal_incidence(>0.05)")
dev.off()

