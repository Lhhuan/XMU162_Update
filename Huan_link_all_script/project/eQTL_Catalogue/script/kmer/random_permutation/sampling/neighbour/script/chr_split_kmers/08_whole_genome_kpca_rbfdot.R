library(unix)
rlimit_as(1e18)
library(kernlab)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure")
load("08_permutation_wilocx_overlap_sig_kmer_0_1000.Rdata")
pca1<- kpca(~.,data=Sorg, kernel="rbfdot",kpar = list(sigma = 0.1),features = 0, th = 1e-4)