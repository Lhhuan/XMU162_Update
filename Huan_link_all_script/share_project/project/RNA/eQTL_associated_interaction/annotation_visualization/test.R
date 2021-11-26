library(GenomicRanges)
library()
library(chipseq)


data("cstest")
ctcfReads <- cstest$ctcf
ctcfFragments <- resize(ctcfReads, 120)
ctcfCoverage <- coverage(ctcfFragments)
ctcfCoverage10 <- ctcfCoverage$chr10
maxPos <- which.max(ctcfCoverage10)


roi <- resize(IRanges(maxPos, width = 1), 5000,"center")
roiCoverage<- ctcfCoverage$chr10[roi]
ctcfPeaks <- slice(ctcfCoverage10, lower = 8)

ctcfMaxs <- viewMaxs(ctcfPeaks)

library(leeBamViews)
library(GenomicFeatures)
bams <- getBamsFromLeeBamViews()


gr <- GRanges(
  seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),#查看基因的名称
  ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),#查看基因的区域
  strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),#查看基因位于哪条链
  score = 1:10,
  GC = seq(1, 0, length=10))
gr

library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
exon_txdb=exons(txdb)
tmp=as.data.frame(exon_txdb)

gene <-as.data.frame(genes(txdb))


chr.rle <- Rle(values = c("Chr1", "Chr2", "Chr3", "Chr1", "Chr3"), lengths = c(3,
    2, 5, 4, 5))
chr.rle