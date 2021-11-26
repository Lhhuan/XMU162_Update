wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp151Common.txt.gz

wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz

gunzip hg19.refGene.gtf.gz

wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes 
#
perl filter_chr1_22_size.pl  #过滤hg19.chrom.sizes中chr1_22得hg19.chrom1_22.sizes
sort -k1,1n hg19.chrom1_22.sizes >hg19.chrom1_22_sizes_sorted.txt

