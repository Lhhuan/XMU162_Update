wget -c https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
perl filter_chr1_22_size.pl
sort -k1,1  hg38.chrom1_22.sizes >hg38.chrom1_22_sizes_sorted.txt