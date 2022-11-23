wget -c https://epd.epfl.ch/ucsc/epdHub/hg19/bigWIG/barski2007/PolII.bw #CD4+ T cells
bigWigToWig PolII.bw PolII.wig
wig2bed < PolII.wig > PolII.bed
liftOver PolII.bed "/home/huanhuan/reference/hg19ToHg38.over.chain.gz" PolII_hg38.bed unmap.bed
gzip PolII.bed
perl 01_filter_chr1_22.pl  #PolII_hg38_chr1_22.bedgraph
less PolII_hg38_chr1_22.bedgraph |sort -k1,1 -k2,2n >PolII_hg38_chr1_22_sorted.bedgraph
less PolII_hg38_chr1_22_sorted.bedgraph |cut -f1-3 |gzip >PolII_hg38_chr1_22_sorted.bed.gz
bedtools merge -i PolII_hg38_chr1_22_sorted.bed.gz |gzip > PolII_hg38_chr1_22_sorted_merge.bed.gz
bedtools intersect -a PolII_hg38_chr1_22_sorted_merge.bed.gz -b PolII_hg38_chr1_22_sorted.bedgraph -wa -wb |gzip >PolII_hg38_chr1_22_sorted_merge_bed_signal.bed.gz


bedGraphToBigWig PolII_hg38_chr1_22.bedgraph "/share/Projects/huanhuan/ref_data/gencode/GRCh38_chrom.sizes" PolII_hg38_CD4_Tcell.bw 