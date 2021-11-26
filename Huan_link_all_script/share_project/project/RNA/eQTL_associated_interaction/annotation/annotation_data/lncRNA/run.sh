wget -c http://www.noncode.org/datadownload/NONCODEv5_hg19.lncAndGene.bed.gz
perl 01_normalized.pl #将NONCODEv5_hg19.lncAndGene.bed.gz中需要的列提出来，得NONCODEv5_hg19.lncAndGene_need.bed.gz
zless NONCODEv5_hg19.lncAndGene_need.bed.gz | sort -k1,1 -k2,2n | gzip >NONCODEv5_hg19.lncAndGene_sorted.bed.gz
perl 03_filter_chr1_22.pl ##过滤出NONCODEv5_hg19.lncAndGene_sorted.bed.gz中的chr1_22,得03_NONCODEv5_hg19.lncAndGene_sorted_chr1_22.bed.gz
#将03_fantom5_enhancers_phase1_phase2_sorted_chr1_22.bed.gz左右各扩500bp得03_NONCODEv5_hg19.lncAndGene_sorted_chr1_22_extend_500bp.bed.gz
#将03_fantom5_enhancers_phase1_phase2_sorted_chr1_22.bed.gz左右各扩100bp得03_NONCODEv5_hg19.lncAndGene_sorted_chr1_22_extend_100bp.bed.gz