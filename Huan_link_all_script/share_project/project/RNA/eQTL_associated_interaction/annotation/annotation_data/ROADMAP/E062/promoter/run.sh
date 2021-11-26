wget -c https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E062-H3K4me3.narrowPeak.gz
perl 01_normal.pl  ##对"E062-H3K4me3.narrowPeak.gz 进行normal和排序得${fo1}sorted.gz
bedtools complement -i 01_normal_E062-H3K4me3.narrowPeaksorted.gz -g /home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt | gzip> 01_promoter_complement.bed.gz
cp 01_normal_E062-H3K4me3.narrowPeaksorted.gz /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/promoter.bed.gz
cp 01_promoter_complement.bed.gz /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/non_promoter.bed.gz



#----------------------
Rscript 02_plot_length_distribution_of_promoter_and_calculate_mean.R #为01_normal_E062-H3K4me1_H3K27ac.narrowPeaksorted_overlap.gz 画长度分布图，得长度文件 02_enhancer_length.txt 和求均数及中位数得02_enhancer_length_statistics.txt


gzip 01_promoter_complement.bed.gz

bedtools makewindows -b  01_promoter_complement.bed.gz -w 1007 | gzip >split_promoter_complement.bed.gz


cp split_promoter_complement.bed.gz /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/non_promoter_split.bed.gz