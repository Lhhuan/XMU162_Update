wget -c https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E096-H3K27ac.narrowPeak.gz
wget -c https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E096-H3K4me1.narrowPeak.gz
#---------------
perl 01_normal.pl #对"E096-H3K27ac.narrowPeak.gz","E096-H3K4me1.narrowPeak.gz" 进行normal和排序得${fo1}sorted.gz
# bedtools intersect -F 0.1 -a 01_normal_E096-H3K4me1.narrowPeaksorted.gz -b 01_normal_E096-H3K27ac.narrowPeaksorted.gz -wo |gzip  > 01_normal_E096-H3K4me1_H3K27ac.narrowPeaksorted.gz
bedtools intersect -F 0.1 -a 01_normal_E096-H3K4me1.narrowPeaksorted.gz -b 01_normal_E096-H3K27ac.narrowPeaksorted.gz  |sort -k1,1 -k2,2n | gzip  > 01_normal_E096-H3K4me1_H3K27ac.narrowPeaksorted_overlap.gz
bedtools  complement -i 01_normal_E096-H3K4me1_H3K27ac.narrowPeaksorted_overlap.gz -g /home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt | gzip> 01_enhancer_complement.bed.gz

#---------------------------------
Rscript 02_plot_length_distribution_of_enhancer_and_calculate_mean.R #为01_normal_E096-H3K4me1_H3K27ac.narrowPeaksorted_overlap.gz 画长度分布图，得长度文件 02_enhancer_length.txt 和求均数及中位数得02_enhancer_length_statistics.txt

bedtools makewindows -b  01_enhancer_complement.bed.gz -w 238 | gzip >split_enhancer_complement.bed.gz
# perl 03_filter_chr1_22_enhancer_complement.pl  #将01_enhancer_complement.bed中的chr1 和chr22提取出来，得01_enhancer_complement_chr1_22.bed
# bedtools makewindows -g  01_enhancer_complement_chr1_22.bed -w 334 | gzip >split_enhancer_complement_chr1_22.bed.gz
ln 01_normal_E096-H3K4me1_H3K27ac.narrowPeaksorted_overlap.gz /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/Lung/enhancer.bed.gz
ln 01_enhancer_complement.bed.gz /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/Lung/non_enhancer.bed.gz
ln split_enhancer_complement.bed.gz /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/Lung/non_enhancer_split.bed.gz

# cd split_complement_chr1_22_files
# zless ../split_enhancer_complement_chr1_22.bed.gz | split -l 1000000 -a 5







