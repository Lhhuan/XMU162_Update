perl 04_extract_cell_name.pl #提取../data/download.sh 中的cell line信息，得./unique_cell_line.txt
#利用"/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cell_line_info/04_existing_ccle_cell_line_info.txt"， "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cell_line_info/04_not_exist_ccle_but_cancer_cell_line_unique.txt";和"/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cell_line_info/04_unique_cell_line_without_info_sort_mannual_find_info.txt"提供cell line 的disease,得有disease ./04_cell_line_in_cistrome.txt,得没有disease文件 ./04_cell_line_without_cistrome.txt

#手动将./04_cell_line_without_cistrome.txt annotation 得./04_cell_line_without_cistrome_manual.txt

cat ./04_cell_line_in_cistrome.txt ./04_cell_line_without_cistrome_manual.txt  > 04_all_cell_line_info.txt

perl 05_filter_normal_cell_line.pl #将04_all_cell_line_info.txt 中对应的normal cell line 从/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/CTCF/data/中选出来，得05_normal_cell_line_ctcf.bed.gz

zless 05_normal_cell_line_ctcf.bed.gz | sort -k1,1 -k2,2n |gzip >05_normal_cell_line_ctcf_sort.bed.gz

#-------------------------
bedtools merge -i 05_normal_cell_line_ctcf_sort.bed.gz | gzip > 05_normal_cell_line_ctcf_sort_union.bed.gz

liftOver 05_normal_cell_line_ctcf_sort_union.bed.gz "/home/huanhuan/reference/hg19ToHg38.over.chain.gz"  ./hg38/05_normal_cell_line_ctcf_sort_union.bed 05_normal_cell_line_ctcf_sort_union_unmap_hg38.bed
gzip ./hg38/05_normal_cell_line_ctcf_sort_union.bed  05_normal_cell_line_ctcf_sort_union_unmap_hg38.bed



zless 05_normal_cell_line_ctcf_sort_union.bed.gz |sort -k1,1 -k2,2n |gzip >05_normal_cell_line_ctcf_sort_union_sort.bed.gz
bedtools merge -i 05_normal_cell_line_ctcf_sort_union_sort.bed.gz | gzip > 05_normal_cell_line_ctcf_sort_union_merge.bed.gz
#------------------------------------------------------unused
#-----------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------



# bedtools complement -i 05_normal_cell_line_ctcf_sort.bed.gz -g /home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt | gzip> 05_CTCF_complement.bed.gz

# Rscript 06_plot_length_distribution_of_CTCF_mean.R #为01_normal_E062-H3K4me1_H3K27ac.narrowPeaksorted_overlap.gz 画长度分布图，得长度文件 02_enhancer_length.txt 和求均数及中位数得02_enhancer_length_statistics.txt
# ln 05_normal_cell_line_ctcf_sort.bed.gz "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/Whole_blood/CTCF.bed.gz"

# ln split_CTCF_complement.bed.gz "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/Whole_blood/non_CTCF_split.bed.gz"

# Rscript 06_plot_length_distribution_of_CTCF_mean_union.R
