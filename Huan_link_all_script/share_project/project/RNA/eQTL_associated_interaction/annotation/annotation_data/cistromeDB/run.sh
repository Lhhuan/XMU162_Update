# the peak calling is from the https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html
# the file header is following the picture
perl 01_merge.pl #将./Human_FACTOR/human_factor, ./HISTONE_MARK_AND_VARIANT/human_hm, ./Human_CHROMATIN_Accessibility/human_ca下的文件merge到一起"${output_dir}/merge_pos_info_sample_narrow_peak.bed.gz"


perl 03_filter_chr1_22.pl
#过滤出${output_dir}/merge_pos_info_narrow_peak_sorted.bed.gz中的chr1_22,得${output_dir}/03_merge_pos_info_narrow_peak_sorted_chr1_22.bed.gz
#将${output_dir}/03_merge_pos_info_narrow_peak_sorted_chr1_22.bed.gz左右各扩500bp得${output_dir}/03_merge_pos_info_narrow_peak_sorted_chr1_22_extend_500bp.bed.gz
#将${output_dir}/03_merge_pos_info_narrow_peak_sorted_chr1_22.bed.gz左右各扩100bp得${output_dir}/03_merge_pos_info_narrow_peak_sorted_chr1_22_extend_100bp.bed.gz



#-----------normal cell line concents
perl 04_class_cell_line.pl #筛选@files = ("./Human_FACTOR/human_factor_full_QC.txt","./HISTONE_MARK_AND_VARIANT/human_hm_full_QC.txt","./Human_CHROMATIN_Accessibility/human_ca_full_QC.txt") 中在"/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/CCLE/CCLE_sample_info_file_2012-10-18.txt"的文件，得存在的./cell_line_info/04_existing_ccle_cell_line_info.txt,得不存在但是cancer cell line ./cell_line_info/04_not_exist_ccle_but_cancer_cell_line_unique.txt，得附加原始信息文件 ./cell_line_info/04_not_exist_ccle_but_cancer_cell_line_unique_info.txt
#得不在ccle中存在，从名字中无法判断cancer cell line的 unique cell line 得./cell_line_info/04_unique_cell_line_without_info.txt 得附加原始信息文件 ./cell_line_info/04_cell_line_without_info_original.txt
cat ./cell_line_info/04_unique_cell_line_without_info.txt | sort > ./cell_line_info/04_unique_cell_line_without_info_sort.txt

cat ./cell_line_info/04_cell_line_without_info_original.txt | awk 'NR >1'| sort -k1 >./cell_line_info/04_cell_line_without_info_original_sort.txt

#手动为./cell_line_info/04_unique_cell_line_without_info_sort.txt 寻找细胞信息得./cell_line_info/04_unique_cell_line_without_info_sort_mannual_find_info.txt
# cat unique_cell_line_without_info.txt | sort > unique_cell_line_without_info_sort.txt
perl 05_filter_normal_chrom_info.pl # #筛选@files = ("./Human_FACTOR/human_factor_full_QC.txt","./HISTONE_MARK_AND_VARIANT/human_hm_full_QC.txt","./Human_CHROMATIN_Accessibility/human_ca_full_QC.txt") 中在"./cell_line_info/04_unique_cell_line_without_info_sort_mannual_find_info.txt"的cell lien，并在相应文件夹提取出文件得${output_dir}/merge_pos_info_sample_narrow_peak.bed.gz，得对于文件及peak信息文件"${output_dir}/merge_pos_info_narrow_peak.bed.gz"

perl 061_union_segment.pl #将normal cell ${dir}/merge_pos_info_narrow_peak_sort.bed.gz bedtools merge -i 为${dir}/merge_pos_info_narrow_peak_sort_union.bed.gz，得互补文件${dir}/merge_pos_info_narrow_peak_union_complement.bed.gz

# Rscript 06_plot_length_distribution_of_factor_and_calculate_mean.R
# perl 07_get_no_factor_and_split.pl 



# #------------------------------
# perl 061_union_segment.pl #将${dir}/merge_pos_info_narrow_peak_sort.bed.gz bedtools merge -i 为${dir}/merge_pos_info_narrow_peak_sort_union.bed.gz，得互补文件${dir}/merge_pos_info_narrow_peak_union_complement.bed.gz
# Rscript 06_plot_length_distribution_of_factor_and_calculate_mean_union.R

#----------------------------
#cancer.txt from lichengxin
perl 08_filter_TCGA_marker.pl #