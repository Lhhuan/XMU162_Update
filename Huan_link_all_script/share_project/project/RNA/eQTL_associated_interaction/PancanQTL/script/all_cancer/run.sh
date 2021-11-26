perl 01_merge_five_race_1kg.pl 
perl 02_Completion_snp_for_PancanQTL_by_1kg.pl
Rscript 03_NHP_big_par.R
perl 04_filter_hotspot_for_interval18.pl
perl 05_bedtools_cancer_normal_total_share_hotspot.pl # # 相同组织，cancer 和 normal -f 1 intersect ，找全部share的hotspot,得my $out1 = "../../output/cancer_total/share/total/${cancer}_contain_${tissue}.bed";my $out2 = "../../output/cancer_total/share/total/${tissue}_contain_${cancer}.bed";
#得汇总文件"../../output/cancer_total/share/total/05_cancer_tissue_intersect_total.bed.gz", cancer 被完全包含的文件"../../output/cancer_total/share/total/05_cancer_share_total.bed.gz"

perl 06_bedtools_cancer_normal_all_overlap.pl # 相同组织，cancer 和 normal  intersect, 找全部overlap(为找cancer specifc做铺垫)的hotspot,得my $out1 = "../../output/cancer_total/share/all/${cancer}_contain_${tissue}.bed";
#得汇总文件"../../output/cancer_total/share/all/05_cancer_tissue_intersect_all.bed.gz",

perl 07_filter_cancer_and_normal_specific.pl #利用../../data/pancanQTL_gtex_eQTL.txt，../../output/cancer_total/share/all/${cancer}_contain_${tissue}.bed，"../../output/${cancer}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${cancer}_segment_hotspot_cutoff_${cutoff}.bed.gz";"../../../GTEx/output/${tissue}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${tissue11}_segment_hotspot_cutoff_${cutoff}.bed.gz" 鉴别组织特异的hotspot"$output_dir/tissue_${cancer}_${tissue}_specific.bed.gz"和cancer特异的hotspot  $output_dir/cancer_${cancer}_${tissue}_specific.bed.gz,cancer 特异汇总文件为../../output/cancer_total/specific/pure/cancer_specific.bed.gz，tissue特异汇总文件../../output/cancer_total/specific/pure/tissue_specific.bed.gz
perl 08_count_cancer_all_and_specific_hotspot.pl #得"../../output/cancer_total/08_number_of_cancer_all_hotspot.txt.gz"，"../../output/cancer_total/08_number_of_cancer_specific_hotspot.txt.gz"

perl 09_count_cancer_related_tissue_specific_hotspot.pl #"../../output/cancer_total/09_count_cancer_related_tissue_specific_hotspot.txt.gz"

Rscript 10_bar_plot_cancer_tissue_specific_and_all_hotspot.R

Rscript 10_bar_plot_cancer_specific_and_all_hotspot.R

perl 11_count_cancer_tissue_share_hotspot.pl  #利用../../output/cancer_total/share/total/05_cancer_share_total.bed.gz得#../../output/cancer_total/11_count_cancer_tissue_share_hotspot.txt.gz
# Rscript 11_cancer_tissue_share_hotspot.R
Rscript 12_barplot_Cleveland_distbution_of_hotspot_in_share_tissues.R 
Rscript 13_Cleveland_dot_plot.R
perl 14_extract_cancer_qtl.pl  #对../../data/cis_eQTLs_all_re.gz进行整理，得cancer 对应的文件"../../output/cancer_total/specific/pure/${cancer}/qtl_gene.bed.gz";
perl 15_find_gene_affected_by_cancer_specific_hotspot.pl  #找 hotspot对应的eGene $output_dir/cancer_${cancer}_${tissue}_specific.bed.gz 和"$output_dir/qtl_gene.bed.gz" 得 "$output_dir/cancer_${cancer}_${tissue}_specific_gene.bed.gz"
#得汇总文件"../../output/cancer_total/15_cancer_specfic_hotspot_gene.txt.gz"

Rscript 16_KEGG_and_go.R
perl 17_bedtools_cancer_mutual_absolute.pl #cancer 两两之间进行 intersect ，找全部share的hotspot,得汇总文件"../../output/cancar_total/share/total/17_absolute_cancer_cancer_intersect.bed.gz"
perl 18_cancer_cancer_share_hotspot.pl #对 "../../output/cancer_total/share/total/17_absolute_cancer_cancer_intersect.bed.gz" 进行 tissue合并得"../../output/cancer_total/share/total/18_all_tissue_cancer_hotspot_total_contain.bed.gz"
perl 19_count_share_cancer_number.pl #对"../../output/cancer_total/share/total/18_all_tissue_cancer_hotspot_total_contain.bed.gz"进行统计,得"../../output/cancer_total/19_share_cancer_number_count.txt.gz";
Rscript 20_barplot_distbution_of_hotspot_in_share_cancers.R
perl 21_count_jaccard_index_similarity_cancer_cancer.pl #得100%绝对overlap文件"../../output/cancer_total/21_cancer_pair_ovelap_absolute.txt.gz"，得没有在"../../output/cancer_total/21_cancer_pair_ovelap_absolute.txt.gz"中的片段"../../output/cancer_total/21_cancer_pair_out_ovelap_absolute.txt.gz"
#得用于计算类似jaccard index 数据"../../output/cancer_total/21_cancer_pair_overlap_index_and_related_number.txt.gz"，得类似jaccard index "../../output/cancer_total/21_cancer_pair_overlap_index.txt.gz"

Rscript 22_heatmap_cancer_tissue_share.R

Rscript 24_pathway_enrichment_plot.R
perl 25_annotation_marker.pl  #annotation marker ,得"/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cancer_cell/HISTONE_MARK_AND_VARIANT/${cancer}/${marker}/merge_pos_info_narrow_peak_sort_union_sort.bed.gz";

bedtools intersect -a ../../output/KIRP/Cis_eQTL/hotspot_cis_eQTL/interval_18/KIRP_segment_hotspot_cutoff_0.176.bed.gz  -b merge_pos_info_narrow_peak_sort_union_sort.bed |gzip > ../../output/KIRP/Cis_eQTL/anno/CHROMATIN_Accessibility_KIRP_segment_hotspot_cutoff_18.bed.gz

perl 26_calculate_jaccard_index.pl #"../../output/${cancer}/Cis_eQTL/anno/${cancer}_segment_hotspot_cutoff_${cutoff}_marker_jaccard_index.txt.gz";
perl 27_annotation_share_and_specific.pl  #"../../output/${cancer}/Cis_eQTL/anno/${cancer}_segment_hotspot_cutoff_${cutoff}_marker_jaccard_index_lable.txt.gz";
Rscript 28_boxplot_marker_jaccard_index.R 
Rscript 29_compare_specific_and_share_marker_fisher_test.R
Rscript 30_geom_pointrange_specific_and_share_marker_fisher_test.R



perl 31_bedtools_cancer_mutual_absolute.pl  ##cancer两两之间进行 intersect ，找全部share的hotspot,得汇总文件"../../output/cancer_total/share/total/31_absolute_cancer_intersect.bed.gz"
perl 32_merge_global_tissue_share_hotspot.pl ##对 "../../output/cancer_total/share/total/31_absolute_cancer_intersect.bed.gz" 进行 tissue合并得../../output/cancer_total/share/total/32_all_tissue_share_hotspot_total_contain.bed.gz

perl 33_count_share_tissue_number.pl #统计 #"../../output/cancer_total/share/total/32_all_tissue_share_hotspot_total_contain.bed.gz"得"../../output/cancer_total/share/total/33_share_cancer_number_count.txt.gz"

Rscript 34_barplot_distbution_of_hotspot_in_share_tissues.R
perl 35_1_extract_max_cancer_share_hotspot.pl
perl 35_2_extract_eqtl_egene.pl
Rscript 35_3_intersect_gene.R
"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_total/11_3_intersect_gene.R"