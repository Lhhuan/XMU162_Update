perl 32_random_genomic_resemble_hotspot.pl  ##产生10000个与../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/whole_blood_segment_${group}_cutoff_${cutoff}.bed.gz相同的resemble hotspot,"$output_dir/${i}_resemble_${input_file_base_name}",
perl 33_random_permutation_hotspot.pl # 对 resemble hotspot,"$output_dir/${i}_resemble_${input_file_base_name} 和../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/whole_blood_segment_${group}_cutoff_${cutoff}.bed.gz 进行 random_permutation，得得../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_random/random_permutation/${cutoff}/${i}_resemble_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz
perl 34_annotation_chrom_state_original_random_and_permutation_hotspot.pl #对1000个original_random和random_permutation进行chromatin_states的注释，得../../../output/${tissue}/Cis_eQTL/annotation/interval_18/ALL/${group}/${cutoff}/${type}/${i}_25_state_resemble_${input_file_base_name}
perl 35_filter_annotation_chrom_state.pl #因为每个片段会被多个chrom state 注释，根据对hotspot的覆盖比例，选出覆盖比例最高的state, interval_18 时，对"../../../output/${tissue}/Cis_eQTL/annotation/interval_18/ALL/${group}/${cutoff}/${type}/${state}_state_resemble_${i}_resemble_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz"进行过滤，得$dir/filter_$name
perl 36_count_per_state_number.pl ##对"$dir/filter_$name"中的chromatin state 进行统计，得"../../../output/${tissue}/Cis_eQTL/enrichment/interval_18/ALL/${type}_${state}_state_count_${tissue}_cutoff_${cutoff}.txt"
#-------------------
mv /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/annotation/interval_18/ALL/hotspot/0.176/random_permutation  /share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/interval_18/annotation/ALL/hotspot/0.176/
mv /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/annotation/interval_18/ALL/hotspot/0.176/original_random  /share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/interval_18/annotation/ALL/hotspot/0.176/

#-------plot figure
Rsript 37_boxplot_plot_compare_hotspot_and_background_15state.R 
Rsript 37_boxplot_plot_compare_hotspot_and_background_15state_count.R
Rsript 37_boxplot_plot_compare_hotspot_and_background_15state_percentage.R
Rsript 37_boxplot_plot_compare_hotspot_and_random_permutation_15state_percentage.R
Rsript 37_boxplot_plot_compare_hotspot_and_emp0_15state_percentage.R
Rsript 37_boxplot_plot_compare_hotspot_and_emp0_0.176_15state_percentage.R
perl 38_transform_36.pl 
Rsript 39_boxplot_plot_compare_hotspot_and_background_25state.R 

#------------histone marker
perl 40_annotation_histone_marker.pl #对10000个original_random histone的注释，得"/share/data0/QTLbase/huan/GTEx/${tissue}/Cis_eQTL/interval_18/annotation/ALL/${group}/${cutoff}/random_permutation/

perl 41_count_annotation_marker.pl ##对/share/data0/QTLbase/huan/GTEx/${tissue}/Cis_eQTL/interval_18/annotation/ALL/${group}/${cutoff}/$type/*的marker进行count,得"$out_dir/${type}_histone_marker.txt.gz";


#------------histone marker
Rsript 42_boxplot_plot_compare_hotspot_and_original_random_histone_marker_count.R
Rsript 42_boxplot_plot_compare_hotspot_and_original_random_histone_marker_percentage.R
Rsript 42_boxplot_plot_compare_hotspot_and_0_histone_marker.R 
Rsript 42_boxplot_plot_compare_hotspot_and_0_0.176_histone_marker_count.R
#-----------------------------------------------------------------------------


perl 43_calculate_jaccard_index_original_random.pl # #对$input_dir/${mark}_${i}_resemble_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz 进行jaccard index 进行计算,得的marker进行计算,得"$out_dir/${type}_jaccard_index_histone_marker.txt.gz";

Rscript 44_boxplot_compare_hotspot_and_0_0.176_histone_marker_jaccard_index.R
Rscript 44_boxplot_compare_hotspot_and_0_0.176_histone_marker_jaccard_index_single.R
Rscript 44_boxplot_compare_hotspot_and_orginal_random_histone_marker_jaccard_index.R
Rscript 44_boxplot_compare_hotspot_and_orginal_random_histone_marker_jaccard_index_single.R
Rscript 44_boxplot_compare_hotspot_and_0_histone_marker_jaccard_index.R

#---------------------------
Rscript 45_compare_hotspot_and_0_0.176_fisher_test.R
Rscript 45_compare_hotspot_and_0_fisher_test.R
Rscript 45_compare_hotspot_and_original_random_fisher_test.R
Rscript 45_compare_hotspot_and_0_0.176_fisher_test_histone_and_tfbs.R
Rscript 45_compare_hotspot_and_0_fisher_test_histone_and_tfbs.R
Rscript 45_compare_hotspot_and_original_random_fisher_test_histone_and_tfbs.R
#------------------------------
Rscript 46_geom_pointrange_0.R
Rscript 46_geom_pointrange_0_log_fdr.R
Rscript 46_geom_pointrange_0_log_fdr_histone_tfbs.R
Rscript 46_geom_pointrange_0_0.176.R
Rscript 46_geom_pointrange_0_0.176_log_fdr.R
Rscript 46_geom_pointrange_0_0.176_log_fdr_histone_tfbs.R
Rscript 46_geom_pointrange_original_random.R
Rscript 46_geom_pointrange_original_random_log_fdr_histone_tfbs.R

cp 45_compare_hotspot_and_original_random_fisher_test.R 45_compare_hotspot_and_original_random_fisher_test_histone_and_tfbs.R





