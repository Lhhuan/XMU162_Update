perl 07_annotation_non_factor_split_interval_18_different_cutoff.pl ##interval_18 时，对../../../output/Liver/Cis_eQTL/hotspot_cis_eQTL/interval_18/Liver_segment_hotspot_cutoff_${cutoff}.bed.gz用annotation_bedtools_intersect_interval18.sh进行annotation,得$output_dir/$factor_$input_file_base_name
perl 08_count_factor_and_non_factor_annotation_interval18_different_cutoff.pl #根据/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/${factor}.bed.gz， /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/non_${factor}.bed.gz，../../../output/Liver/Cis_eQTL/annotation/interval_18/ALL/factor/hotspot/${cutoff}/${factor}_Liver_segment_hotspot_cutoff_${cutoff}.bed.gz，../../../output/Liver/Cis_eQTL/annotation/interval_18/ALL/non_factor_split/hotspot/${cutoff}/non_${factor}_Liver_segment_hotspot_cutoff_${cutoff}.bed.gz 以factor为基数，准备计算ROC的四格表得"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Liver/Cis_eQTL/ROC/interval_18/ALL/08_prepare_number_ROC_refine_factor_count.txt"
perl 09_merge_qtl_factor_ratio_and_hotspot_factor_ratio.pl
Rscript 10_plot_point_factor_non_factor.R
Rscript 10_fisher_exact_test_factor_and_non_factor.R
perl 11_F-score.pl


#--------------------------更细的cutoff,from=0.1,to=0.3,by=0.01

perl  07_annotation_non_factor_split_interval_18_different_cutoff_re.pl
perl 08_count_factor_and_non_factor_annotation_interval18_different_cutoff_re.pl

perl 09_merge_qtl_factor_ratio_and_hotspot_factor_ratio_re.pl #将/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Liver/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/08_Liver_prepare_number_ROC_factor_count_re.txt和



