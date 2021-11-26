perl  07_annotation_non_factor_split_interval_18_different_cutoff.pl #interval_18 时，对."/share/data0/QTLbase/huan/GTEx/random_select/${tissue}/${group}/${number}/interval_18/${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz"用annotation_${type}_${tissue}_bedtools_intersect_interval18.sh进行annotation,得${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz
perl 08_count_factor_and_non_factor_annotation_interval18_different_cutoff.pl #根据/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/${factor}.bed.gz， /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/non_${factor}.bed.gz，"/share/data0/QTLbase/huan/GTEx/random_select/${tissue}/annotation/interval_18/ALL/${number}/factor/${group}/${cutoff}/${factor}_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz"，"/share/data0/QTLbase/huan/GTEx/random_select/${tissue}/annotation/interval_18/ALL/${number}/non_factor_split/${group}/${cutoff}/non_${factor}_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz" 以factor为基数，准备计算ROC的四格表得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample_exclude_eQTL/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/08_prepare_number_ROC_factor_count.txt

Rscript 10_plot_point_factor_non_factor.R
Rscript 10_fisher_exact_test_factor_and_non_factor.R
Rscript 11_F-score.R

#--------------------------
perl 081_sum_random_result.pl #将/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample_exclude_eQTL/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/08_prepare_number_ROC_factor_count.txt中random的数据合在一起，得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample_exclude_eQTL/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/081_sum_random_result.txt
Rscript 10_plot_point_factor_non_factor_sum.R
Rscript 10_fisher_exact_test_factor_and_non_factor_summary.R

#-----
perl 11_merge_random_and_original_prepare_number_ROC_factor_count.pl #将#/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample_exclude_eQTL/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/081_sum_random_result.txt和
# /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/08_prepare_number_ROC_factor_count.txt 
# /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/08_prepare_number_ROC_factor_count.txt merge 在一起
#得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample_exclude_eQTL/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/11_merge_random_and_original_prepare_number_ROC_factor_count.txt
perl 12_merge_random_and_original_fisher_exact_test.pl  #将#/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample_exclude_eQTL/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/10_fisher_exact_test_result_factor_summary.txt和
# /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/10_fisher_exact_test_result_factor.txt 
# /home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/10_fisher_exact_test_result_factor.txt merge 在一起
#得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample_exclude_eQTL/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/11_merge_random_and_original_fisher_exact_test.txt