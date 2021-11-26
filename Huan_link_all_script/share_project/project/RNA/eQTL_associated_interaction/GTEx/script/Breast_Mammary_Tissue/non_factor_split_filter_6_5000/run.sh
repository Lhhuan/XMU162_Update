perl 06_filter_segment_6_5000.pl #过滤../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz中6-5000 bp的片段

perl 061_count_hotspot_in_different_cutoff.pl #过滤../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz中6-5000 bp的片段,得../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_filter/6_5000/${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz



# perl 07_annotation_non_factor_split_interval_18_different_cutoff.pl ##interval_18 时，对../../../output/Lung/Cis_eQTL/hotspot_cis_eQTL/interval_18/Lung_segment_hotspot_cutoff_${cutoff}.bed.gz用annotation_bedtools_intersect_interval18.sh进行annotation,得$output_dir/$factor_$input_file_base_name
# perl 08_count_factor_and_non_factor_annotation_interval18_different_cutoff.pl #根据/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/${factor}.bed.gz， /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/non_${factor}.bed.gz，../../../output/Lung/Cis_eQTL/annotation/interval_18/ALL/factor/hotspot/${cutoff}/${factor}_Lung_segment_hotspot_cutoff_${cutoff}.bed.gz，../../../output/Lung/Cis_eQTL/annotation/interval_18/ALL/non_factor_split/hotspot/${cutoff}/non_${factor}_Lung_segment_hotspot_cutoff_${cutoff}.bed.gz 以factor为基数，准备计算ROC的四格表得"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/ROC/interval_18/ALL/08_prepare_number_ROC_refine_factor_count.txt"

# Rscript 10_plot_point_factor_non_factor.R
# Rscript 10_fisher_exact_test_factor_and_non_factor.R
# perl 11_F-score.pl


#--------------------------更细的cutoff,from=0.1,to=0.3,by=0.01

perl  07_annotation_non_factor_split_interval_18_different_cutoff_re.pl
perl 08_count_factor_and_non_factor_annotation_interval18_different_cutoff_re.pl

perl  20_count_hotspot_in_factor_filter_6_5000.pl #根据"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/ROC/interval_18/ALL/split_non_factor_filter/6_5000/061_hotspot_in_different_cutoff_interval_18.txtt记录的hotspot数目，和annotation信息得每个cutoff下hit住factor的hotspot的比例文件"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/ROC/interval_18/ALL/split_non_factor_filter/6_5000/20_${tissue}_count_hotspot_in_factor.txt.gz"; 得每个cutoff下被每种factor hit住hotspot的比例文件"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/ROC/interval_18/ALL/split_non_factor_filter/6_5000/20_${tissue}_count_hotspot_in_per_factor_ratio.txt.gz";得每个cutoff下被每种factor hit住hotspot的文件"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Lung/Cis_eQTL/ROC/interval_18/ALL/split_non_factor_filter/6_5000/20_${tissue}_count_hotspot_in_factor_pos.txt.gz"