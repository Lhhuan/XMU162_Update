# perl 01_transform_varint_ID_hg38_to_hg19.pl #利用"../data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz" 将../data/GTEx_Analysis_v8_eQTL/*.v8.signif_variant_gene_pairs.txt.gz转换为hg19, 得../data/GTEx_Analysis_v8_eQTL/*.v8.signif_variant_gene_pairs.txt.gz
perl 02_Completion_snp_for_xQTL_by_1kg.pl # # 用"/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz" 补全"${dir}/${tissue}${suffix}"; 得"../output/Pancreas_cis_eQTL_1kg_Completion.txt.gz"
Rscript 03_NHP_big_par.R
Rscript 04_histgram_density_interval.R
Rscript 05_histgram_density_interval_chr_specific.R
#-------------------以hotspot 为基准计数
# perl 06_filter_hotspot_for_interval18.pl ###../../output/${tissue}/Cis_eQTL/NHP/NHPoisson_emplambda_interval_${j}_cutoff_7.3_${tissue}.txt.gz 时得不同cutoff下的hotspot(segment),"../../output/${tissue}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${tissue}_segment_hotspot_cutoff_${cutoff}.bed.gz"
# perl 06_filter_non_hotspot_for_interval18.pl ###../output/Pancreas/Cis_eQTL/NHP/NHPoisson_emplambda_interval_${j}_cutoff_7.3_Pancreas.txt.gz 时得不同cutoff下的non hotspot(segment),../output/${tissue}/Cis_eQTL/non_hotspot_cis_eQTL/interval_${j}/${tissue}_segment_non_hotspot_cutoff_${cutoff}.bed.gz
# #----------------------
# perl 061_count_hotspot_in_different_cutoff.pl #统计不同cutoff下hotspot的数目/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Pancreas/Cis_eQTL/hotspot_cis_eQTL/interval_18/，得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Pancreas/Cis_eQTL/061_hotspot_in_different_cutoff_interval_18.txt


# Rscript 062_histgram_plot_the_length_of_hotspot.R #看各种cutoff 下面hotspot长度的分布
# # perl 07_annotation_hotspot_interval_18_different_cutoff.pl #interval_18 时，对../output/Pancreas/Cis_eQTL/hotspot_cis_eQTL/interval_18/Pancreas_segment_hotspot_cutoff_${cutoff}.bed.gz用annotation_bedtools_intersect_interval18.sh进行annotation,得$output_dir/$factor_$input_file_base_name

# Rscript 08_merge_cis_eQTL_hotspot_annotation_interval18.R

# Rscript 08_test.R

# perl 09_prepare_for_roc.pl ##看各种ROC准备数据，利用"~/project/RNA/eQTL_associated_interaction/GTEx/output/Pancreas/Cis_eQTL/annotation_out/interval_18/ALL/hotspot/08_annotation_merge_interval_18_overlap_cutoff_${cutoff}.txt.gz" "~/project/RNA/eQTL_associated_interaction/GTEx/output/Pancreas/Cis_eQTL/annotation_out/interval_18/ALL/non_hotspot/08_annotation_merge_interval_18_overlap_cutoff_${cutoff}.txt.gz"计算roc曲线的准备数据得，
# #得"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Pancreas/Cis_eQTL/ROC/interval_18/ALL/09_prepare_number_ROC.txt"
# perl 091_ratio_of_factor_in_hotspot.pl ###看各种factor在"~/project/RNA/eQTL_associated_interaction/GTEx/output/Pancreas/Cis_eQTL/annotation_out/interval_18/ALL/hotspot/08_annotation_merge_interval_18_overlap_cutoff_${cutoff}.txt.gz"中hotspot 在factor中的比例得"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Pancreas/Cis_eQTL/ROC/interval_18/ALL/091_factor_ratio_of_hotspot.txt"
# Rscript 10_fisher_exact_test.R
# Rscript 10_plot_point.R
# #-------------------------以factor 为基准计数
# perl 07_annotation_factor_and_non_factor_interval_18_different_cutoff.pl #interval_18 时，对../output/Pancreas/Cis_eQTL/hotspot_cis_eQTL/interval_18/Pancreas_segment_hotspot_cutoff_${cutoff}.bed.gz用annotation_bedtools_intersect_interval18.sh进行annotation,得$output_dir/$factor_$input_file_base_name
# perl 08_count_factor_and_non_factor_annotation_interval18_different_cutoff.pl #根据/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/${factor}.bed.gz， /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/non_${factor}.bed.gz，../output/Pancreas/Cis_eQTL/annotation/interval_18/ALL/factor/hotspot/${cutoff}/${factor}_Pancreas_segment_hotspot_cutoff_${cutoff}.bed.gz，../output/Pancreas/Cis_eQTL/annotation/interval_18/ALL/non_factor/hotspot/${cutoff}/non_${factor}_Pancreas_segment_hotspot_cutoff_${cutoff}.bed.gz 以factor为基数，准备计算ROC的四格表得"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Pancreas/Cis_eQTL/ROC/interval_18/ALL/08_prepare_number_ROC_refine_factor_count.txt"

# Rscript 10_plot_point_factor_non_factor.R
# Rscript 10_fisher_exact_test_factor_and_non_factor.R

# #---------------------------------#不同fraction 下进行bedtools intersect 
# perl 07_annotation_factor_and_non_factor_interval_18_different_cutoff_fraction.pl #interval_18 时，对../output/Pancreas/Cis_eQTL/hotspot_cis_eQTL/interval_18/Pancreas_segment_hotspot_cutoff_${cutoff}.bed.gz用annotation_bedtools_intersect_interval18.sh进行annotation,得$output_dir/$factor_$input_file_base_name
# perl 08_count_factor_and_non_factor_annotation_interval18_different_cutoff_fraction.pl #根据/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/${factor}.bed.gz， /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/non_${factor}.bed.gz，../output/Pancreas/Cis_eQTL/annotation/interval_18/ALL/factor/hotspot/${cutoff}/${fraction}/${factor}_Pancreas_segment_hotspot_cutoff_${cutoff}.bed.gz，../output/Pancreas/Cis_eQTL/annotation/interval_18/ALL/non_factor/hotspot/${cutoff}/${fraction}/non_${factor}_Pancreas_segment_hotspot_cutoff_${cutoff}.bed.gz 以factor为基数，准备计算ROC的四格表得"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Pancreas/Cis_eQTL/ROC/interval_18/ALL/08_prepare_number_ROC_refine_factor_count_fraction.txt"
# #------------------------------



# #--------
# perl 09_filter_significant_eQTL.pl #

# #------------------------更细的cutoff

perl 06_filter_hotspot_for_interval18_re.pl
perl 061_count_hotspot_in_different_cutoff_re.pl #统计不同cutoff下hotspot的数目/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Pancreas/Cis_eQTL/hotspot_cis_eQTL/interval_18/，得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Pancreas/Cis_eQTL/061_hotspot_in_different_cutoff_interval_18_re.txt
Rscript 062_histgram_plot_the_length_of_hotspot_re.R

#----以hotspot数目为基准
perl 20_count_hotspot_in_factor.pl ##根据"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/061_hotspot_in_different_cutoff_interval_18_re.txt记录的hotspot数目，和annotation信息得每个cutoff下hit住factor的hotspot的比例文件"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/20_${tissue}_count_hotspot_in_factor.txt"; 得每个cutoff下被每种factor hit住hotspot的比例文件"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/20_${tissue}_count_hotspot_in_per_factor_ratio.txt";得每个cutoff下被每种factor hit住hotspot的文件"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/20_count_hotspot_in_factor_pos.txt.gz"

perl 21_NHP_background_intersect_hotspot.pl #将"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/${tissue}_cis_eQTL_1kg_Completion.txt.gz” 转换为bed得 "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/${tissue}_cis_eQTL_1kg_Completion.bed.gz,排序得$sort_fo1，用$sort_fo1 和"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/${tissue}_cis_eQTL_1kg_Completion.bed.gz bedtools intersect得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/hotspot_cis_eQTL_overlap_SNP/${tissue}_segment_hotspot_cutoff_${cutoff}_SNP.bed.gz

perl 22_count_eqtl_in_hotspot_ratio.pl ##统计/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/hotspot_cis_eQTL_overlap_SNP/${tissue}_segment_hotspot_cutoff_${cutoff}_SNP.bed.gz中eQTL在SNP中的比例得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/${tissue}_cis_eQTL_in_hotspot_ratio.gz