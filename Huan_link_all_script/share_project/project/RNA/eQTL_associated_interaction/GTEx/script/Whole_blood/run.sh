perl 01_transform_varint_ID_hg38_to_hg19.pl #利用"../data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz" 将../data/GTEx_Analysis_v8_eQTL/*.v8.signif_variant_gene_pairs.txt.gz转换为hg19, 得../data/GTEx_Analysis_v8_eQTL/*.v8.signif_variant_gene_pairs.txt.gz
perl 02_Completion_snp_for_xQTL_by_1kg.pl # # 用"/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz" 补全"${dir}/${tissue}${suffix}"; 得"../output/Whole_Blood_cis_eQTL_1kg_Completion.txt.gz"
Rscript 03_NHP_big_par.R
Rscript 04_histgram_density_interval.R
Rscript 04_violin_plot_emplamda.R
Rscript 05_histgram_density_interval_chr_specific.R
perl 06_filter_hotspot_for_interval18_chr1_22.pl ######../../output/Whole_Blood/Cis_eQTL/NHP/NHPoisson_emplambda_interval_${j}_cutoff_7.3_Whole_Blood.txt.gz 时得不同cutoff下的hotspot(segment),"../../output/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/whole_blood_segment_hotspot_cutoff_${cutoff}_chr1_22.bed.gz"
#----------------
perl 06_filter_non_hotspot_for_interval18_chr1_22.pl ###../../output/Whole_Blood/Cis_eQTL/NHP/NHPoisson_emplambda_interval_${j}_cutoff_7.3_Whole_Blood.txt.gz 时得不同cutoff下的hotspot(segment),../../output/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/whole_blood_segment_hotspot_cutoff_${cutoff}_chr1_22.bed.gz
#------------------------
perl 07_annotation_hotspot_interval_18_different_cutoff_chr1_22.pl #interval_18 时，对../../output/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18/whole_blood_segment_hotspot_cutoff_${cutoff}_chr1_22.bed.gz用annotation_bedtools_intersect_interval18.sh进行annotation,得$output_dir/$factor_$input_file_base_name
Rscript 08_merge_cis_eQTL_hotspot_annotation_interval18_chr1_22.R
perl 09_prepare_for_roc_chr1_22.pl ##看各种ROC准备数据，利用"~/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/annotation_out/interval_18/ALL/hotspot/08_annotation_merge_interval_18_overlap_cutoff_${cutoff}_chr1_22.txt.gz" "~/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/annotation_out/interval_18/ALL/non_hotspot/08_annotation_merge_interval_18_overlap_cutoff_${cutoff}_chr1_22.txt.gz"计算roc曲线的准备数据得，
#得"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/ALL/09_prepare_number_ROC.txt"
#-------------------以hotspot 为基准计数
perl 06_filter_hotspot_for_interval18.pl ###../../output/Whole_Blood/Cis_eQTL/NHP/NHPoisson_emplambda_interval_${j}_cutoff_7.3_Whole_Blood.txt.gz 时得不同cutoff下的hotspot(segment),"../../output/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/whole_blood_segment_hotspot_cutoff_${cutoff}_chr1_22.bed.gz"
perl 06_filter_non_hotspot_for_interval18.pl ###../../output/Whole_Blood/Cis_eQTL/NHP/NHPoisson_emplambda_interval_${j}_cutoff_7.3_Whole_Blood.txt.gz 时得不同cutoff下的hotspot(segment),../../output/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/whole_blood_segment_hotspot_cutoff_${cutoff}_chr1_22.bed.gz
#----------------------
perl 061_count_hotspot_in_different_cutoff.pl #统计不同cutoff下hotspot的数目/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18/，得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/061_hotspot_in_different_cutoff_interval_18.txt
Rscript 062_histgram_plot_the_length_of_hotspot.R #看各种cutoff 下面hotspot长度的分布
perl 07_annotation_hotspot_interval_18_different_cutoff.pl #interval_18 时，对../../output/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18/whole_blood_segment_hotspot_cutoff_${cutoff}.bed.gz用annotation_bedtools_intersect_interval18.sh进行annotation,得$output_dir/$factor_$input_file_base_name

Rscript 08_merge_cis_eQTL_hotspot_annotation_interval18.R


perl 09_prepare_for_roc.pl ##看各种ROC准备数据，利用"~/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/annotation_out/interval_18/ALL/hotspot/08_annotation_merge_interval_18_overla


Rscript 10_fisher_exact_test.R
Rscript 10_fisher_exact_test_factor_and_non_factor.R
Rscript 10_plot_point.R
Rscript 10_plot_point_factor_non_factor.R

#--------------更细的cutoff
perl 06_filter_hotspot_for_interval18_re.pl
perl 061_count_hotspot_in_different_cutoff_re.pl
062_count_the_length_of_hotspot_re.R


perl 20_count_hotspot_in_factor.pl #根据"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/061_hotspot_in_different_cutoff_interval_18_re.txt记录的hotspot数目，和annotation信息得每个cutoff下hit住factor的hotspot的比例文件"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/20_count_hotspot_in_factor.txt.gz"; 得每个cutoff下被每种factor hit住hotspot的比例文件"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/20_count_hotspot_in_per_factor_ratio.txt.gz";得每个cutoff下被每种factor hit住hotspot的文件"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/20_count_hotspot_in_factor_pos.txt.gz"

perl 21_NHP_background_intersect_hotspot.pl #将"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/${tissue}_cis_eQTL_1kg_Completion.txt.gz” 转换为bed得 "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/${tissue}_cis_eQTL_1kg_Completion.bed.gz,排序得$sort_fo1，用$sort_fo1 和"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/${tissue}_cis_eQTL_1kg_Completion.bed.gz bedtools intersect得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/hotspot_cis_eQTL_overlap_SNP/${tissue}_segment_hotspot_cutoff_${cutoff}_SNP.bed.gz

perl 22_count_eqtl_in_hotspot_ratio.pl ##统计/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/hotspot_cis_eQTL_overlap_SNP/${tissue}_segment_hotspot_cutoff_${cutoff}_SNP.bed.gz中eQTL在SNP中的比例得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/${tissue}_cis_eQTL_in_hotspot_ratio.gz

Rscript 23_plot_boxplot.R
#--------------------
perl 30_annotation_marker.pl # #interval_18 时，对../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/whole_blood_segment_${group}_cutoff_${cutoff}.bed.gz用annotation_chromatin_states_interval18.sh进行annotation,得$output_dir/$factor_$input_file_base_name  #包括histone marker #my @markers = ("H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9ac","H3K9me3","CHROMATIN_Accessibility","TFBS","CTCF");
perl 31_filter_annotation_state.pl ##因为每个片段会被多个chrom state 注释，根据对hotspot的覆盖比例，选出覆盖比例最高的state, interval_18 时，对"../../output/${tissue}/Cis_eQTL/annotation/interval_18/ALL/${group}/${cutoff}/${state}_state_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz"进行过滤，得$dir/filter_$name
perl 31_1_count_chrom_state.pl ###对"$dir/filter_$name"中的chromatin state 进行统计，得"../../../output/${tissue}/Cis_eQTL/enrichment/interval_18/ALL/${type}_${state}_state_count_${tissue}_cutoff_${cutoff}.txt"

perl 32_count_anno_histone_mark.pl #对../../output/${tissue}/Cis_eQTL/annotation/interval_18/ALL/${group}/${cutoff}*的marker进行count,得$out_dir/${group}_histone_marker.txt.gz




#-----------
perl 33_calculate_jaccard_index_mark.pl #对 $input_dir/${mark}_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz 计算 jaccard index得$out_dir/${group}_cutoff_${cutoff}_histone_marker_jaccard_index.txt.gz



Rscript 34_qqman.R 

