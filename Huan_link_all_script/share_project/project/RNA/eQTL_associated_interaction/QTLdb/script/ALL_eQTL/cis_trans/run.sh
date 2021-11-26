perl 01_Completion_snp_for_cis_trans_eQTL_by_1kg.pl # # 用../../../output/all_1kg_phase3_v5_hg19_snp.txt.gz 补全../../../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz 得/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/cis_1MB_eQTL_pop_1kg_Completion.txt.gz
#trans_1MB_eQTL_pop_1kg_Completion.txt.gz, trans_1MB_eQTL_pop_1kg_Completion.txt.gz, cis_10MB_eQTL_pop_1kg_Completion.txt.gz,
# Rscript 02_huan_NHPoisson.R
# Rscript 02_huan_NHPoisson_single.R
Rscript 02_huan_NHPoisson_par.R 
# #----------------------------------
#     perl 03_filter_emplambda_in_eQTL.pl #将../../../output/ALL_eQTL/cis_trans/NHPoisson_emplambda_interval_1000_cutoff_7.3_${name}_eQTL.txt.gz 中的QTLbase中../../../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz的eQTL数据过滤出来，
#     #得../../../output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_${name}_eQTL.txt.gz
#     #并得汇总文件"../../../output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL.txt.gz"
#     perl 04_normalize.pl #将"../../../output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL.txt.gz" 处理成$cis_or_trans:$emplambda的格式，并将位置处理为bed
#     #得../../../output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL_normalized.bed.gz
#     zless ../../../output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL_normalized.bed.gz |awk 'NR>1'|sort -k1,1 -k2,2n > ../../../output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL_normalized_sorted.bed
#     gzip ../../../output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL_normalized_sorted.bed 



#     # perl 04_merge_all_cis_trans.pl #将 ../../../output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_${name}_eQTL.txt.gz "cis_1MB","trans_1MB","cis_10MB","trans_10MB" merge到一起得
#     # #../../../output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_cis_trans.txt.gz
#     perl 04_cis_trans_normalized_bedgraph.pl

# #----------------------------------
Rscript 021_histgram_density_interval.R


perl 03_filter_hotspot_cis_trans_eQTL.pl # ###@QTLs =("eQTL")在@interval = (6,7,8,9,12,15)时的hotspot(segment),得"../../../output/ALL_${QTL}/cis_trans/hotspot/interval_${j}_cutoff_7.3_${type}_${QTL}_segment_hotspot.txt.gz"; 
#得point 热点"../../../output/ALL_${QTL}/cis_trans/hotspot/interval_${j}_cutoff_7.3_${type}_${QTL}_point_hotspot.txt.gz"
perl 03_filter_non_hotspot_cis_trans_eQTL.pl ####QTLs =eQTL,在@interval = (6,7,8,9,12,15)时的non hotspot(segment),得"../../../output/ALL_${QTL}/cis_trans/non_hotspot/interval_${j}_cutoff_7.3_${type}_${QTL}_segment_non_hotspot.txt.gz";
        # 得point 非热点"../../../output/ALL_${QTL}/cis_trans/non_hotspot/interval_${j}_cutoff_7.3_${type}_${QTL}_point_non_hotspot.txt.gz"
# perl 03_filter_hotspot_cis_trans_eQTL_interval_18.pl # ###@QTLs =("eQTL")在@interval = 18时的hotspot(segment),得"../../../output/ALL_${QTL}/cis_trans/hotspot/interval_${j}_cutoff_7.3_${type}_${QTL}_segment_hotspot.txt.gz"; 
# #得point 热点"../../../output/ALL_${QTL}/cis_trans/hotspot/interval_${j}_cutoff_7.3_${type}_${QTL}_point_hotspot.txt.gz"
# perl 03_filter_non_hotspot_cis_trans_eQTL_interval_18.pl ####QTLs =eQTL,在@interval = 18时的non hotspot(segment),得"../../../output/ALL_${QTL}/cis_trans/non_hotspot/interval_${j}_cutoff_7.3_${type}_${QTL}_segment_non_hotspot.txt.gz"
#及seegent_non_hotspot.bed.gz,得point 非热点"../../../output/ALL_${QTL}/cis_trans/non_hotspot/interval_${j}_cutoff_7.3_${type}_${QTL}_point_non_hotspot.txt.gz"
Rscript 031_histgram_density_length_of_segment.R #plot distribution of segment in hotspot and non-hotspot 
Rscript 031_histgram_density_length_of_segment_part.R 
perl 041_filter_hotspot_segment.pl #将../../../output/ALL_${QTL}/cis_trans/hotspot/interval_${j}_cutoff_7.3_${type}_${QTL}_segment_hotspot.bed.gz 长度>=6的片段，得
#../../../output/ALL_${QTL}/cis_trans/hotspot/interval_${j}_cutoff_7.3_${type}_${QTL}_segment_hotspot_length_more_than6.bed.gz
# Rscript 032_violin_boxplot_length_of_segment.R
#-------------------------------- interval 15
#         perl 04_annotation_interval_15.pl ##@types=("cis_1MB","cis_10MB","trans_1MB","trans_10MB")和@groups = ("hotspot","non_hotspot")时,用annotation_bedtools_intersect.sh进行annotation,得$output_dir/RBP_$input_file_base_name
#         Rscript 05_merge_trans_cis_eQTL_hotspot_interval15_annotation.R #将/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/$group/$type/ 下的文件进行annotation,得
#         #"~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/05_annotation_merge_${group}_${type}.txt.gz"
#         Rscript 05_merge_trans_cis_eQTL_hotspot_interval15_annotation_par.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_interval15_annotation_1.R 
#         Rscript 05_merge_trans_cis_eQTL_hotspot_interval15_annotation_2.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_interval15_annotation_3.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_interval15_annotation_4.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_interval15_annotation_5.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_interval15_annotation_6.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_interval15_annotation_7.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_interval15_annotation_8.R


#         #--------------------
#         perl 06_interval_15_count_factor_in_and_out_hotspot.pl #看各种factor在#"~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/01_annotation_merge_${group}_${type}.txt.gz"中hotspot 和 non-hotspotfactor的个数，得
#         #"~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/factors_score/06_annotation_merge_factors_score_${group}_${type}.txt.gz"
#         perl 07_interval_15_prepare_for_fisher.pl ##看各种factor在#"~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/01_annotation_merge_${group}_${type}.txt.gz 的fisher exact test 准备数据,
#         #得"/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/fisher_exact_test/prepare/07_interval_15_prepare_fisher_test.txt"
#         Rscript 08_interval_15_fisher_exact_test.R #对"/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/fisher_exact_test/prepare/07_interval_15_prepare_fisher_test.txt"
#         #进行fisher exact test,得"/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/fisher_exact_test/prepare/07_interval_15_prepare_fisher_test.txt"
#         Rscript 09_interval_15_boxplot_for_factor_score.R #对/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/factors_score/ 中同一种type进行boxplot比较
# #------------------------------@interval = (6,7,8,9,12,15,18)
#         perl 041_overlap_hotspot.pl ##@types=("cis_1MB","cis_10MB","trans_1MB","trans_10MB")和@groups = ("hotspot","non_hotspot")时,，对interval为15_6_7_8_9_12_18_时，用041_bedtools_intersect_hotspot_overlap.sh进行overlap进行overlap,
#         #得../../../output/ALL_eQTL/cis_trans/$group/segment_overlap/15_6_7_8_9_12_18_cutoff_7.3_${type}_eQTL_segment_$group.bed.gz
#         perl 042_annotation_overlap_interval.pl #@types=("cis_1MB","cis_10MB","trans_1MB","trans_10MB")和@groups = ("hotspot","non_hotspot")时, "../../../output/ALL_${QTL}/cis_trans/${group}/segment_overlap/${interval}_cutoff_7.3_${type}_${QTL}_segment_${group}.bed.gz"
#         #用annotation_bedtools_intersect.sh进行annotation,得$output_dir/RBP_$input_file_base_name
#         #-----------------------------------------------------------------------
#                 ##将/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/overlap/$group/$type/ 下的文件进行annotation,得
#                 #"~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/overlap/annotation_out/05_annotation_merge_${group}_${type}.txt.gz"
#         Rscript 05_merge_trans_cis_eQTL_hotspot_overlap_annotation_1.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_overlap_annotation_2.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_overlap_annotation_3.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_overlap_annotation_4.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_overlap_annotation_5.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_overlap_annotation_6.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_overlap_annotation_7.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_overlap_annotation_8.R
#         #----------------------------------------------------------------
#         perl 06_overlap_count_factor_in_and_out_hotspot.pl #看各种factor在#"~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/overlap/annotation_out/05_annotation_merge_${group}_${type}.txt.gz"中hotspot 和 non-hotspotfactor的个数，得
#         #"~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/overlap/annotation_out/factors_score/06_annotation_merge_factors_score_${group}_${type}.txt.gz"
#         perl 07_overlap_prepare_for_fisher.pl ##看各种factor在#"~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/overlap/annotation_out/05_annotation_merge_${group}_${type}.txt.gz 的fisher exact test 准备数据,
#         #得"/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/overlap/annotation_out/fisher_exact_test/prepare/07_interval_15_6_7_8_9_12_18_prepare_fisher_test.txt"
#         Rscript 08_overlap_fisher_exact_test.R #对"/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/overlap/annotation_out/fisher_exact_test/07_interval_15_6_7_8_9_12_18_prepare_fisher_test.txt"
#         #进行fisher exact test,得/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/overlap/annotation_out/fisher_exact_test/08_fisher_exact_test_result_interval_overlap.txt
#         Rscript 09_overlap_boxplot_for_factor_score.R #对/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/factors_score/ 中同一种type进行boxplot比较
#         #-------------------------------------------------
# #------------------------------interval 18
#         perl 04_annotation_interval_18.pl ##@types=("cis_1MB","cis_10MB","trans_1MB","trans_10MB")和@groups = ("hotspot","non_hotspot")时,用annotation_bedtools_intersect.sh进行annotation,得$output_dir/RBP_$input_file_base_name
#         # Rscript 05_merge_trans_cis_eQTL_hotspot_interval15_annotation.R #将/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/$group/$type/ 下的文件进行annotation,得
#         #"~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/05_annotation_merge_${group}_${type}.txt.gz"
#         Rscript 05_merge_trans_cis_eQTL_hotspot_interval18_annotation_1.R 
#         Rscript 05_merge_trans_cis_eQTL_hotspot_interval18_annotation_2.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_interval18_annotation_3.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_interval18_annotation_4.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_interval18_annotation_5.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_interval18_annotation_6.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_interval18_annotation_7.R
#         Rscript 05_merge_trans_cis_eQTL_hotspot_interval18_annotation_8.R


#         #--------------------
#         perl 06_interval_18_count_factor_in_and_out_hotspot.pl #看各种factor在#"~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_18/annotation_out/01_annotation_merge_${group}_${type}.txt.gz"中hotspot 和 non-hotspotfactor的个数，得
#         #"~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_18/annotation_out/factors_score/06_annotation_merge_factors_score_${group}_${type}.txt.gz"
#         perl 07_interval_18_prepare_for_fisher.pl ##看各种factor在#"~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_18/annotation_out/01_annotation_merge_${group}_${type}.txt.gz 的fisher exact test 准备数据,
#         #得"/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_18/annotation_out/fisher_exact_test/prepare/07_interval_18_prepare_fisher_test.txt"
#         Rscript 08_interval_18_fisher_exact_test.R #对"/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_18/annotation_out/fisher_exact_test/prepare/07_interval_18_prepare_fisher_test.txt"
#         #进行fisher exact test,得"/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_18/annotation_out/fisher_exact_test/prepare/07_interval_18_prepare_fisher_test.txt"
#         Rscript 09_interval_18_boxplot_for_factor_score.R #对/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_18/annotation_out/factors_score/ 中同一种type进行boxplot比较       
# #--------------------
#         perl 10_bedtools_fisher_test.pl 

#         bedtools fisher -F 0.05 -a ../../../output/ALL_eQTL/cis_trans/hotspot/interval_12_cutoff_7.3_cis_10MB_eQTL_segment_hotspot_length_of_segment.bed.gz  -b ../../../annotation_data/circRNA/03_needed_hsa_hg19_circRNA_sorted_chr1_22.bed.gz -g  /home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt




#         bedtools fisher -F 0.9 -a ../../../output/ALL_eQTL/cis_trans/hotspot/interval_12_cutoff_7.3_cis_10MB_eQTL_segment_hotspot_length_of_segment.bed.gz  -b ../../../annotation_data/enhancer/fantom5_enhancers_phase1_phase2_sorted.bed.gz -g  /home/huanhuan/ref_data/UCSC/hg19.chrom1_22_sizes_sorted.txt 

#------------------------------
        perl 04_annotation_different_interval.pl ##@types=("cis_1MB","cis_10MB","trans_1MB","trans_10MB")和@groups = ("hotspot","non_hotspot"),@intervals = (6,7,8,9,12,15,18),@fractions = (0.5,0.6,0.7,0.8,0.9,1)时,用annotation_bedtools_intersect.sh进行annotation,得$output_dir/RBP_$input_file_base_name
        perl 04_annotation_different_interval_18.pl
        Rscript 05_merge_trans_cis_eQTL_hotspot_annotation_interval6.R
        Rscript 05_merge_trans_cis_eQTL_hotspot_annotation_interval7.R
        Rscript 05_merge_trans_cis_eQTL_hotspot_annotation_interval8.R
        Rscript 05_merge_trans_cis_eQTL_hotspot_annotation_interval9.R
        Rscript 05_merge_trans_cis_eQTL_hotspot_annotation_interval12.R
        Rscript 05_merge_trans_cis_eQTL_hotspot_annotation_interval15.R
        Rscript 05_merge_trans_cis_eQTL_hotspot_annotation_interval15_test.R
        Rscript 05_merge_trans_cis_eQTL_hotspot_annotation_interval18.R
        Rscript 05_merge_trans_cis_eQTL_hotspot_annotation_interval18_test.R
        perl 07_prepare_for_fisher.pl ##看各种factor在#"~/project/RNA/eQTL_associated_interaction/QTLdb/output/ALL_eQTL/cis_trans/fisher_exact_test/annotation_out/hotspot/05_annotation_merge_${type}_interval_${interval}_overlap_fraction_${fraction}.txt.gz" 的fisher exact test 准备数据,
#得"/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/output/ALL_eQTL/cis_trans/fisher_exact_test/fisher_exact_test/interval_${interval}/07_prepare_fisher_test.txt"
        Rscript 08_fisher_exact_test.R
        Rscript 09_interval_18_factor_boxplot.R
        #------------------------
        perl 041_annotation_extend_different_interval.pl ##@types=("cis_1MB","cis_10MB","trans_1MB","trans_10MB")和@groups = ("hotspot","non_hotspot"),@intervals = (6,7,8,9,12,15,18),@fractions = (0.5,0.6,0.7,0.8,0.9,1)时,用annotation_bedtools_intersect_v2.sh进行annotation,得extend100bp时，$output_dir/RBP_$input_file_base_name
        Rscript 05_merge_trans_cis_eQTL_hotspot_annotation_extend_100bp_interval18.R
        Rscript 05_merge_trans_cis_eQTL_hotspot_annotation_extend_100bp_interval15.R
        perl 07_prepare_for_fisher_extend_100.pl
        Rscript 08_fisher_exact_test_extend_100.R
#----------------------------------------------
        perl 041_annotation_extend_500bp_different_interval.pl ##@types=("cis_1MB","cis_10MB","trans_1MB","trans_10MB")和@groups = ("hotspot","non_hotspot"),@intervals = (6,7,8,9,12,15,18),@fractions = (0.5,0.6,0.7,0.8,0.9,1)时,用annotation_bedtools_intersect_v2.sh进行annotation,得extend500bp时，$output_dir/RBP_$input_file_base_name
