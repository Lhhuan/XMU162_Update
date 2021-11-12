perl 07_filter_chr1_11.pl #过滤"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz 中的chr1_11,得"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/chr1_11/Chr1_11_Tissue_merge_0.176_extend_sorted_merge.bed.gz";

bash 6kmer.sh 


Rscript find_optimal_k.R 


perl 08_transform_kmer_result_other_n_community.pl
#go to /public/home/huanhuan/project/GTEx/homer/script run homer
# perl 11_split_cluster_and_homer_other_n_community_par.pl # "/public/home/huanhuan/project/GTEx/homer/script/11_split_cluster_and_homer_other_n_community_par.pl" 
perl 09_community_motif_overlap.pl #统计每个community中各个部分显著motif"${output_dir}/knowResult_merge.txt"， overlap情况得"${output_dir}/knowResult_merge_class_merge.txt",overlap ratio 文件"./6kmer/community_motif_overlap_ratio.txt"
perl 091_community_motif_overlap_community5_2_3.pl #
Rscript 09_3community_pie.R
Rscript 09_4community_pie.R
Rscript 09_5community_pie.R
Rscript 09_6community_pie.R
Rscript 09_7community_pie.R
Rscript 09_8community_pie.R
Rscript 09_9community_pie.R
Rscript 10_venn_community_tf.R 
#--------------------------------------
Rscript 10_venn_community3_tf.R
Rscript 10_venn_community4_tf.R
Rscript 10_venn_community5_tf.R
Rscript 10_venn_community6_tf.R
Rscript 10_venn_community_tf_5_union_23.R
#-------------------------------------
Rscript 11_diff_venn_community1_5_3.R

Rscript 12_plot_distribution_length_of_c1_and_other_com.R 

#-------------------testkkl_kmer.sh 
bash test_kkl.sh 
