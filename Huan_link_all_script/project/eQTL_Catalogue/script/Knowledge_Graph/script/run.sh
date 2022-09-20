Rscript demo_plot.R 
perl 01_merge_target_gene_reactomeFI.pl #"../../../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz" 与reactomeFI "/home/huanhuan/project/link_database/reactome_FI/output/04_ReactomeFI_part1_start_end.txt.gz", "/home/huanhuan/project/link_database/reactome_FI/output/04_ReactomeFI_part2_end_start.txt.gz","/home/huanhuan/project/link_database/reactome_FI/output/04_ReactomeFI_part3_both.txt.gz", "/home/huanhuan/project/link_database/reactome_FI/output/04_ReactomeFI_part3_both.txt.gz" merge在一起得 ./output/01_hotspot_target_gene_reactomeFI.bed.gz 
perl 02_annotation_co-expression.pl #为./output/01_hotspot_target_gene_reactomeFI.bed.gz  annotation co-expression数据，得./ENSG_G16808_S85825/${file_entrezgene}_${file_ensembl}，得./output/01_hotspot_target_gene_reactomeFI_co-expression.bed.gz 
zless ../output/02_hotspot_target_gene_reactomeFI_co-expression.bed.gz  |cut -f1-4,7-8 |gzip > ../output/02_hotspot_target_gene_bidirection_reactomeFI_co-expression.bed.gz 
                                                #-------------------------------------nodes annotation-------------------------------------

perl 03_nodes_annotation_markers.pl ##用带有样本信息的10种marker及enhancer 对"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" 进行注释得../output/nodes_annotation/${markers}_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz

bedtools intersect -a "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_5e_8.bed.gz" -b "/home/huanhuan/project/link_database/JEME/output/hg38/01_merge_enhancer_target_sample_sorted.bed.gz" -wa -wb |gzip >../output/nodes_annotation/JEME_enhancer.bed.gz
perl 04_adjust_nodes_annotation_markers_format.pl #对../output/nodes_annotation/${markers}_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz 格式进行调整，得../output/nodes_annotation/Adjust_${markers}_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz
perl 0401_adjust_enhancer_format.pl ##对合并enhancer annotation"../output/nodes_annotation/EnhancerAtlas_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz"，"../output/nodes_annotation/JEME_enhancer.bed.gz"信息得"../output/nodes_annotation/Adjust_${marker}_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz"
bedtools intersect -a "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" -b "/home/huanhuan/project/eQTL_Catalogue/output/01_merge_all_tissue_cis_eQTLtissue_label_sorted.bed.gz" -wo |gzip > ../output/nodes_annotation/hotspot_tissue_anno.bed.gz
perl 041_adjust_hotspot_tissue_anno.pl # 调整../output/nodes_annotation/hotspot_tissue_anno.bed.gz 的格式得../output/nodes_annotation/Adjust_hotspot_tissue_anno.bed.gz
Rscript 05_merge_node_markers_annotation.R  #../output/nodes_annotation/nodes_marker_enhancers_eqtl_annotation.txt
gzip ../output/nodes_annotation/nodes_marker_enhancers_eqtl_annotation.txt
# zless ../output/nodes_annotation/nodes_marker_enhancers_eqtl_annotation.txt.gz |cut -f1-3,5-15|gzip >../output/nodes_annotation/hotspot_annotation.txt.gz
Rscript 051_adjust_hotspot_annotation.R # ../output/nodes_annotation/hotspot_annotation.txt
gzip ../output/nodes_annotation/hotspot_annotation.txt

                                                #-------------------------------------edges annotation---------------------------------------

# bedtools intersect -a "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_5e_8.bed.gz" -b "/home/huanhuan/project/link_database/EnhancerAtlas/02_merge_enhancer_target_sample_sorted.bed.gz" -wo |gzip > ../output/edges_annotation/enhancer_target_anno.bed.gz
perl 06_judge_enhancer_target_edge.pl #筛选 ../output/nodes_annotation/JEME_enhancer.bed.gz中enhancer-target gene和 egene是否是相同基因，得../output/edges_annotation/success_enhancer_target_anno.bed.gz

perl 07_annotation_hotspot_gene_tissue.pl # #/home/huanhuan/project/eQTL_Catalogue/output/02_gene_chr_split/02_merge_all_tissue_cis_eQTL_eur_egene_sig_5e_8_sorted_chr${i}.bed.gz 对"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz"进行 annotation得../output/edges_annotation/07_annotation_hotspot_egene_tissue.bed.gz, 然后汇总调整格式得"../output/edges_annotation/07_annotation_hotspot_egene_tissue.bed.gz"

perl 08_annotation_enhancer_promoter_interaction.pl ## 用/home/huanhuan/project/link_database/OncoBase/output/hg38/01_${type}_sorted.bed.gz 对"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_5e_8.bed.gz"进行annotation,得"../output/edges_annotation/success_${type}_hotspot_egene.bed.gz"，adjust format得"../output/edges_annotation/${type}_hotspot_egene.bed.gz"
perl 09_merge_egene_and_pos.pl #"/home/huanhuan/reference/grch38_ensg_pos_from_ensembl106.txt.gz" 为"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_5e_8.bed.gz" 注释ensg位置得../output/09_hotspot_5e_8egene_pos.bed.gz
#与TAD /home/huanhuan/project/link_database/ENCOED/output/hg38/01_merge_TAD_sample_sorted.bed.gz intersect得"../output/edges_annotation/09_TAD_hotspot.bed.gz"， 将gene位置放前面得"../output/edges_annotation/09_Adjust_TAD_hotspot.bed.gz"，用gene 位置与TAD intersect得"../output/edges_annotation/09_TAD_egene.bed.gz"， 判断hotspot-egene在相同TAD得"../output/edges_annotation/09_success_egene_hotspot_TAD.bed.gz"

perl 091_adjust_TAD_format.pl #Adjust "../output/edges_annotation/09_success_egene_hotspot_TAD.bed.gz"得 "../output/edges_annotation/09_Adjust_success_egene_hotspot_TAD.bed.gz"
#----------------
Rscript 10_nodes_and_edges_annotation.R #
Rscript 11_plot_nodes_edges_features_hotspot_as_background.R #



#-------------------statistic by enhancer-target(background)
perl 12_statistic_ratio_enhancer_target_as_background.pl ## 统计("TSS_TSS","ENH_ENH","TSS_ENH"),enhancer-target 被hotspot cover 的比例，得../output/edges_annotation/12_cover_enhancer_target_type_ratio.txt.gz
Rscript 13_plot_edges_cover_enhancer_target_ratio.R #