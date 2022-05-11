Rscript demo_plot.R 
perl 01_merge_target_gene_reactomeFI.pl #"../../../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz" 与reactomeFI "/home/huanhuan/project/link_database/reactome_FI/output/04_ReactomeFI_part1_start_end.txt.gz", "/home/huanhuan/project/link_database/reactome_FI/output/04_ReactomeFI_part2_end_start.txt.gz","/home/huanhuan/project/link_database/reactome_FI/output/04_ReactomeFI_part3_both.txt.gz", "/home/huanhuan/project/link_database/reactome_FI/output/04_ReactomeFI_part3_both.txt.gz" merge在一起得 ./output/01_hotspot_target_gene_reactomeFI.bed.gz 
perl 02_annotation_co-expression.pl #为./output/01_hotspot_target_gene_reactomeFI.bed.gz  annotation co-expression数据，得./ENSG_G16808_S85825/${file_entrezgene}_${file_ensembl}，得./output/01_hotspot_target_gene_reactomeFI_co-expression.bed.gz 

                                                #-------------------------------------nodes annotation-------------------------------------

perl 03_nodes_annotation_markers.pl ##用带有样本信息的10种marker及enhancer 对"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" 进行注释得../output/nodes_annotation/${markers}_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz
perl 04_adjust_nodes_annotation_markers_format.pl #对../output/nodes_annotation/${markers}_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz 格式进行调整，得../output/nodes_annotation/Adjust_${markers}_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz
bedtools intersect -a "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" -b "/home/huanhuan/project/eQTL_Catalogue/output/01_merge_all_tissue_cis_eQTLtissue_label_sorted.bed.gz" -wo |gzip > ../output/nodes_annotation/hotspot_tissue_anno.bed.gz
perl 041_adjust_hotspot_tissue_anno.pl # 调整../output/nodes_annotation/hotspot_tissue_anno.bed.gz 的格式得../output/nodes_annotation/Adjust_hotspot_tissue_anno.bed.gz
Rscript 05_merge_node_markers_annotation.R  #../output/nodes_annotation/nodes_marker_enhancers_eqtl_annotation.txt
gzip ../output/nodes_annotation/nodes_marker_enhancers_eqtl_annotation.txt

                                                #-------------------------------------edges annotation---------------------------------------

bedtools intersect -a "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz" -b "/home/huanhuan/project/link_database/EnhancerAtlas/02_merge_enhancer_target_sample_sorted.bed.gz" -wo |gzip > ../output/edges_annotation/enhancer_target_anno.bed.gz
perl 06_judge_enhancer_target_edge.pl #筛选 ../output/edges_annotation/enhancer_target_anno.bed.gz中enhancer-target gene和 egene是否是相同基因，得../output/edges_annotation/success_enhancer_target_anno.bed.gz

bedtools intersect -a "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" -b "/home/huanhuan/project/eQTL_Catalogue/output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_0.05_sorted.bed.gz" -wa -wb |cut -f1-3,7-9|gzip >../output/edges_annotation/hotspot_egene_tissue_anno.bed.gz
bedtools intersect -a 1234.bed.gz -b 345.bed.gz -wo |gzip >re.bed.gz

#----------------