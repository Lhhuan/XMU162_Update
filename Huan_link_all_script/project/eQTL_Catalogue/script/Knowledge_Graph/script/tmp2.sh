# perl 03_nodes_annotation_markers.pl
bedtools intersect -a "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_5e_8.bed.gz" -b "/home/huanhuan/project/link_database/JEME/output/hg38/01_merge_enhancer_target_sample_sorted.bed.gz" -wa -wb |gzip >../output/nodes_annotation/JEME_enhancer.bed.gz
perl 04_adjust_nodes_annotation_markers_format.pl
echo "04"
perl 0401_adjust_enhancer_format.pl
echo "0401"
bedtools intersect -a "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" -b "/home/huanhuan/project/eQTL_Catalogue/output/01_merge_all_tissue_cis_eQTLtissue_label_sorted.bed.gz" -wo |gzip > ../output/nodes_annotation/hotspot_tissue_anno.bed.gz
perl 041_adjust_hotspot_tissue_anno.pl 
echo "041"
Rscript 05_merge_node_markers_annotation.R
echo "05"
gzip ../output/nodes_annotation/nodes_marker_enhancers_eqtl_annotation.txt
zless ../output/nodes_annotation/nodes_marker_enhancers_eqtl_annotation.txt.gz |cut -f1-3,5-15|gzip >../output/nodes_annotation/hotspot_annotation.txt.gz
Rscript 051_adjust_hotspot_annotation.R
echo "051"
gzip ../output/nodes_annotation/hotspot_annotation.txt
perl 06_judge_enhancer_target_edge.pl
echo "06"
perl 07_annotation_hotspot_gene_tissue.pl 
echo "07"
perl 08_annotation_enhancer_promoter_interaction.pl
echo "08"
perl 09_merge_egene_and_pos.pl
echo "09"
perl 091_adjust_TAD_format.pl
echo "091"
Rscript 10_nodes_and_edges_annotation.R
echo "10"
Rscript 11_plot_nodes_edges_features_hotspot_as_background.R
echo "11"
perl 12_statistic_ratio_enhancer_target_as_background.pl 
echo "12"
