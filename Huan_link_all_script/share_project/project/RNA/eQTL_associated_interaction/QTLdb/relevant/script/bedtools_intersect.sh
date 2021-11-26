bedtools intersect -a ../../output/ALL_eQTL/hotspot/segment_overlap/15_6_7_8_9_12_18_segment_hotspot.bed.gz -b ../../annotation_data/enhancer/fantom5_enhancers_phase1_phase2_sorted.bed.gz -wa -wb | gzip > ../output/annotation/QTLbase_all_eQTL_15_6_7_8_9_12_18_enhancer.bed.gz
bedtools intersect -a ../../output/ALL_eQTL/hotspot/segment_overlap/15_6_7_8_9_12_18_segment_hotspot.bed.gz -b ../../annotation_data/promoter/fantom5_promoter_phase1_phase2_sorted.bed.gz -wa -wb | gzip > ../output/QTLbase_all_eQTL_15_6_7_8_9_12_18_promoter.bed.gz

bedtools intersect -a ../../output/ALL_eQTL/hotspot/segment_overlap/15_6_7_8_9_12_18_segment_hotspot.bed.gz -b ../../annotation_data/circRNA/needed_hsa_hg19_circRNA_sorted.bed.gz -wa -wb | gzip > ../output/annotation/QTLbase_all_eQTL_15_6_7_8_9_12_18_circRNA.bed.gz
bedtools intersect -a ../../output/ALL_eQTL/hotspot/segment_overlap/15_6_7_8_9_12_18_segment_hotspot.bed.gz -b ../../annotation_data/hic/02_all_unique_hic_loops_tran_cis_1_sorted.bed.gz -wa -wb | gzip > ../output/annotation/QTLbase_all_eQTL_15_6_7_8_9_12_18_hic_loops_1.bed.gz
bedtools intersect -a ../../output/ALL_eQTL/hotspot/segment_overlap/15_6_7_8_9_12_18_segment_hotspot.bed.gz -b ../../annotation_data/hic/02_all_unique_hic_loops_tran_cis_2_sorted.bed.gz -wa -wb | gzip > ../output/annotation/QTLbase_all_eQTL_15_6_7_8_9_12_18_hic_loops_2.bed.gz
bedtools intersect -a ../../output/ALL_eQTL/hotspot/segment_overlap/15_6_7_8_9_12_18_segment_hotspot.bed.gz -b ../../annotation_data/lncRNA/NONCODEv5_hg19.lncAndGene_sorted.bed.gz -wa -wb | gzip > ../output/annotation/QTLbase_all_eQTL_15_6_7_8_9_12_18_lncRNA.bed.gz
bedtools intersect -a ../../output/ALL_eQTL/hotspot/segment_overlap/15_6_7_8_9_12_18_segment_hotspot.bed.gz -b ../../annotation_data/CTCF/02_all_unique_merge_ctcf_binding_site_narrow_peak_sorted.bed.gz -wa -wb | gzip >  ../output/annotation/QTLbase_all_eQTL_15_6_7_8_9_12_18_CTCF.bed.gz

bedtools intersect -a ../../output/ALL_eQTL/hotspot/segment_overlap/15_6_7_8_9_12_18_segment_hotspot.bed.gz -b ../../annotation_data/cistromeDB/Human_FACTOR/merge_pos_info_narrow_peak_sorted.bed.gz -wa -wb | gzip > ../output/annotation/QTLbase_all_eQTL_15_6_7_8_9_12_18_cistrome_TFBS_narrow_peak.bed.gz
bedtools intersect -a ../../output/ALL_eQTL/hotspot/segment_overlap/15_6_7_8_9_12_18_segment_hotspot.bed.gz -b ../../annotation_data/cistromeDB/Human_CHROMATIN_Accessibility/merge_pos_info_narrow_peak_sorted.bed.gz -wa -wb | gzip > ../output/annotation/QTLbase_all_eQTL_15_6_7_8_9_12_18_cistrome_CHROMATIN_Accessibility_narrow_peak.bed.gz
bedtools intersect -a ../../output/ALL_eQTL/hotspot/segment_overlap/15_6_7_8_9_12_18_segment_hotspot.bed.gz -b ../../annotation_data/cistromeDB/HISTONE_MARK_AND_VARIANT/merge_pos_info_narrow_peak_sorted.bed.gz -wa -wb | gzip > ../output/annotation/QTLbase_all_eQTL_15_6_7_8_9_12_18_cistrome_HISTONE_narrow_peak.bed.gz
bedtools intersect -a ../../output/ALL_eQTL/hotspot/segment_overlap/15_6_7_8_9_12_18_segment_hotspot.bed.gz -b ../../annotation_data/RBP/02_all_unique_merge_RBP_narrow_peak_sorted.bed.gz -wa -wb | gzip > ../output/annotation/QTLbase_all_eQTL_15_6_7_8_9_12_18_RBP_narrow_peak.bed.gz

bedtools intersect -a ../../output/ALL_eQTL/hotspot/segment_overlap/15_6_7_8_9_12_18_segment_hotspot.bed.gz -b ../../annotation_data/conservation/phastCons100way/merge_phastCons100way_sorted.bed.gz -wa -wb | gzip > ../output/annotation/QTLbase_all_eQTL_15_6_7_8_9_12_18_phastCons100way.bed.gz

# for i in $( seq 22);
# do
# 	echo $i;
#     bedtools intersect -a ../../output/ALL_eQTL/hotspot/segment_overlap/15_6_7_8_9_12_18_segment_hotspot.bed.gz -b ../../annotation_data/conservation/phastCons100way/normalized_per_chr/chr$i.sorted_merge_phastCons100way.bed.gz -wa -wb | gzip > ../output/conservation_per_chr/QTLbase_chr$i.all_eQTL_phastCons100way.bed.gz
# done
# zcat ../output/conservation_per_chr/* | gzip > ../output/annotation/QTLbase_all_eQTL_15_6_7_8_9_12_18_phastCons100way.bed.gz