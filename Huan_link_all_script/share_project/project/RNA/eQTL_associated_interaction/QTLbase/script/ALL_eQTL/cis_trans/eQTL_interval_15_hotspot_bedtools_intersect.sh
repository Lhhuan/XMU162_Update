bedtools intersect -a ../../output/ALL_eQTL/hotspot/interval_15_segment_hotspot.bed.gz -b ../../annotation_data/enhancer/fantom5_enhancers_phase1_phase2_sorted.bed.gz -wa -wb | gzip > ../output/annotation/hotspot/QTLbase_all_eQTL_interval_hotspot_15_enhancer.bed.gz
bedtools intersect -a ../../output/ALL_eQTL/hotspot/interval_15_segment_hotspot.bed.gz -b ../../annotation_data/promoter/fantom5_promoter_phase1_phase2_sorted.bed.gz -wa -wb | gzip > ../output/annotation/hotspot/QTLbase_all_eQTL_interval_hotspot_15_promoter.bed.gz

bedtools intersect -a ../../output/ALL_eQTL/hotspot/interval_15_segment_hotspot.bed.gz -b ../../annotation_data/circRNA/needed_hsa_hg19_circRNA_sorted.bed.gz -wa -wb | gzip > ../output/annotation/hotspot/QTLbase_all_eQTL_interval_hotspot_15_circRNA.bed.gz
# bedtools intersect -a ../../output/ALL_eQTL/hotspot/interval_15_segment_hotspot.bed.gz -b ../../annotation_data/hic/02_all_unique_hic_loops_tran_cis_1_sorted.bed.gz -wa -wb | gzip > ../output/annotation/hotspot/QTLbase_all_eQTL_interval_hotspot_15_hic_loops_1.bed.gz
bedtools intersect -a ../../output/ALL_eQTL/hotspot/interval_15_segment_hotspot.bed.gz -b ../../annotation_data/hic/02_all_unique_hic_loops_tran_cis_1_2_sorted.bed.gz -wa -wb | gzip > ../output/annotation/hotspot/QTLbase_all_eQTL_interval_hotspot_15_hic_loops_1_2.bed.gz
bedtools intersect -a ../../output/ALL_eQTL/hotspot/interval_15_segment_hotspot.bed.gz -b ../../annotation_data/lncRNA/NONCODEv5_hg19.lncAndGene_sorted.bed.gz -wa -wb | gzip > ../output/annotation/hotspot/QTLbase_all_eQTL_interval_hotspot_15_lncRNA.bed.gz
bedtools intersect -a ../../output/ALL_eQTL/hotspot/interval_15_segment_hotspot.bed.gz -b ../../annotation_data/CTCF/02_all_unique_merge_ctcf_binding_site_narrow_peak_sorted.bed.gz -wa -wb | gzip >  ../output/annotation/hotspot/QTLbase_all_eQTL_interval_hotspot_15_CTCF.bed.gz

bedtools intersect -a ../../output/ALL_eQTL/hotspot/interval_15_segment_hotspot.bed.gz -b ../../annotation_data/cistromeDB/Human_FACTOR/merge_pos_info_narrow_peak_sorted.bed.gz -wa -wb | gzip > ../output/annotation/hotspot/QTLbase_all_eQTL_interval_hotspot_15_cistrome_TFBS_narrow_peak.bed.gz
bedtools intersect -a ../../output/ALL_eQTL/hotspot/interval_15_segment_hotspot.bed.gz -b ../../annotation_data/cistromeDB/Human_CHROMATIN_Accessibility/merge_pos_info_narrow_peak_sorted.bed.gz -wa -wb | gzip > ../output/annotation/hotspot/QTLbase_all_eQTL_interval_hotspot_15_cistrome_CHROMATIN_Accessibility_narrow_peak.bed.gz
bedtools intersect -a ../../output/ALL_eQTL/hotspot/interval_15_segment_hotspot.bed.gz -b ../../annotation_data/cistromeDB/HISTONE_MARK_AND_VARIANT/merge_pos_info_narrow_peak_sorted.bed.gz -wa -wb | gzip > ../output/annotation/hotspot/QTLbase_all_eQTL_interval_hotspot_15_cistrome_HISTONE_narrow_peak.bed.gz
bedtools intersect -a ../../output/ALL_eQTL/hotspot/interval_15_segment_hotspot.bed.gz -b ../../annotation_data/RBP/02_all_unique_merge_RBP_narrow_peak_sorted.bed.gz -wa -wb | gzip > ../output/annotation/hotspot/QTLbase_all_eQTL_interval_hotspot_15_RBP_narrow_peak.bed.gz

# bedtools intersect -a ../../output/ALL_eQTL/hotspot/interval_15_segment_hotspot.bed.gz -b ../../annotation_data/conservation/phastCons100way/merge_phastCons100way_sorted.bed.gz -wa -wb | gzip > ../output/annotation/hotspot/QTLbase_all_eQTL_interval_hotspot_15_phastCons100way.bed.gz

# mkdir -p ../output/15_hot_conservation_per_chr/
for i in $( seq 22);
do
	echo $i;
    bedtools intersect -a ../../output/ALL_eQTL/hotspot/interval_15_segment_hotspot.bed.gz -b ../../annotation_data/conservation/phastCons100way/normalized_per_chr/chr$i.sorted_merge_phastCons100way.bed.gz -wa -wb | gzip > ../output/15_hot_conservation_per_chr/QTLbase_chr$i.all_eQTL_phastCons100way.bed.gz
done
zcat ../output/15_hot_conservation_per_chr/* | gzip > ../output/annotation/hotspot/QTLbase_all_eQTL_interval_hotspot_15_phastCons100way.bed.gz