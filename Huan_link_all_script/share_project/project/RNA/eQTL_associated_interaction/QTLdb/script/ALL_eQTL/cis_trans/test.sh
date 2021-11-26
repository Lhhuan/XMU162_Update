fraction=0.1
QTL="eQTL"
group="hotspot"
interval="18"
type="cis_1MB"
input_file="../../../output/ALL_${QTL}/cis_trans/${group}/interval_${interval}_cutoff_7.3_${type}_${QTL}_segment_${group}.bed.gz"
input_file_base_name="interval_${interval}_cutoff_7.3_${type}_${QTL}_segment_${group}.bed.gz"
output_dir="test_output"
bedtools intersect -F $fraction -a $input_file -b ../../../annotation_data/cistromeDB/HISTONE_MARK_AND_VARIANT/03_merge_pos_info_narrow_peak_sorted_chr1_22_extend_100bp.bed.gz -wa -wb | gzip > $output_dir/HISTONE_$input_file_base_name



bedtools intersect -F $fraction -a $input_file -b ../../../annotation_data/RBP/03_all_unique_merge_RBP_narrow_peak_sorted_chr1_22_extend_100bp.bed.gz -wa -wb | gzip > $output_dir/RBP_$input_file_base_name