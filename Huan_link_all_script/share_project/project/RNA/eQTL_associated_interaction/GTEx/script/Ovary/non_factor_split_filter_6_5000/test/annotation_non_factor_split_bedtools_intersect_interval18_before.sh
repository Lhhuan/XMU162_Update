bedtools intersect   -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/ROADMAP/Lung/enhancer/split_enhancer_complement.bed.gz" -wa -wb | gzip > $output_dir/non_enhancer_$input_file_base_name
bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/ROADMAP/Lung/promoter/split_promoter_complement.bed.gz" -wa -wb | gzip > $output_dir/non_promoter_$input_file_base_name

bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/Human_FACTOR/merge_pos_info_narrow_peak_complement_split.bed.gz" -wa -wb | gzip > $output_dir/non_TFBS_$input_file_base_name

bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/Human_CHROMATIN_Accessibility/merge_pos_info_narrow_peak_complement_split.bed.gz" -wa -wb | gzip > $output_dir/non_CHROMATIN_Accessibility_$input_file_base_name

bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/HISTONE_MARK_AND_VARIANT/merge_pos_info_narrow_peak_complement_split.bed.gz" -wa -wb | gzip > $output_dir/non_HISTONE_modification_$input_file_base_name


bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/CTCF/normal_cell_line/split_CTCF_complement.bed.gz" -wa -wb | gzip > $output_dir/non_CTCF_$input_file_base_name