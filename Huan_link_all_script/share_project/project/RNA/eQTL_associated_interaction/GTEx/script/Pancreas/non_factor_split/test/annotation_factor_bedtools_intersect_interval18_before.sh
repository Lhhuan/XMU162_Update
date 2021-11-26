bedtools intersect   -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/ROADMAP/Lung/enhancer/01_normal_E096-H3K4me1_H3K27ac.narrowPeaksorted_overlap.gz" -wa -wb | gzip > $output_dir/enhancer_$input_file_base_name
bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/ROADMAP/Lung/promoter/01_normal_E096-H3K4me3.narrowPeaksorted.gz" -wa -wb | gzip > $output_dir/promoter_$input_file_base_name
bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/Human_FACTOR/merge_pos_info_narrow_peak_sort.bed.gz" -wa -wb | gzip > $output_dir/TFBS_$input_file_base_name

bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/Human_CHROMATIN_Accessibility/merge_pos_info_narrow_peak_sort.bed.gz" -wa -wb | gzip > $output_dir/CHROMATIN_Accessibility_$input_file_base_name

bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/HISTONE_MARK_AND_VARIANT/merge_pos_info_narrow_peak_sort.bed.gz" -wa -wb | gzip > $output_dir/HISTONE_modification_$input_file_base_name

bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/CTCF/normal_cell_line/05_normal_cell_line_ctcf_sort.bed.gz" -wa -wb | gzip > $output_dir/CTCF_$input_file_base_name