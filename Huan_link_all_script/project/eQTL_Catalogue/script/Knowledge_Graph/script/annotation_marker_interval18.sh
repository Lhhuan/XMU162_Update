bedtools intersect -a $sorted_input_file -b "$anno_dir/H3K27ac_sorted_sample.bed.gz" -wo | gzip > $output_dir/H3K27ac_$input_file_base_name
bedtools intersect -a $sorted_input_file -b "$anno_dir/H3K27me3_sorted_sample.bed.gz" -wo | gzip > $output_dir/H3K27me3_$input_file_base_name

bedtools intersect -a $sorted_input_file -b "$anno_dir/H3K36me3_sorted_sample.bed.gz" -wo | gzip > $output_dir/H3K36me3_$input_file_base_name
bedtools intersect -a $sorted_input_file -b "$anno_dir/H3K4me1_sorted_sample.bed.gz" -wo | gzip > $output_dir/H3K4me1_$input_file_base_name
bedtools intersect -a $sorted_input_file -b "$anno_dir/H3K4me3_sorted_sample.bed.gz" -wo | gzip > $output_dir/H3K4me3_$input_file_base_name
bedtools intersect -a $sorted_input_file -b "$anno_dir/H3K9ac_sorted_sample.bed.gz" -wo | gzip > $output_dir/H3K9ac_$input_file_base_name
bedtools intersect -a $sorted_input_file -b "$anno_dir/H3K9me3_sorted_sample.bed.gz" -wo | gzip > $output_dir/H3K9me3_$input_file_base_name


#hg38
bedtools intersect -a $sorted_input_file -b "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/Human_CHROMATIN_Accessibility/merge_pos_info_sample_narrow_peak_sorted.bed.gz" -wo | gzip > $output_dir/CHROMATIN_Accessibility_$input_file_base_name
bedtools intersect -a $sorted_input_file -b "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/Human_FACTOR/merge_pos_info_sample_narrow_peak_sorted.bed.gz" -wo | gzip > $output_dir/TFBS_$input_file_base_name


#hg38
bedtools intersect -a $sorted_input_file -b "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/CTCF/normal_cell_line/hg38/05_normal_cell_line_ctcf_sample_sorted.bed.gz" -wo | gzip > $output_dir/CTCF_$input_file_base_name

bedtools intersect -a $sorted_input_file -b "/home/huanhuan/project/link_database/EnhancerAtlas/01_merge_enahcner_sample_sorted.bed.gz" -wo | gzip > $output_dir/Enhancer_$input_file_base_name