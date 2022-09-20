bedtools intersect -a $sorted_input_file -b "$anno_dir/H3K27ac_merge_max_signalvalue.bedgraph.gz" -wo | gzip > $output_dir/H3K27ac_$input_file_base_name
bedtools intersect -a $sorted_input_file -b "$anno_dir/H3K27me3_merge_max_signalvalue.bedgraph.gz" -wo | gzip > $output_dir/H3K27me3_$input_file_base_name

bedtools intersect -a $sorted_input_file -b "$anno_dir/H3K36me3_merge_max_signalvalue.bedgraph.gz" -wo | gzip > $output_dir/H3K36me3_$input_file_base_name
bedtools intersect -a $sorted_input_file -b "$anno_dir/H3K4me1_merge_max_signalvalue.bedgraph.gz" -wo | gzip > $output_dir/H3K4me1_$input_file_base_name
bedtools intersect -a $sorted_input_file -b "$anno_dir/H3K4me3_merge_max_signalvalue.bedgraph.gz" -wo | gzip > $output_dir/H3K4me3_$input_file_base_name
bedtools intersect -a $sorted_input_file -b "$anno_dir/H3K9ac_merge_max_signalvalue.bedgraph.gz" -wo | gzip > $output_dir/H3K9ac_$input_file_base_name
bedtools intersect -a $sorted_input_file -b "$anno_dir/H3K9me3_merge_max_signalvalue.bedgraph.gz" -wo | gzip > $output_dir/H3K9me3_$input_file_base_name


#hg38
bedtools intersect -a $sorted_input_file -b "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/Human_CHROMATIN_Accessibility/merge_max_signalvalue.bedgraph.gz" -wo | gzip > $output_dir/CHROMATIN_Accessibility_$input_file_base_name
bedtools intersect -a $sorted_input_file -b "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/Human_FACTOR/merge_max_signalvalue.bedgraph.gz" -wo | gzip > $output_dir/TFBS_$input_file_base_name


#hg38
bedtools intersect -a $sorted_input_file -b "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/CTCF/normal_cell_line/hg38/normal_cell_line_ctcf_merge_max_signalvalue.bedgraph.gz" -wo | gzip > $output_dir/CTCF_$input_file_base_name

