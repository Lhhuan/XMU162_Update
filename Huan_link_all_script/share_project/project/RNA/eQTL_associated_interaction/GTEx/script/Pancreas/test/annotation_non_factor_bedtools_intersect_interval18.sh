bedtools intersect  -F 0.1 -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/ROADMAP/enhancer/01_enhancer_complement.bed.gz" -wa -wb | gzip > $output_dir/non_enhancer_$input_file_base_name
bedtools intersect -F 0.1 -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/ROADMAP/promoter/01_promoter_complement.bed.gz" -wa -wb | gzip > $output_dir/non_promoter_$input_file_base_name
