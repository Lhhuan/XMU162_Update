bedtools intersect   -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/${tissue}/enhancer_union.bed.gz" -wa -wb | gzip > $output_dir/enhancer_$input_file_base_name
bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/${tissue}/promoter_union.bed.gz" -wa -wb | gzip > $output_dir/promoter_$input_file_base_name
bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/${tissue}/TFBS_union.bed.gz" -wa -wb | gzip > $output_dir/TFBS_$input_file_base_name

bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/${tissue}/CHROMATIN_Accessibility_union.bed.gz" -wa -wb | gzip > $output_dir/CHROMATIN_Accessibility_$input_file_base_name

bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/${tissue}/HISTONE_modification_union.bed.gz" -wa -wb | gzip > $output_dir/HISTONE_modification_$input_file_base_name

bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/${tissue}/CTCF_union.bed.gz" -wa -wb | gzip > $output_dir/CTCF_$input_file_base_name