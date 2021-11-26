bedtools intersect   -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/ROADMAP/Lung/enhancer/split_enhancer_complement.bed.gz" -wa -wb | gzip > $output_dir/non_enhancer_$input_file_base_name
bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/ROADMAP/Lung/promoter/split_promoter_complement.bed.gz" -wa -wb | gzip > $output_dir/non_promoter_$input_file_base_name

bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/Lung/non_TFBS_split_union.bed.gz" -wa -wb | gzip > $output_dir/non_TFBS_$input_file_base_name

bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/Lung/non_CHROMATIN_Accessibility_split_union.bed.gz" -wa -wb | gzip > $output_dir/non_CHROMATIN_Accessibility_$input_file_base_name

bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/Lung/non_HISTONE_modification_split_union.bed.gz" -wa -wb | gzip > $output_dir/non_HISTONE_modification_$input_file_base_name


bedtools intersect -a $input_file -b "/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/Lung/non_CTCF_split_union.bed.gz" -wa -wb | gzip > $output_dir/non_CTCF_$input_file_base_name