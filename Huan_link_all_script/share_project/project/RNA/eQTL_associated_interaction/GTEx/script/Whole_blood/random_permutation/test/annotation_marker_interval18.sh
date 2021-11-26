bedtools intersect -a $sorted_input_file -b "$anno_dir/$roadmap_biospecimen-H3K27ac.narrowPeak.gz" -wo | gzip > $output_dir/H3K27ac_$input_file_base_name
bedtools intersect -a $sorted_input_file -b "$anno_dir/$roadmap_biospecimen-H3K27me3.narrowPeak.gz" -wo | gzip > $output_dir/H3K27me3_$input_file_base_name

bedtools intersect -a $sorted_input_file -b "$anno_dir/$roadmap_biospecimen-H3K36me3.narrowPeak.gz" -wo | gzip > $output_dir/H3K36me3_$input_file_base_name
bedtools intersect -a $sorted_input_file -b "$anno_dir/$roadmap_biospecimen-H3K4me1.narrowPeak.gz" -wo | gzip > $output_dir/H3K4me1_$input_file_base_name
bedtools intersect -a $sorted_input_file -b "$anno_dir/$roadmap_biospecimen-H3K4me3.narrowPeak.gz" -wo | gzip > $output_dir/H3K4me3_$input_file_base_name
bedtools intersect -a $sorted_input_file -b "$anno_dir/$roadmap_biospecimen-H3K9ac.narrowPeak.gz" -wo | gzip > $output_dir/H3K9ac_$input_file_base_name
bedtools intersect -a $sorted_input_file -b "$anno_dir/$roadmap_biospecimen-H3K9me3.narrowPeak.gz" -wo | gzip > $output_dir/H3K9me3_$input_file_base_name

# echo  -e "$anno_dir/$roadmap_biospecimen-H3K9ac.narrowPeak.gz\n"
# echo  -e "$sorted_input_file\n"