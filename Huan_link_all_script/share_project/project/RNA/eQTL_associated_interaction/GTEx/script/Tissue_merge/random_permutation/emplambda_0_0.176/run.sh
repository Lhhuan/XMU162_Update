perl 01_filter_emplambda_0_0.176_segment_fact.pl ###过滤出"../../../../output/${tissue}/Cis_eQTL/NHP/NHPoisson_emplambda_interval_${j}_cutoff_7.3_${tissue}.txt.gz" 中emplambda$emplambda <0.176，得 "../../../../output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/sampling/Tissue_merge_segment_hotspot_cutoff_0_0.176.bed.gz"
perl 02_random_genomic_resemble_hotspot_fact.pl #产生10000个与"../../../../output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz"相同的resemble hotspot,"$output_dir/${i}_resemble_${input_file_base_name}"
perl 02_random_genomic_resemble_hotspot_fact_filter_length_3103833.pl #"../../../../output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge_filter_3103833.bed.gz"; #hotspot filter  3103833 ,该片段太长，随机在背景中取不到，去掉离群点再random 
#----------------------此后都是对 filter length_3103833 进行处理
perl 03_annotation_marker.pl
perl 04_count_annotation_marker.pl
perl 05_calculate_jaccard_index_1000.pl