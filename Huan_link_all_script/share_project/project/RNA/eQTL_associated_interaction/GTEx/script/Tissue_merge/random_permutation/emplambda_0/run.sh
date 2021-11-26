perl 01_filter_emplambda_0_segment.pl ###过滤出../../../../../output/Tissue_merge/Cis_eQTL/NHP/NHPoisson_emplambda_interval_${j}_cutoff_7.3_Tissue_merge.txt.gz 中emplambda$emplambda ==0，得 "../../../../../output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/sampling/Tissue_merge_segment_hotspot_cutoff_0.bed.gz"

zless "../../../../../output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/sampling/Tissue_merge_segment_hotspot_cutoff_0.bed.gz" |sort -k1,1 -k2,2n |sort -u |gzip "../../../../../output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/sampling/Tissue_merge_segment_hotspot_cutoff_0_sorted.bed.gz"
perl 02_random_genomic_resemble_hotspot.pl
perl 02_random_genomic_resemble_hotspot_filter_length_3103833.pl #hotspot filter  3103833 ,改片段太长，在随机中取不到
#----------------------此后都是对 filter length_3103833 进行处理
perl 03_annotation_marker.pl
perl 04_count_annotation_marker.pl
perl 05_calculate_jaccard_index_1000.pl