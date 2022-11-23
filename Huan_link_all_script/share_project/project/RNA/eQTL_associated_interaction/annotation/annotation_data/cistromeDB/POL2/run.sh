perl 05_filter_normal_chrom_info.pl
bedtools merge -i ../normal_cell/Human_FACTOR/POL2/merge_pos_info_narrow_peak_sorted.bed.gz |gzip > ../normal_cell/Human_FACTOR/POL2/merge_pos_info_narrow_peak_sort_union.bed.gz
Rscript 07_intersect_merge_bed_and_signalValue.R