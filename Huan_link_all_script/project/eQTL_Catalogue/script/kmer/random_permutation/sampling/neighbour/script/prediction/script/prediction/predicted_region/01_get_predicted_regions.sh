zcat "/home/huanhuan/project/eQTL_Catalogue/output/sampling/Tissue_merge_segment_hotspot_cutoff_0_0.176.bed.gz" "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz"|sort -k1,1 -k2,2n |gzip > all_region_sorted.bed.gz
bedtools merge -i all_region_sorted.bed.gz |gzip >all_region_sorted_merge.bed.gz

bedtools subtract -a "/home/huanhuan/project/eQTL_Catalogue/output/sampling/Tissue_merge_segment_hotspot_cutoff_0.bed.gz" -b "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" |gzip > pure_cold_region.bed.gz

bedtools subtract -a all_region_sorted_merge.bed.gz -b pure_cold_region.bed.gz |sort -k1,1 -k2,2n |gzip > all_region_without_cold_sorted.bed.gz

bedtools merge -i all_region_without_cold_sorted.bed.gz |gzip >all_region_without_cold_sorted_merge.bed.gz

bedtools makewindows -b all_region_without_cold_sorted_merge.bed.gz -w 5000 |gzip > predicted_regions_win5000.bed.gz