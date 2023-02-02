zless ../../../output/warm_region.bed.gz >warm_hotspot_regions_in_pantissue.bed
zless "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" >> warm_hotspot_regions_in_pantissue.bed
less warm_hotspot_regions_in_pantissue.bed |sort -k1,1 -k2,2n |gzip >warm_hotspot_regions_in_pantissue_sorted.bed.gz
bedtools merge -i warm_hotspot_regions_in_pantissue_sorted.bed.gz |gzip > warm_hotspot_regions_in_pantissue_sorted_merge.bed.gz
bedtools makewindows -b warm_hotspot_regions_in_pantissue_sorted_merge.bed.gz -w 4518 |gzip > warm_hotspot_regions_in_pantissue_sorted_merge_win4518.bed.gz

Rscript 0601_filter_lsmall6.R

bash 06_get_kmer_for_predicted_gc.sh

Rscript 07_exact_gc_kmer_for_predicted.R