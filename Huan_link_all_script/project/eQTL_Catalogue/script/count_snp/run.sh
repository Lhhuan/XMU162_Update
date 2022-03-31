perl 01_transform_nhp_bed.pl #将../../output/all_tissue_status/NHP/NHPoisson_emplambda_interval_18_cutoff_7.3_Tissue_merge.txt.gz 转化为bed格式，得./output/01_nhp.bed.gz,排序得./output/01_nhp_sorted.bed.gz
echo "01 finish\n"
bedtools intersect -a "../../output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176.bed.gz" -b ./output/01_nhp.bed.gz  -wa -wb |gzip > ./output/hotspot_snp.bed.gz
echo "02 finish\n"
perl 02_count_snp_in_hotspot.pl #统计"./output/hotspot_snp.bed.gz" 包含的snps数目，得"./output/02_count_snp_in_hotspot.bed.gz"
zless ./output/02_count_snp_in_hotspot.bed.gz| sort -k1,1 -k2,2n |gzip >./output/02_count_snp_in_hotspot_sorted.bed.gz




Rscript 03_barplot_snp_number_in_hotspot_distribution.R 