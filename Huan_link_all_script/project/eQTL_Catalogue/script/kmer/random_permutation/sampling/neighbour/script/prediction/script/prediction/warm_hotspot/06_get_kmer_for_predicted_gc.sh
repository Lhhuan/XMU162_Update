
#==========================
bedtools getfasta -fi "/share/Projects/huanhuan/ref_data/gencode/GRCh38.primary_assembly.genome.fa" -bed warm_hotspot_region_win4518_large_than6.bed.gz -fo "06_predicted_warm_hotspot_regions.fa"
# cd ../output/

source activate seekr_source
seekr_kmer_counts 06_predicted_warm_hotspot_regions.fa  -o  06_predicted_warm_hotspot_regions_6mers_uc_us_no_log.csv --log2 none -uc -us
gzip 06_predicted_warm_hotspot_regions_6mers_uc_us_no_log.csv
rm 06_predicted_warm_hotspot_regions.fa

conda deactivate 
bedtools nuc -fi /share/Projects/huanhuan/ref_data/gencode/GRCh38.primary_assembly.genome.fa -bed warm_hotspot_region_win4518_large_than6.bed.gz |cut -f1-3,5 |gzip >06_predicted_warm_hotspot_regions.gc_content.bed.gz