
#==========================
bedtools getfasta -fi "/share/Projects/huanhuan/ref_data/gencode/GRCh38.primary_assembly.genome.fa" -bed predicted_regions_win5000_large_than6.bed -fo "01_predicted_regions.fa"
# cd ../output/

source activate seekr_source
seekr_kmer_counts 01_predicted_regions.fa  -o  03_predicted_regions_6mers_uc_us_no_log.csv --log2 none -uc -us
gzip 03_predicted_regions_6mers_uc_us_no_log.csv
rm 01_predicted_regions.fa

conda deactivate 
bedtools nuc -fi /share/Projects/huanhuan/ref_data/gencode/GRCh38.primary_assembly.genome.fa -bed predicted_regions_win5000_large_than6.bed |cut -f1-3,5 |gzip >03_predicted_regions_gc_content.bed.gz