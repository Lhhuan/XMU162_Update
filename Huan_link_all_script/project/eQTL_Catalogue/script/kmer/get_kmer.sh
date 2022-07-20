
source activate seekr_source
seekr_kmer_counts $output_fa  -o ${i}_6mers_uc_us_no_log.csv --log2 none -uc -us
gzip ${i}_6mers_uc_us_no_log.csv
rm $output_fa