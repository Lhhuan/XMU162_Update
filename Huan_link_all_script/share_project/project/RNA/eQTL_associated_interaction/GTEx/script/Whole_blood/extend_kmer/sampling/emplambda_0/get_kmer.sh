
source activate /home/huanhuan/miniconda3
seekr_kmer_counts $input  -o ${i}_6mers_uc_us_no_log.csv --log2 none -uc -us
gzip ${i}_6mers_uc_us_no_log.csv
source deactivate