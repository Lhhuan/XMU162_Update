perl 01_extract_data.pl #
bedtools getfasta -fi ~/ref_data/gencode/GRCh38.primary_assembly.genome.fa -bed 01_filter_nodes.bed.gz -fo nodes.fa



source activate /home/huanhuan/miniconda3
seekr_kmer_counts nodes.fa -o 6mers_uc_us_no_log.csv --log2 none -uc -us