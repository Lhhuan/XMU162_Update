# perl 01_extend_5_fold_neighbour.pl 
perl 01_extend_10_fold_neighbour.pl
perl 02_random_10_fold_neighbour1000.pl


####get kmer count 
    bedtools getfasta -fi "/share/Projects/huanhuan/ref_data/gencode/GRCh38.primary_assembly.genome.fa" -bed "../output/1_10_fold_neighbour_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" -fo "../output/1_10_fold_neighbour_resemble.fa"
    cd ../output/
    source activate seekr_source
    seekr_kmer_counts 1_10_fold_neighbour_resemble.fa  -o  6mers_uc_us_no_log.csv --log2 none -uc -us
    gzip 6mers_uc_us_no_log.csv
    cd ../script
#----------
perl 02_random_10_fold_neighbour100.pl
perl 021_get_kmer_count.pl
perl 021_get_kmer_count1_100.pl
perl 021_get_kmer_count101_200.pl
perl 021_get_kmer_count201_300.pl
perl 021_get_kmer_count301_400.pl
perl 021_get_kmer_count401_500.pl
perl 021_get_kmer_count501_600.pl
perl 021_get_kmer_count601_700.pl
perl 021_get_kmer_count701_800.pl
perl 021_get_kmer_count801_900.pl
perl 021_get_kmer_count901_1000.pl

perl 021_judge_kmer_exists_and_create.pl #

Rscript 03_identify_wilcox_significant_kmer.R
Rscript 03_identify_wilcox_significant_kmer_two_side.R
Rscript 04_count_kmer_count.R 
Rscript 05_plot_top30_kmer_distribution.R
Rscript 06_identify_permutation_significant_kmer.R
Rscript 07_plot_top30_kmer_distribution.R

