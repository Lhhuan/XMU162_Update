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

# perl 021_judge_kmer_exists_and_create.pl #

# Rscript 03_identify_wilcox_significant_kmer.R
Rscript 03_identify_wilcox_significant_kmer_two_side.R
Rscript 04_count_kmer_count.R 
Rscript 05_plot_top30_kmer_distribution.R
Rscript 06_identify_permutation_significant_kmer.R
Rscript 07_plot_top30_kmer_distribution.R
Rscript 08_exact_sigkmer_whole_genome.R
Rscript 09_Cluster_hotspot_by_sig_kmer_kpca_chr_per.R
Rscript 10_chr1_split_rbfdot.R
Rscript 11_chr1_kpca_louvain.R
source activate huan_py3 
perl 12_chr1_homer.pl 
perl 12_chr1_homer_200.pl
perl 12_chr1_all_homer.pl
conda deactivate 
perl 12_01_filter_homer200_result.pl #
source activate homer_hg38
PATH=$PATH:/home/huanhuan/tools/homer_hg38/.//bin/
perl 12_chr1_homer_200_new38.pl
conda deactivate 
perl 12_1_get_fa_for_meme.pl 



Rscript 13_count_cluster_specific_gene.R

perl 14_chr1_gc_content.pl #
Rscript 15_boxplot_gc_content.R
Rscript 16_cluster_marker.R
perl 16_1_plot_marker_peak.pl #

Rscript 16_2_marker_fisher_test.R
echo -e "chr1\t248956422" > /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/markers/bedtools_fisher/t.genome
Rscript 16_2_marker_bedtools_fisher_test.R

perl 17_merge_gwas_and_hotspot.pl 
Rscript 18_boxplot_gwas.R 

