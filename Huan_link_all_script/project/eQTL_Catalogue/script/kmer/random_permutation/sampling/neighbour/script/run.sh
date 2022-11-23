# perl 01_extend_5_fold_neighbour.pl 
perl 01_extend_10_fold_neighbour.pl
perl 02_random_10_fold_neighbour1000.pl
#----------
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
scp -r huan@202.113.53.212:"/h/huan/huan_k/whole_genome_rbfdot_pca.Rdata" ./
cp whole_genome_rbfdot_pca.Rdata /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/

Rscript 11_whole_genome_kpca_louvain_leiden.R
source activate huan_py3 
perl 12_whole_genome_homer_given.pl 
perl 12_whole_genome_homer_200.pl
perl 12_umap_whole_genome_homer_200.pl
conda deactivate 
perl 12_01_filter_homer200_result.pl #
Rscript 12_02_venn_homer200.R
Rscript 12_03_cluster_specific_homer.R


Rscript 13_count_cluster_specific_gene.R

perl 14_gc_content.pl #
Rscript 15_boxplot_gc_content.R
Rscript 16_cluster_marker.R
perl 16_1_plot_marker_peak_point.pl
perl 16_1_umap_plot_marker_peak_point.pl

Rscript 16_2_marker_fisher_test.R
echo -e "chr1\t248956422" > /home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/markers/bedtools_fisher/t.genome
Rscript 16_2_marker_bedtools_fisher_test.R

perl 17_merge_gwas_and_hotspot.pl 
Rscript 18_boxplot_gwas.R 
Rscript 19_barplot_gc_island.R 




