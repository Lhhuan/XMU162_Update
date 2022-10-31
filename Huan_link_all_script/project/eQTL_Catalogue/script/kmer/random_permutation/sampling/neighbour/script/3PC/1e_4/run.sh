perl 12_3pcn_whole_genome_homer_200.pl
perl 12_01_filter_homer200_result.pl
Rscript 12_02_venn_homer200.R
# Rscript 12_03_cluster_specific_homer.R

Rscript 13_count_cluster_specific_eqtl_gene.R
Rscript 13_1_cluster_specific_egene.R
Rscript 13_2_cluster_tissue-speific_analysis.R 

perl 14_gc_content.pl
Rscript 15_boxplot_gc_content.R
perl 16_1_plot_marker_peak_point.pl
perl 17_merge_gwas_and_hotspot.pl

Rscript 18_boxplot_gwas.R
Rscript 19_barplot_gc_island.R
Rscript 20_methylation.R
Rscript 20_methylation_huan_mean.R
Rscript 21_whole_genome_kpca_louvain_leiden_pcn3_adjust_umap.R

bedtools intersect -a ../../../output/figure/whole_genome/3pca_1e_4/GWAS/17_whole_genome_leiden_pca3_k50_resolution1e-04_sorted.bed.gz -b "/share/data0/QTLbase/huan/Cis_Regulatory_Elements/encodeCcreCombined_sorted.bed.gz" -wa -wb |gzip > ../../../output/figure/whole_genome/3pca_1e_4/Cis_Regulatory_Elements/hotspot_cluster_cCREs.bed.gz
Rscript 22_barplot_cCREs.R 