
perl 12_3pcn_whole_genome_homer_200.pl
perl 12_01_filter_homer200_result.pl
Rscript 12_02_venn_homer200.R
# Rscript 12_03_cluster_specific_homer.R

Rscript 13_count_cluster_specific_eqtl_gene.R
Rscript 13_1_cluster_specific_egene.R
Rscript 13_2_cluster_tissue-speific_analysis.R 
Rscript 13_3_eqtl_count.R 

perl 14_gc_content.pl
Rscript 15_boxplot_gc_content.R
perl 16_1_plot_marker_peak_point.pl
perl 17_merge_gwas_and_hotspot.pl

Rscript 18_boxplot_gwas_trait.R
Rscript 18_boxplot_gwas_loci.R
Rscript 19_barplot_gc_island.R
Rscript 20_methylation.R
Rscript 20_methylation_huan_mean.R
# Rscript 21_whole_genome_kpca_louvain_leiden_pcn3_adjust_umap.R

bedtools intersect -a ../../../output/figure/whole_genome/3pca_1e_4/GWAS/17_whole_genome_leiden_pca3_k50_resolution1e-04_sorted.bed.gz -b "/share/data0/QTLbase/huan/Cis_Regulatory_Elements/encodeCcreCombined_sorted.bed.gz" -wa -wb |gzip > ../../../output/figure/whole_genome/3pca_1e_4/Cis_Regulatory_Elements/hotspot_cluster_cCREs.bed.gz
Rscript 22_barplot_cCREs.R 
Rscript 22_1_upsetplot_cCREs.R 
Rscript 22_2_heatmap_cCREs.R 
# scp "/share/data0/QTLbase/huan/Conservation/hg38.phastCons100way_transform.bed.gz" huan@202.113.53.212:/h/huan/project/huan_tmp/


bedtools intersect -a "/share/data0/QTLbase/huan/Conservation/hg38.phastCons100way_transform.bed.gz" -b ../../../output/figure/whole_genome/3pca_1e_4/GWAS/17_whole_genome_leiden_pca3_k50_resolution1e-04_sorted.bed.gz  -wa -wb |gzip > ../../../output/figure/whole_genome/3pca_1e_4/Conservation/hotspot_cluster_phastCons100way.bed.gz
Rscript 23_barplot_Conservation_phastCons100way.R 
bedtools intersect -a "/share/data0/QTLbase/huan/Conservation/hg38.phyloP100way_transform.bed.gz" -b ../../../output/figure/whole_genome/3pca_1e_4/GWAS/17_whole_genome_leiden_pca3_k50_resolution1e-04_sorted.bed.gz  -wa -wb |gzip > ../../../output/figure/whole_genome/3pca_1e_4/Conservation/hotspot_cluster_phyloP100way.bed.gz

Rscript 23_boxplot_Conservation_phastCons100way.R

perl 24_adjust_hic_format_intra_chr.pl
perl 24_adjust_hic_format_inter_chr.pl
Rscript 25_boxplot_chr_intra.R
Rscript 25_boxplot_chr_inter.R
#-----------------------------------
perl 26_adjust_gwas_catalogy.pl 
Rscript 27_boxplot_gwas.R
Rscript 27_boxplot_gwas_loci.R

Rscript 28_circos_density.R


wget -c http://www.oreganno.org/dump/ORegAnno_Combined_2016.01.19.tsv


