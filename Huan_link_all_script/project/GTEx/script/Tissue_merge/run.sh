Rscript 01_fiter_hotspot_for_homer_and_bar_plot_distribution.R 
#./figure/01_hotspot_distribution.txt,../../output/Tissue_merge/Cis_eQTL/01_hotspot_length_small_than_2290.bed
#---------homer
    source activate huan_py3
    findMotifsGenome.pl ../../output/Tissue_merge/Cis_eQTL/01_hotspot_length_small_than_2290.bed hg19 /home/huanhuan/project/GTEx/output/Tissue_merge/Cis_eQTL/homer/ -size given

#----------

perl 02_calculate_marker_annotation_fraction.pl #"../../output/Tissue_merge/Cis_eQTL/02_hotspot_mark_annotatiob_fraction.txt.gz"

bedtools nuc -fi "/share/Projects/huanhuan/ref_data/gencode/GRCh37.primary_assembly.genome.fa" -bed /share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge_filter_3103833.bed.gz |gzip >/home/huanhuan/project/GTEx/output/Tissue_merge/Cis_eQTL/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge_filter_3103833_GC.bed.gz

Rscript 03_plot_gc_content_distrbution.R 
Rscript 04_pca_cluster.R  #04_10_markers_euclidean_dist.Rdata
Rscript 04_kpca_cluster.R
Rscript 04_Kernel_pca_cluster.R
Rscript 05_kmeans.R 
Rscript 05_kkl.R 
/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge_filter_3103833.bed.gz