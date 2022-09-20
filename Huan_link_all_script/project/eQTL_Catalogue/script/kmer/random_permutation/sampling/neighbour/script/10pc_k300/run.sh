source activate homer_hg38
PATH=$PATH:/home/huanhuan/tools/homer_hg38/.//bin/
perl 12_chr1_homer_200_new38.pl
conda deactivate 
perl 12_1_get_fa_for_meme.pl 



Rscript 13_count_cluster_specific_gene.R

perl 14_chr1_gc_content.pl #
Rscript 15_boxplot_gc_content.R
Rscript 16_cluster_marker.R
perl 17_merge_gwas_and_hotspot.pl 
Rscript 18_boxplot_gwas.R 

