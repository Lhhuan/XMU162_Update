perl 01_merge_all_tissue_eQTL.pl #
Rscript 02_NHP_big_par.R
perl 03_filter_hotspot_for_interval18.pl
Rscript 04_circos_density.R
mv *.pdf ./figure
Rscript 05_find_chr1_max_density.R 
perl 06_merge_merge_hotspot_and_egene.pl 
