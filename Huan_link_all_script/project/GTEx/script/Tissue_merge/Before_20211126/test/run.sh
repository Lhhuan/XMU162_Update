perl 07_test_tissue_qtl_diff.pl 
zless "./tmp_output/07_chr_muliti_tissue_eqtl.bed.gz" | sort -k1,1 -k2,2n |gzip >"./tmp_output/07_chr_muliti_tissue_eqtl_sorted.bed.gz"

bedtools intersect -a "./tmp_output/07_chr_hotspot.bed.gz" -b "./tmp_output/07_chr_muliti_tissue_eqtl_sorted.bed.gz" -wa -wb |gzip > ./tmp_output/07_chr_muliti_tissue_eqtl_hotspot.bed.gz


perl 08_check_loci_diff_cross_tissue.pl 

Rscript 09_find_chr1_snp_density.R 
