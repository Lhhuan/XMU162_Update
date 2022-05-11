# echo "01 start"
# perl 01_merge_all_tissue_egene.pl
# echo "01 finish"
# perl 02_filter_sig_egene.pl #
# echo "02 finish"
zless "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_0.05.bed.gz" |sort -k1,1 -k2,2n |gzip > "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_0.05_sorted.bed.gz"
echo "sort 0.05 finish"
zgrep -v "SNP_chr" "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_5e_8.bed.gz" |sort -k1,1 -k2,2n |gzip > "../output/02_merge_all_tissue_cis_eQTL_eur_egene_sig_5e_8_sorted.bed.gz"
echo "sort 5e-8 finish"