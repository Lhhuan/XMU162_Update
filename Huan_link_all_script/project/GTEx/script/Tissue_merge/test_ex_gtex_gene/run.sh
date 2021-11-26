perl 01_extract_eqtl_egene.pl
Rscript 01_extract_eqtl_egene.R
less R_Whole_Blood_cis_sig_eQTL_egene.txt |sort -k1,1 -k2,2n > R_Whole_Blood_cis_sig_eQTL_egene_sorted.txt

perl 02_filter_no_egene_hotspot.pl

perl 03_merge_emplamda_and_p.pl


zless ./output/03_qtl_p_emplamda.bed.gz |sort -k1,1 -k2,2n |gzip > ./output/03_qtl_p_emplamda_sorted.bed.gz

bedtools intersect -a ./output/hotspot_without_egene_sorted.bed.gz -b ./output/03_qtl_p_emplamda_sorted.bed.gz -wa -wb  > ./output/hotspot_without_egene_sorted_eQTL.bed