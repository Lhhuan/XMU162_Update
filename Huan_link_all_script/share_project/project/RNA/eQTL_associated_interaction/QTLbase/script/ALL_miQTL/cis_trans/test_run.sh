# perl 01_Completion_snp_for_cis_trans_miQTL_by_1kg.pl | sort -u
# echo -e "01_Completion_snp_for_cis_trans_miQTL_by_1kg\n"
Rscript 02_huan_NHPoisson.R
echo -e "huan_NHPoisson\n"
perl 03_filter_emplambda_in_miQTL.pl
echo -e "03_filter_emplambda_in_miQTL\n"