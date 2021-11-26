# Rscript 11_Point_plot_emplambda_relevance_cis_pdf.R #画../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_cis_10MB_${QTL1}_${QTL2}.txt.gz中QTL1和QTL2的关联
# echo -e "11_Point_plot_emplambda_relevance_cis_pdf\n"
# Rscript 11_Point_plot_emplambda_relevance_cis_png.R
# echo -e "11_Point_plot_emplambda_relevance_cis_png\n"
# Rscript 11_Point_plot_emplambda_relevance_trans_pdf.R
# echo -e "11_Point_plot_emplambda_relevance_trans_pdf\n"
# Rscript 11_Point_plot_emplambda_relevance_cis_trans_pdf.R
# echo -e "11_Point_plot_emplambda_relevance_cis_trans_pdf\n"


perl 10_overlap_emplambda_eQTL_eQTL.pl
echo -e "10_overlap_emplambda_eQTL_eQTL.pl\n"
Rscript 11_Point_plot_eQTL_cis_trans_emplambda_relevance_pdf.R
echo -e "11_Point_plot_eQTL_cis_emplambda_relevance_pdf\n"
