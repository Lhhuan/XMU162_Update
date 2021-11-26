# perl 01_merge_all_kinds_QTL.pl 
# echo -e "01_merge_all_kinds_QTL.pl\n"
perl 0014_merge_pop_tissue.pl #画../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_${QTL1}_${QTL2}.txt.gz中QTL1和QTL2的关联
echo -e "merge_pop_tissue.pl\n"
perl 0015_judge_xQTL_cis_trans.pl #分染色体画../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_${QTL1}_${QTL2}.txt.gz中QTL1和QTL2的关联
echo -e "judge_xQTL_cis_trans.pl\n"
# Rscript 11_Point_plot_emplambda_relevance_per_chr_png.R
# echo -e "11_Point_plot_emplambda_relevance_per_chr_png\n"
# Rscript 11_Point_plot_emplambda_relevance_cis_pdf.R #画../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_cis_10MB_${QTL1}_${QTL2}.txt.gz中QTL1和QTL2的关联
# echo -e "11_Point_plot_emplambda_relevance_cis_pdf\n"
# Rscript 11_Point_plot_emplambda_relevance_cis_png.R
# echo -e "11_Point_plot_emplambda_relevance_cis_png\n"
# Rscript 11_Point_plot_emplambda_relevance_trans_pdf.R
# echo -e "11_Point_plot_emplambda_relevance_trans_pdf\n"
# Rscript 11_Point_plot_emplambda_relevance_cis_trans_pdf.R
# echo -e "11_Point_plot_emplambda_relevance_cis_trans_pdf\n"