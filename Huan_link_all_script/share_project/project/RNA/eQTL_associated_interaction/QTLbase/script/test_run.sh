# Rscript Point_plot_emplambda_relevance_cis_trans_glm_png.R
# echo -e "Point_plot_emplambda_relevance_cis_trans_glm_png\n"
# Rscript Point_plot_emplambda_relevance_cis_trans_glm_pdf.R
# echo -e "Point_plot_emplambda_relevance_cis_trans_glm_pdf\n"
# Rscript Point_plot_emplambda_relevance_cis_trans_loess_pdf.R
# echo -e "Point_plot_emplambda_relevance_cis_trans_loess_pdf\n"
# Rscript Point_plot_emplambda_relevance_trans_pdf_fitting.R
# echo -e "Point_plot_emplambda_relevance_trans_pdf_fitting\n"
# Rscript Point_plot_emplambda_relevance_trans_png_fitting.R
# echo -e "Point_plot_emplambda_relevance_trans_png_fitting\n"

# Rscript Point_plot_emplambda_relevance_cis_png_fitting.R
# echo -e "Point_plot_emplambda_relevance_cis_png_fitting\n"

# Rscript Point_plot_emplambda_relevance_cis_trans_loess_png.R
# echo -e "Point_plot_emplambda_relevance_cis_trans_loess_png\n"
# #---------------------------------------------------------




# Rscript Point_plot_emplambda_relevance_cis_pdf_fitting.R
# echo -e "Point_plot_emplambda_relevance_cis_pdf_fitting\n"
# Rscript Point_plot_emplambda_relevance_cis_trans_loess_pdf.R
# echo -e "Point_plot_emplambda_relevance_cis_trans_loess_pdf\n"


# perl 09_merge_emplambda_and_org_xQTL.pl 
# echo -e "09_merge_emplambda_and_org_xQTL\n"
# perl 09_merge_emplambda_and_org_all_QTLbase.pl
# echo -e "09_merge_emplambda_and_org_all_QTLbase\n"  
# perl 10_overlap_emplambda_xQTL.pl  
# echo -e "10_overlap_emplambda_xQTL\n"
# perl 10_overlap_emplambda_xQTL_cis.pl 
# echo -e "10_overlap_emplambda_xQTL_cis\n"
# perl 10_overlap_emplambda_xQTL_cis_trans.pl 
# echo -e "10_overlap_emplambda_xQTL_cis_trans\n"
# perl 10_overlap_emplambda_xQTL_trans_cis.pl 
# echo -e "10_overlap_emplambda_xQTL_trans_cis\n"
# perl 10_overlap_emplambda_xQTL_trans.pl
# echo -e "10_overlap_emplambda_xQTL_trans\n"


Rscript 11_Point_plot_emplambda_relevance_pdf.R 
echo -e "11_Point_plot_emplambda_relevance_pdf\n"
Rscript 11_Point_plot_emplambda_relevance_png.R #画../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_${QTL1}_${QTL2}.txt.gz中QTL1和QTL2的关联
echo -e "11_Point_plot_emplambda_relevance_png\n"
Rscript 11_Point_plot_emplambda_relevance_per_chr_pdf.R #分染色体画../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_${QTL1}_${QTL2}.txt.gz中QTL1和QTL2的关联
echo -e "11_Point_plot_emplambda_relevance_per_chr_pdf\n"
Rscript 11_Point_plot_emplambda_relevance_per_chr_png.R
echo -e "11_Point_plot_emplambda_relevance_per_chr_png\n"
Rscript 11_Point_plot_emplambda_relevance_cis_pdf.R #画../output/xQTL_merge/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_cis_10MB_${QTL1}_${QTL2}.txt.gz中QTL1和QTL2的关联
echo -e "11_Point_plot_emplambda_relevance_cis_pdf\n"
Rscript 11_Point_plot_emplambda_relevance_cis_png.R
echo -e "11_Point_plot_emplambda_relevance_cis_png\n"
Rscript 11_Point_plot_emplambda_relevance_trans_pdf.R
echo -e "11_Point_plot_emplambda_relevance_trans_pdf\n"
Rscript 11_Point_plot_emplambda_relevance_cis_trans_pdf.R
echo -e "11_Point_plot_emplambda_relevance_cis_trans_pdf\n"