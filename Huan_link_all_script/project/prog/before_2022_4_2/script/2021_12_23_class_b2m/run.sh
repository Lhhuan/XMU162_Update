Rscript filter.R #./data/rawdata_filter.csv
Rscript 01_add_feature_and_filter_grade.R #save(dat,file="01_add_age_raw_pfs_os_filter_grade.Rdata")
#---------------
Rscript 03_cox_variable_analysis_not_fill_na_3a.R #01_add_age_raw_pfs_os_filter_grade_mergrPod_3A.txt
Rscript 03_cox_variable_analysis_fill_3a.R #01_add_age_raw_pfs_os_filter_grade_mergrPod_3A_fillna.txt

#----------------------------

Rscript 03_plot_survival_curve.R
Rscript 04_refine_ldh_b2m.R #04_add_age_raw_pfs_os_filter_grade_refine_ldh_b2m.Rdata
#--------------------------------
Rscript 05_feature_chi_test_na_fill.R
Rscript 06_cox_variable_and_plot.R  #06_3a_file.Rdata
Rscript 070_calculate_and_plot_3a_fpb.R  ##07_add_model_index_age_raw_pfs_os_filter_grade_mergrPod_3A.Rdata

#--------------------------------------------
Rscript 090_extract_train_dataset_for_model.R  #09_test_train_dataset.Rdata 09_train_dataset.txt




#----------------------- 2class

python 10_parameter_optimization_2XGB_not_fill_na.py
python 10_1_2XGB_not_fill_gain_weight_shap.py
Rscript 10_2_plot_feature_importance_not_fill_gain_weight_shap.R #10_2_feature_ratio.txt
python 10_03_feature_selection_2XGB_not_fill_gain_weight_shap.py
Rscript 10_04_plot_feature_importance_not_fill_gain_weight_shap.R

Rscript 10_3_new_Pod_total_test_dataset.R #importance 10_3_train_dataset_new_cutoff.txt