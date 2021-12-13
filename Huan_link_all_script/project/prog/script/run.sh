Rscript 01_add_feature_and_filter_grade.R #save(dat,file="01_add_age_raw_pfs_os_filter_grade.Rdata")

Rscript 01_plot_feature_distrubution.R
#Rscript 02_generate_feature_table_before_2021.12.9.R 
Rscript 02_generate_feature_table.R

#---------------
Rscript 03_cox_variable_analysis_not_fill_na_3a.R #01_add_age_raw_pfs_os_filter_grade_mergrPod_3A.txt
Rscript 03_cox_variable_analysis_fill_3a.R
Rscript 03_cox_variable_analysis_not_fill_na_all.R
Rscript 03_cox_variable_analysis_fill_all.R
#----------------------------

Rscript 03_plot_survival_curve.R
python 04_multiClassifier_fill.py
python 04_test_cv.py
python 05_parameter_optimization_XGB_fill.py
python 05_parameter_optimization_XGB_not_fill_na.py




python 06_XGBClassifier.py
python 06_XGBClassifier_fill.py
python 06_XGBClassifier_not_fill_na.py
python 06_XGBClassifier_not_fill_test.py
python 061_XGBClassifier_fill_feature_selection.py
python 061_XGBClassifier_not_fill_na_feature_selection.py
python 062_XGBClassifier_not_fill_gain_weight_shap.py



#------------------
Rscript 063_plot_feature_importance_not_fill_gain_weight_shap.R
python 064_XGBClassifier_not_fill_na_feature_selection_weight.py #---i
Rscript 065_plot_feature_importance_not_fill_gain_weight_shap_step3.R #i 
Rscript 065_plot_feature_importance_not_fill_gain_weight_shap_step2.R
Rscript 066_plot_feature_distrbution.R
Rscript 066_plot_feature_distrbution_new_Pod_total.R
Rscript 066_plot_feature_distrbution_new_Pod_total_refine.R
#--------------------------------
python 067_XGBClassifier_not_fill_gain_weight_shap.py #iiiiiii
Rscript 068_plot_feature_importance_not_fill_gain_weight_shap.R #iiiii
Rscript 069_plot_feature_distrbution_new_Pod_total.R
Rscript 069_plot_feature_distrbution_new_Pod_total_train_test.R
#
Rscript 069_plot_feature_distrbution_new_Pod_total_train_test_refine.R #部分修图

Rscript 070_calculate_and_plot_3a_fpb.R  ##07_add_model_index_age_raw_pfs_os_filter_grade_mergrPod_3A.Rdata
Rscript 071_calculate_and_plot_all_fpb.R #071_add_model_index_age_raw_pfs_os_filter_grade_mergrPod_all.Rdata
Rscript 080_xxxxx.R
Rscript 080_predict_and_plot_roc_survival.R
#--------------------------------------------muli
Rscript 090_exact_train_dataset_for_model.R  #09_test_train_dataset.Rdata 09_train_dataset.txt
python 091_XGBClassifier_not_fill_gain_weight_shap.py
Rscript 092_plot_feature_importance_not_fill_gain_weight_shap.R #----

Rscript 093_new_Pod_total_test_dataset.R 
Rscript 094_new_Pod_total_test_dataset_cutoff2.R
#----------------------- 2class


#
python 10_parameter_optimization_2XGB_not_fill_na.py
python 10_1_2XGB_not_fill_gain_weight_shap.py
python 10_1_2XGB_not_fill_gain_weight_shap_all_features.py
Rscript 10_2_plot_feature_importance_not_fill_gain_weight_shap.R #10_2_feature_ratio.txt
Rscript 10_2_plot_feature_importance_not_fill_gain_weight_shap_refine_figure.R #修图
Rscript 10_2_plot_feature_importance_not_fill_gain_weight_shap_all_features.R
Rscript 10_3_new_Pod_total_test_dataset.R
Rscript 10_3_new_Pod_total_test_dataset_refine_figure.R  #修图
Rscript 10_3_new_Pod_total_test_dataset_refine_figure_2.R #importance 10_3_train_dataset_new_cutoff.txt
Rscript 10_4_cox_variable_and_plot.R 

Rscript 10_5_feature_chi_test_fill.R
Rscript 10_5_feature_chi_test_na.R





Rscript check_FLIPI.R 
check_and_plot_FLIPI1_all_sample.R
original_FLIPI.R
test2.R
test.R
test_plip.R