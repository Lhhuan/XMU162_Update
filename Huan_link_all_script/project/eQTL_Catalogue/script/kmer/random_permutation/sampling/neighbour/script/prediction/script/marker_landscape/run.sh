cp "/share/data0/QTLbase/huan/eQTL_Catalogue/original_random/extend_18snp/emp0.176/0/1_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" ./output/
gunzip ./output/1_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz

cp "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz" ./output/

gunzip ./output/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz
perl 01_computer_marker_matrix.pl #
perl 01_computer_marker_matrix_na0.pl #

Rscript 02_exact_training_feature.R
Rscript 02_exact_training_feature_miss_as_0.R
Rscript 02_exact_training_feature_miss_as_0_mean.R

source activate huan_py3
python 03_marker_matrix_miss_as_na.py
python 03_marker_matrix_miss_as_0.py
python 03_marker_matrix_miss_as_0_mean.py
python 03_all_feature_marker_matrix_miss_as_0_mean.py

#====================
python 031_AdaBoostClassifier.py
python 031_DecisionTreeClassifier.py
python 031_RidgeClassifier.py
python 031_SGDClassifier.py
python 031_randomforest.py

python 051_final_XGB_class_feature.py
python 05_final_XGB_model_features.py