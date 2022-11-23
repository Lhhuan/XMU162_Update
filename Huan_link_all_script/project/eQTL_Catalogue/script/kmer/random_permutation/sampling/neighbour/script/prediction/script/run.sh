bash 01_get_kmer_for_negative_gc.sh
perl 011_annotation_for_negative.pl #
Rscript 012_mean_marker_signalValue_for_negative_hotspot.R
Rscript 02_exact_training_feature.R
    python 03_AdaBoostClassifier.py
    python 03_SGDClassifier.py
    python 03_DecisionTreeClassifier.py
    python 03_RidgeClassifier.py
    python 03_XGBoost_class.py
    python 03_randomforest.py
python 04_XGB_parameter_optimization.py
python 05_final_XGB_CV.py

#=====================

#=====================

