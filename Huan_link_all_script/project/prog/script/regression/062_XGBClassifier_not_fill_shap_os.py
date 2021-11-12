import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from sklearn.metrics import plot_roc_curve
from sklearn import tree
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_validate, cross_val_predict
from xgboost import XGBRegressor
import xgboost as xgb
from sklearn.metrics import accuracy_score
from sklearn.feature_selection import SelectFromModel
from numpy import sort
import shap
#----------------------------------------------------

# dataset = pd.read_table("../output/01_add_age_raw_pfs_os_filter_grade_mergrPod_3A_fillna.txt")
dataset = pd.read_table("../../output/01_add_age_raw_pfs_os_filter_grade_mergrPod_3A.txt")

X = dataset.loc[:, ['gender','Ki.67','stage','Bsym','LN_num','site0','extend','BM','spleen','extend_num','BM_extend','LN6','SUVmax','SPD','X150b2mg_ldh','b2mg_LDH','ECOG','B2mg','LDH','LDH0','WBC','HGB','HGB0','PLT','Mono','Lym','age_raw','age_raw_60','ki67_20','LN_num_6','extend_num_0','SUVmax_2','LDH_300','Lym_Mono_10','B2mg_3.4','SPD_0']]
y = dataset.loc[:, 'OS_month']


model = XGBRegressor(booster='gbtree',
                    eval_metric='rmse',
                    gamma = 0.1,
                    min_child_weight= 1,
                    max_depth= 5,
                    subsample= 0.8,
                    colsample_bytree= 0.8,
                    tree_method= 'exact',
                    learning_rate=0.05,
                    n_estimators=100,
                    nthread=10,
                    reg_alpha = 2e-05,
                    seed=27)

model.fit(X,y)


explainer = shap.TreeExplainer(model)
shap_values = explainer(X)
shap.plots.bar(shap_values)
# shap.plots.waterfall(shap_values)
plt.savefig('output/figure/062_os_feature_shap_bar.pdf')
plt.show()
# fig, ax = plt.subplots(figsize=(15, 15))

feature_shap =pd.DataFrame(shap_values.values)
feature_shap.columns =X.columns
feature_shap.to_csv("output/062_not_fill_feature_importance_shap_os.txt",sep="\t",index=None)