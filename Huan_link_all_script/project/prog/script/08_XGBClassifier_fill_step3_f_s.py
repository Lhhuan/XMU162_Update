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
from xgboost import XGBClassifier
import xgboost as xgb
from sklearn.metrics import accuracy_score

#----------------------------------------------------

dataset = pd.read_table("../output/01_add_age_raw_pfs_os_filter_grade_mergrPod_3A_fillna.txt")
# dataset = pd.read_table("../output/01_add_age_raw_pfs_os_filter_grade_mergrPod_3A.txt")

X = dataset.loc[:, ['gender','Ki.67','stage','Bsym','LN_num','extend','BM','extend_num','BM_extend','LN6','SUVmax','SPD','b2mg_LDH','ECOG','B2mg','LDH','LDH0','WBC','HGB','PLT','Mono','Lym','interm_res','end_res','Rmaintain','trans','age_raw']].values

y = dataset.loc[:, 'pod_total'].values

# model=XGBClassifier(reg_alpha = 1e-05, 
#     colsample_bytree= 0.65, 
#     subsample= 0.5,
#     gamma = 1.7,
#     min_child_weight =3,
#     max_depth= 3, 
#     objective = "multi:softprob",
#     eval_metric= "mlogloss",
#     use_label_encoder=False,
#     learning_rate = 0.1,
#     n_estimators=50, 
#     seed=1024)

model =XGBClassifier(reg_alpha=0.01,
    colsample_bytree = 0.95,
    subsample = 0.9, 
    gamma = 0.9,
    min_child_weight =3,
    max_depth= 3, 
    objective = "multi:softprob",
    eval_metric= "mlogloss",
    use_label_encoder=False,
    learning_rate = 0.1,
    n_estimators=50, 
    seed=1024)


cv = KFold(n_splits=10, shuffle=True, random_state=0)
y_pred = cross_val_predict(model, X, y,cv=cv)
dataset['predict_pod_total']=pd.DataFrame(y_pred)
accuracy_score(y,y_pred)
dataset.to_csv("../output/08_XGBClassifier_predict_3a_fill_CV_step3_f_s.txt",sep="\t",index=None)

# import xgboost as xgb