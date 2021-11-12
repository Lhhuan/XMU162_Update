#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np

# dataset = pd.read_table("../output/01_add_age_raw_pfs_os_filter_grade_mergrPod_3A.txt")
dataset = pd.read_table("../output/01_add_age_raw_pfs_os_filter_grade_mergrPod_3A_fillna.txt")
dataset.head
# y = dataset.iloc[:, -1].values
X = dataset.loc[:, ['gender','Ki.67','stage','Bsym','LN_num','site0','extend','BM','spleen','extend_num','BM_extend','LN6','SUVmax','SPD','X150b2mg_ldh','b2mg_LDH','ECOG','B2mg','LDH','LDH0','WBC','HGB','HGB0','PLT','Mono','Lym','interm_res','end_res','CR','Rmaintain','trans','age_raw','age_raw_60','ki67_20','LN_num_6','extend_num_0','stage_III_IV','stage_I_II','SUVmax_2','LDH_300','Lym_Mono_10','B2mg_3.4','SPD_0']].values


y = dataset.loc[:, 'pod_total'].values

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
from sklearn.metrics import accuracy_score
cv = StratifiedKFold(n_splits=5)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

# model = RandomForestClassifier(n_estimators=100,random_state=1)
# model =tree.DecisionTreeClassifier()

# model.fit(X_train,y_train)
# y_predict=model.predict(X_test)
# accuracy_score(y_test,y_predict)

# model = XGBRegressor()
# cv = KFold(n_splits=5, shuffle=True, random_state=0)
# y_pred = cross_val_predict(model, X, y,cv=cv)
# accuracy_score(y_test,y_predict)
# model = RandomForestClassifier(n_estimators=100,random_state=1)
# auc_scores = np.mean(cross_val_score(model, X, y, scoring='roc_auc', cv=5))
model = RandomForestClassifier(n_estimators=100,random_state=1)
accuracy_scores = np.mean(cross_val_score(model, X, y, scoring='accuracy', cv=5))
#0.8379343816734615
#0.782644146767618
model =tree.DecisionTreeClassifier()
accuracy_scores = np.mean(cross_val_score(model, X, y, scoring='accuracy', cv=5))
#0.8043486701611338
#0.6026596777324792

model =XGBClassifier(use_label_encoder=False,objective = "multi:softprob",eval_metric= "mlogloss")
accuracy_scores = np.mean(cross_val_score(model, X, y, scoring='accuracy', cv=5))
#0.8418753640069889
#0.7569209862162687