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
from sklearn.metrics import mean_squared_error,r2_score
from sklearn.feature_selection import SelectFromModel
from numpy import sort

#----------------------------------------------------

# dataset = pd.read_table("../output/01_add_age_raw_pfs_os_filter_grade_mergrPod_3A_fillna.txt")
dataset = pd.read_table("../../output/01_add_age_raw_pfs_os_filter_grade_mergrPod_3A.txt")

X = dataset.loc[:, ['gender','Ki.67','stage','Bsym','LN_num','site0','extend','BM','spleen','extend_num','BM_extend','LN6','SUVmax','SPD','X150b2mg_ldh','b2mg_LDH','ECOG','B2mg','LDH','LDH0','WBC','HGB','HGB0','PLT','Mono','Lym','age_raw','age_raw_60','ki67_20','LN_num_6','extend_num_0','SUVmax_2','LDH_300','Lym_Mono_10','B2mg_3.4','SPD_0']]
y = dataset.loc[:, 'pfs_month']
X_COL = X.columns
# X= dataset.iloc[:,1:540]
# y = dataset.loc[:, 'AUC'].values
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=7)

model = XGBRegressor(booster='gbtree',
                    eval_metric='rmse',
                    gamma = 0.4,
                    min_child_weight= 1,
                    max_depth= 9,
                    subsample= 0.9,
                    colsample_bytree= 0.2,
                    tree_method= 'exact',
                    learning_rate=0.1,
                    n_estimators=100,
                    nthread=10,
                    reg_alpha=1e-06,
                    seed=27)

model.fit(X_train,y_train)
fig, ax = plt.subplots(1, 1, figsize=(8, 13))
xgb.plot_importance(model, max_num_features=30, height=0.5, ax=ax)
plt.savefig('output/figure/06_feature_importance_3a_not_fill_pfs.pdf')

# make predictions for test data and evaluate
predictions = model.predict(X_test)
mse = mean_squared_error(y_test, predictions)
print("MSE:",mse)
print("r2:",r2_score(y_test, predictions))
#----------------------------------------------------------------
# Fit model using each importance as a threshold
X=X.values
y=y.values
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=7)
#----------------
thrs=[]
x_train_s=[]
scores=[]
thresholds = sort(model.feature_importances_)
for thresh in thresholds:
	# select features using threshold
	selection = SelectFromModel(model, threshold=thresh, prefit=True)
	select_X_train = selection.transform(X_train)
	# train model
	selection_model = XGBRegressor(booster='gbtree',
                    eval_metric='rmse',
                    gamma = 0.4,
                    min_child_weight= 1,
                    max_depth= 9,
                    subsample= 0.9,
                    colsample_bytree= 0.2,
                    tree_method= 'exact',
                    learning_rate=0.1,
                    n_estimators=100,
                    nthread=10,
                    reg_alpha=1e-06,
                    seed=27)
	selection_model.fit(select_X_train, y_train)
	# eval model
	select_X_test = selection.transform(X_test)
	predictions = selection_model.predict(select_X_test)
	mse = mean_squared_error(y_test, predictions)
	print("Thresh=%.3f, n=%d, mse: %.2f" % (thresh, select_X_train.shape[1], mse)),thrs.append(thresh),x_train_s.append(select_X_train.shape[1]),scores.append(mse)

f_sel=pd.DataFrame([thrs,x_train_s,scores]).T
f_sel.columns=['thresh','feature_n','mse']
f_weight= pd.DataFrame([X_COL,model.feature_importances_]).T
f_weight.columns=['feature','weight']
f_weight_step1 = f_weight.sort_values(by='weight',ascending=True)
f_sel_step1= f_sel
f_weight_step1.to_csv("output/06_step1_feature_importance_3a_not_fill_pfs.txt",sep="\t",index=None)
f_sel_step1.to_csv("output/06_step1_feature_selection_3a_not_fill_pfs.txt",sep="\t",index=None)
#------------------------------------------------------------------------------------------------step2
