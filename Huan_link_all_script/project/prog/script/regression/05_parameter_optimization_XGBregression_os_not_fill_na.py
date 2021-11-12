#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np

dataset = pd.read_table("../../output/01_add_age_raw_pfs_os_filter_grade_mergrPod_3A.txt")
# dataset = pd.read_table("../output/01_add_age_raw_pfs_os_filter_grade_mergrPod_3A_fillna.txt")
dataset.head
# y = dataset.iloc[:, -1].values

X = dataset.loc[:, ['gender','Ki.67','stage','Bsym','LN_num','site0','extend','BM','spleen','extend_num','BM_extend','LN6','SUVmax','SPD','X150b2mg_ldh','b2mg_LDH','ECOG','B2mg','LDH','LDH0','WBC','HGB','HGB0','PLT','Mono','Lym','age_raw','age_raw_60','ki67_20','LN_num_6','extend_num_0','SUVmax_2','LDH_300','Lym_Mono_10','B2mg_3.4','SPD_0']].values

y = dataset.loc[:, 'OS_month'].values

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from sklearn.metrics import plot_roc_curve
from scipy.stats import spearmanr
from pandas import read_csv
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_validate, cross_val_predict
from sklearn import metrics
from xgboost import XGBRegressor
import xgboost as xgb
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from sklearn.model_selection import GridSearchCV
# %matplotlib inline
from matplotlib.pylab import rcParams
rcParams['figure.figsize'] = 12, 4
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import fbeta_score, make_scorer
from sklearn.metrics import mean_squared_error
from math import sqrt
#---------------------
train_X,test_X,train_Y,test_Y =train_test_split(X,y,test_size=0.3,random_state=23)
def modelfit(alg,useTrainCV=True, cv_folds=10, early_stopping_rounds=50):
    if useTrainCV:
        xgb_param = alg.get_xgb_params()
        xgtrain =xgb.DMatrix(train_X,label=train_Y)
        xgtest = xgb.DMatrix(test_X)
        cvresult = xgb.cv(xgb_param, xgtrain, num_boost_round=alg.get_params()['n_estimators'], nfold=cv_folds,
                          early_stopping_rounds=early_stopping_rounds,show_stdv=False)
        alg.set_params(n_estimators=cvresult.shape[0])#cvresult.shape[0]和alg.get_params()['n_estimators']值一样
    #Fit the algorithm on the data
    alg.fit(train_X, train_Y,eval_metric='rmse')
    #Predict training set:
    dtrain_predictions = alg.predict(train_X)
    #Print model report:
    print(" Score (Train): %f" % metrics.mean_squared_error(train_Y, dtrain_predictions))
    #Predict on testing data:
    dtest_predictions = alg.predict(test_X)
    print("Score (Test): %f" % metrics.mean_squared_error(test_Y, dtest_predictions))

xgb1 = XGBRegressor(booster='gbtree',
                    eval_metric='rmse',
                    gamma = 0.1,
                    min_child_weight= 1.1,
                    max_depth= 5,
                    subsample= 0.8,
                    colsample_bytree= 0.8,
                    tree_method= 'exact',
                    learning_rate=0.1,
                    n_estimators=100,
                    nthread=10,
                    seed=27)
modelfit(xgb1)


#max_depth 和 min_child_weight
param_test1 = {
    'max_depth':[3,5,7,9],
    'min_child_weight':[1,3,5]
}
gsearch1 = GridSearchCV(estimator = XGBRegressor(booster='gbtree',
                    eval_metric='rmse',
                    gamma = 0.1,
                    min_child_weight= 1.1,
                    max_depth= 5,
                    subsample= 0.8,
                    colsample_bytree= 0.8,
                    tree_method= 'exact',
                    learning_rate=0.1,
                    n_estimators=100,
                    nthread=10,
                    seed=27),
                       param_grid = param_test1, scoring='neg_mean_squared_error',n_jobs=10, cv=10)
gsearch1.fit(train_X,train_Y)
# gsearch1.grid_scores_,
gsearch1.best_params_, gsearch1.best_score_
#{'max_depth': 5, 'min_child_weight': 1}

#-----------------------------gamma
%%time
param_test3 = {
    'gamma':[i/10.0 for i in range(0,5)]
}
gsearch3 = GridSearchCV(estimator = XGBRegressor(booster='gbtree',
                    eval_metric='rmse',
                    gamma = 0.1,
                    min_child_weight= 1,
                    max_depth= 5,
                    subsample= 0.8,
                    colsample_bytree= 0.8,
                    tree_method= 'exact',
                    learning_rate=0.1,
                    n_estimators=100,
                    nthread=10,
                    seed=27),
                       param_grid = param_test3, scoring='neg_mean_squared_error',n_jobs=10, cv=10)
gsearch3.fit(train_X,train_Y)
gsearch3.best_params_, gsearch3.best_score_
#{'gamma': 0.1}
#------------------subsample 和 colsample_bytree 
%%time
#Grid seach on subsample and max_features
#Choose all predictors except target & IDcols
param_test4 = {
    'subsample':[i/10.0 for i in range(6,10)],
    'colsample_bytree':[i/10.0 for i in range(6,10)]
}
gsearch4 = GridSearchCV(estimator = XGBRegressor(booster='gbtree',
                    eval_metric='rmse',
                    gamma = 0.1,
                    min_child_weight= 1,
                    max_depth= 5,
                    subsample= 0.8,
                    colsample_bytree= 0.8,
                    tree_method= 'exact',
                    learning_rate=0.1,
                    n_estimators=100,
                    nthread=10,
                    seed=27),
                       param_grid = param_test4, scoring='neg_mean_squared_error',n_jobs=10, cv=10)
gsearch4.fit(train_X,train_Y)
gsearch4.best_params_, gsearch4.best_score_
#({'colsample_bytree': 0.8, 'subsample': 0.8}
#--------------------reg_alpha、reg_lambda
#粗调
%%time
#Grid seach on subsample and max_features
#Choose all predictors except target & IDcols
param_test6 = {
    'reg_alpha':[1e-5, 1e-2, 0.1, 1, 100]
}
gsearch6 = GridSearchCV(estimator = XGBRegressor(booster='gbtree',
                    eval_metric='rmse',
                    gamma = 0.1,
                    min_child_weight= 1,
                    max_depth= 5,
                    subsample= 0.8,
                    colsample_bytree= 0.8,
                    tree_method= 'exact',
                    learning_rate=0.1,
                    n_estimators=100,
                    nthread=10,
                    seed=27),
                       param_grid = param_test6, scoring='neg_mean_squared_error',n_jobs=10, cv=10)
gsearch6.fit(train_X,train_Y)
gsearch6.best_params_, gsearch6.best_score_
#{'reg_alpha': 1e-05},
#微调
%%time
#Grid seach on subsample and max_features
#Choose all predictors except target & IDcols
param_test7 = {
    'reg_alpha':[1e-6,1e-5,1e-4] #----------修改
}
gsearch7 = GridSearchCV(estimator = XGBRegressor(booster='gbtree',
                    eval_metric='rmse',
                    gamma = 0.1,
                    min_child_weight= 1,
                    max_depth= 5,
                    subsample= 0.8,
                    colsample_bytree= 0.8,
                    tree_method= 'exact',
                    learning_rate=0.1,
                    n_estimators=100,
                    nthread=10,
                    seed=27),
                       param_grid = param_test7, scoring='neg_mean_squared_error',n_jobs=10, cv=10)
gsearch7.fit(train_X,train_Y)
gsearch7.best_params_, gsearch7.best_score_
#{'reg_alpha': 1e-05},
param_test8 = {
    'reg_alpha':[9e-6,1e-5,2e-5,3e-5] #----------修改
}
gsearch8 = GridSearchCV(estimator = XGBRegressor(booster='gbtree',
                    eval_metric='rmse',
                    gamma = 0.1,
                    min_child_weight= 1,
                    max_depth= 5,
                    subsample= 0.8,
                    colsample_bytree= 0.8,
                    tree_method= 'exact',
                    learning_rate=0.1,
                    n_estimators=100,
                    nthread=10,
                    seed=27),
                       param_grid = param_test8, scoring='neg_mean_squared_error',n_jobs=10, cv=10)
gsearch8.fit(train_X,train_Y)
gsearch8.best_params_, gsearch8.best_score_
#{'reg_alpha': 2e-05}
#--------------------learning_rate、n_estimators
param_test9 = {
    'n_estimators':[50, 100, 200, 500,1000],
    'learning_rate':[0.001, 0.01, 0.05, 0.1,0.2]
}
gsearch9 = GridSearchCV(estimator = XGBRegressor(booster='gbtree',
                    eval_metric='rmse',
                    gamma = 0.1,
                    min_child_weight= 1,
                    max_depth= 5,
                    subsample= 0.8,
                    colsample_bytree= 0.8,
                    tree_method= 'exact',
                    learning_rate=0.1,
                    n_estimators=100,
                    nthread=10,
                    reg_alpha = 2e-05,
                    seed=27),
                       param_grid = param_test9, scoring='neg_mean_squared_error',n_jobs=10, cv=10)
gsearch9.fit(train_X,train_Y)
gsearch9.best_params_, gsearch9.best_score_
#{'learning_rate': 0.05, 'n_estimators': 100}




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


model1=XGBRegressor()
model1.fit(train_X,train_Y)
y_predict=model1.predict(test_X)
mean_squared_error(test_Y,y_predict)
#1127.9463454265656
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

model.fit(train_X,train_Y)
y_predict=model.predict(test_X)
mean_squared_error(test_Y,y_predict)
#953.0560365109402
