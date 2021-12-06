#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np

dataset = pd.read_table("../output/09_train_dataset.txt")
dataset.head


X = dataset.loc[:, ['B2mg','LN_num','LDH','age_raw','Lym_Mono','HGB','Ki.67','SPD','SUVmax','Bsym','BM']].values
y = dataset.loc[:, 'new_pod_total'].values

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
from sklearn.model_selection import GridSearchCV

train_X,test_X,train_Y,test_Y =train_test_split(X,y,test_size=0.3,random_state=23)

def modelfit(alg,useTrainCV=True, cv_folds=5, early_stopping_rounds=50):
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
    print(" Score (Train): %f" % metrics.accuracy_score(train_Y, dtrain_predictions))
    #Predict on testing data:
    dtest_predictions = alg.predict(test_X)
    print("Score (Test): %f" % metrics.accuracy_score(test_Y, dtest_predictions))


param_test1 = {
    'max_depth':range(3,10,2),
    'min_child_weight':range(1,6,2)}

gsearch1 = GridSearchCV(estimator = XGBClassifier(objective = "binary:logistic",eval_metric= "logloss",use_label_encoder=False,learning_rate = 0.1,n_estimators=50, seed=1024,importance_type='weight'), 
                       param_grid = param_test1, scoring='accuracy',cv=5)
gsearch1.fit(train_X,train_Y)
gsearch1.best_params_, gsearch1.best_score_
#({'max_depth': 3, 'min_child_weight': 5}, 0.8127659574468085)

param_test3 = {'gamma':[i/10.0 for i in range(0,20)]}
gsearch3 = GridSearchCV(estimator = XGBClassifier(max_depth =3, min_child_weight =5,objective = "binary:logistic",eval_metric= "logloss",use_label_encoder=False,learning_rate = 0.1,n_estimators=50, seed=1024,importance_type='weight'), 
                       param_grid = param_test3, scoring='accuracy',cv=5)
gsearch3.fit(train_X,train_Y)
gsearch3.best_params_, gsearch3.best_score_
#({'gamma': 1.8}, 0.825531914893617)
#-----
param_test4 = {
    'subsample':[i/10.0 for i in range(6,10)],
    'colsample_bytree':[i/10.0 for i in range(6,10)]
}
gsearch4 = GridSearchCV(estimator = XGBClassifier(gamma = 1.8, max_depth =3, min_child_weight =5,objective = "binary:logistic",eval_metric= "logloss",use_label_encoder=False,learning_rate = 0.1,n_estimators=50, seed=1024,importance_type='weight'), 
                       param_grid = param_test4, scoring='accuracy',cv=5)
gsearch4.fit(train_X,train_Y)
gsearch4.best_params_, gsearch4.best_score_
#({'colsample_bytree': 0.9, 'subsample': 0.9}, 0.8170212765957446)

param_test5 = {
    'subsample':[i/100.0 for i in range(80,100,5)],
    'colsample_bytree':[i/100.0 for i in range(80,100,5)]
}
gsearch5 = GridSearchCV(estimator = XGBClassifier(gamma = 1.8, max_depth =3, min_child_weight =5,objective = "binary:logistic",eval_metric= "logloss",use_label_encoder=False,learning_rate = 0.1,n_estimators=50, seed=1024,importance_type='weight'), 
                       param_grid = param_test5, scoring='accuracy',cv=5)
gsearch5.fit(train_X,train_Y)
gsearch5.best_params_, gsearch5.best_score_
#({'colsample_bytree': 0.85, 'subsample': 0.9}, 0.8170212765957446)
#------------
param_test6 = {
    'reg_alpha':[1e-5, 1e-2, 0.1, 1, 100]
}
gsearch6 = GridSearchCV(estimator = XGBClassifier(colsample_bytree=0.85, subsample =0.9, gamma = 1.8, max_depth =3, min_child_weight =5,objective = "binary:logistic",eval_metric= "logloss",use_label_encoder=False,learning_rate = 0.1,n_estimators=50, seed=1024,importance_type='weight'), 
                       param_grid = param_test6, scoring='accuracy',cv=5)
gsearch6.fit(train_X,train_Y)
gsearch6.best_params_, gsearch6.best_score_
#-----------------------
#{'reg_alpha': 1e-05}

param_test7 = {
    'n_estimators':[50, 100, 200, 500,1000],
    'learning_rate':[0.001, 0.01, 0.05, 0.1,0.2]
}
gsearch7 = GridSearchCV(estimator = XGBClassifier(reg_alpha = 1e-05,colsample_bytree=0.85, subsample =0.9, gamma = 1.8, max_depth =3, min_child_weight =5,objective = "binary:logistic",eval_metric= "logloss",use_label_encoder=False,learning_rate = 0.1,n_estimators=50, seed=1024,importance_type='weight'), 
                       param_grid = param_test7, scoring='accuracy',cv=5)
gsearch7.fit(train_X,train_Y)
gsearch7.best_params_, gsearch7.best_score_
#({'learning_rate': 0.001, 'n_estimators': 50}, 0.8212765957446809)




model = XGBClassifier(reg_alpha = 1e-05,
    colsample_bytree=0.85, 
    subsample =0.9, 
    gamma = 1.8, 
    max_depth =3, 
    min_child_weight =5,
    objective = "binary:logistic",
    eval_metric= "logloss",
    use_label_encoder=False,
    learning_rate = 0.001,
    n_estimators=50, 
    seed=1024,
    importance_type='weight')
