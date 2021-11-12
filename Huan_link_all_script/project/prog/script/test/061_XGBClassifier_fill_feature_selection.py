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
from sklearn.feature_selection import SelectFromModel
from numpy import sort
#----------------------------------------------------

dataset = pd.read_table("../output/01_add_age_raw_pfs_os_filter_grade_mergrPod_3A_fillna.txt")
# dataset = pd.read_table("../output/01_add_age_raw_pfs_os_filter_grade_mergrPod_3A.txt")

X = dataset.loc[:, ['gender','Ki.67','stage','Bsym','LN_num','site0','extend','BM','spleen','extend_num','BM_extend','LN6','SUVmax','SPD','X150b2mg_ldh','b2mg_LDH','ECOG','B2mg','LDH','LDH0','WBC','HGB','HGB0','PLT','Mono','Lym','interm_res','end_res','CR','Rmaintain','trans','age_raw','age_raw_60','ki67_20','LN_num_6','extend_num_0','stage_III_IV','stage_I_II','SUVmax_2','LDH_300','Lym_Mono_10','B2mg_3.4','SPD_0']]

y = dataset.loc[:, 'pod_total']

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
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=7)
# fit model on all training data
model.fit(X_train, y_train)
# make predictions for test data and evaluate
predictions = model.predict(X_test)
accuracy = accuracy_score(y_test, predictions)
print("Accuracy: %.2f%%" % (accuracy * 100.0))
# Fit model using each importance as a threshold
thrs=[]
x_train_s=[]
accs=[]
thresholds = sort(model.feature_importances_)
for thresh in thresholds:
	# select features using threshold
	selection = SelectFromModel(model, threshold=thresh, prefit=True)
	select_X_train = selection.transform(X_train)
	# train model
	selection_model = XGBClassifier()
	selection_model.fit(select_X_train, y_train)
	# eval model
	select_X_test = selection.transform(X_test)
	predictions = selection_model.predict(select_X_test)
	accuracy = accuracy_score(y_test, predictions)
	print("Thresh=%.3f, n=%d, Accuracy: %.2f%%" % (thresh, select_X_train.shape[1], accuracy*100.0)),thrs.append(thresh),x_train_s.append(select_X_train.shape[1]),accs.append(accuracy*100.0)

f_sel=pd.DataFrame([thrs,x_train_s,accs]).T
f_sel.columns=['thresh','feature_n','Accuracy']

f_weight= pd.DataFrame([X_train.columns,model.feature_importances_]).T
f_weight.columns=['feature','weight']
f_weight_step1 = f_weight.sort_values(by='weight',ascending=True)
f_sel_step1= f_sel
f_weight_step1.to_csv("../output/06_step1_feature_importance_3a_fill.txt",sep="\t",index=None)
f_sel_step1.to_csv("../output/06_step1_feature_selection_3a_fill.txt",sep="\t",index=None)
#-----------------------------------step2
X = dataset.loc[:, ['gender','Ki.67','stage','Bsym','LN_num','extend','BM','extend_num','BM_extend','LN6','SUVmax','SPD','b2mg_LDH','ECOG','B2mg','LDH','LDH0','WBC','HGB','PLT','Mono','Lym','interm_res','end_res','CR','Rmaintain','trans','age_raw','LN_num_6','extend_num_0']]
y = dataset.loc[:, 'pod_total']

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
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=7)
# fit model on all training data
model.fit(X_train, y_train)
# make predictions for test data and evaluate
predictions = model.predict(X_test)
accuracy = accuracy_score(y_test, predictions)
print("Accuracy: %.2f%%" % (accuracy * 100.0))
# Fit model using each importance as a threshold
thrs=[]
x_train_s=[]
accs=[]
thresholds = sort(model.feature_importances_)
for thresh in thresholds:
	# select features using threshold
	selection = SelectFromModel(model, threshold=thresh, prefit=True)
	select_X_train = selection.transform(X_train)
	# train model
	selection_model = XGBClassifier()
	selection_model.fit(select_X_train, y_train)
	# eval model
	select_X_test = selection.transform(X_test)
	predictions = selection_model.predict(select_X_test)
	accuracy = accuracy_score(y_test, predictions)
	print("Thresh=%.3f, n=%d, Accuracy: %.2f%%" % (thresh, select_X_train.shape[1], accuracy*100.0)),thrs.append(thresh),x_train_s.append(select_X_train.shape[1]),accs.append(accuracy*100.0)

f_sel=pd.DataFrame([thrs,x_train_s,accs]).T
f_sel.columns=['thresh','feature_n','Accuracy']

f_weight= pd.DataFrame([X_train.columns,model.feature_importances_]).T
f_weight.columns=['feature','weight']
f_weight_step2 = f_weight.sort_values(by='weight',ascending=True)
f_sel_step2= f_sel
#------------------------------------------------------------------------step3
X = dataset.loc[:, ['gender','Ki.67','stage','Bsym','LN_num','extend','BM','extend_num','BM_extend','LN6','SUVmax','SPD','b2mg_LDH','ECOG','B2mg','LDH','LDH0','WBC','HGB','PLT','Mono','Lym','interm_res','end_res','Rmaintain','trans','age_raw']]
y = dataset.loc[:, 'pod_total']

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
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=7)
# fit model on all training data
model.fit(X_train, y_train)
# make predictions for test data and evaluate
predictions = model.predict(X_test)
accuracy = accuracy_score(y_test, predictions)
print("Accuracy: %.2f%%" % (accuracy * 100.0))
# Fit model using each importance as a threshold
thrs=[]
x_train_s=[]
accs=[]
thresholds = sort(model.feature_importances_)
for thresh in thresholds:
	# select features using threshold
	selection = SelectFromModel(model, threshold=thresh, prefit=True)
	select_X_train = selection.transform(X_train)
	# train model
	selection_model = XGBClassifier()
	selection_model.fit(select_X_train, y_train)
	# eval model
	select_X_test = selection.transform(X_test)
	predictions = selection_model.predict(select_X_test)
	accuracy = accuracy_score(y_test, predictions)
	print("Thresh=%.3f, n=%d, Accuracy: %.2f%%" % (thresh, select_X_train.shape[1], accuracy*100.0)),thrs.append(thresh),x_train_s.append(select_X_train.shape[1]),accs.append(accuracy*100.0)

f_sel=pd.DataFrame([thrs,x_train_s,accs]).T
f_sel.columns=['thresh','feature_n','Accuracy']

f_weight= pd.DataFrame([X_train.columns,model.feature_importances_]).T
f_weight.columns=['feature','weight']
f_weight_step3 = f_weight.sort_values(by='weight',ascending=True)
f_sel_step3= f_sel
f_weight_step3.to_csv("../output/06_step3_feature_importance_3a_fill.txt",sep="\t",index=None)
f_sel_step3.to_csv("../output/06_step3_feature_selection_3a_fill.txt",sep="\t",index=None)
#--------------------------------------step4
