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
import shap
#----------------------------------------------------

# dataset = pd.read_table("../output/01_add_age_raw_pfs_os_filter_grade_mergrPod_3A_fillna.txt")
dataset = pd.read_table("../output/01_add_age_raw_pfs_os_filter_grade_mergrPod_3A.txt")

X = dataset.loc[:, ['gender','Ki.67','stage','Bsym','LN_num','site0','extend','BM','spleen','extend_num','BM_extend','LN6','SUVmax','SPD','X150b2mg_ldh','b2mg_LDH','ECOG','B2mg','LDH','LDH0','HGB','HGB0','Mono','Lym','age_raw','age_raw_60','ki67_20','LN_num_6','extend_num_0','SUVmax_2','LDH_300','Lym_Mono','B2mg_3.4','SPD_0']]
y = dataset.loc[:, 'pod_total']
X_COL = X.columns
# X= dataset.iloc[:,1:540]
# y = dataset.loc[:, 'AUC'].values
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=7)

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

# model.fit(X_train,y_train)




# # fig, ax = plt.subplots(1, 1, figsize=(8, 13))
# # xgb.plot_importance(model, max_num_features=30, height=0.5, ax=ax)
# # plt.savefig('../output/figure/06_feature_importance_3a_not_fill.pdf')

# # make predictions for test data and evaluate
# predictions = model.predict(X_test)
# accuracy = accuracy_score(y_test, predictions)
# print("Accuracy: %.2f%%" % (accuracy * 100.0))
# #76.32%
#-----------------------
model = XGBClassifier(learning_rate=0.001, 
    n_estimators = 50,
    reg_alpha= 1e-5,
    colsample_bytree = 0.5, 
    subsample= 0.65,
    gamma= 2.4, 
    min_child_weight =5,
    max_depth= 3, 
    objective = "multi:softprob",
    eval_metric= "mlogloss",
    use_label_encoder=False,
    seed=1024,
    importance_type='weight')
model.fit(X_train,y_train)
predictions = model.predict(X_test)
accuracy = accuracy_score(y_test, predictions)
print("Accuracy: %.2f%%" % (accuracy * 100.0))
#78.95%
#------------
booster = model.get_booster()
importance_weight = booster.get_score(importance_type="weight")
im_weight_df = pd.DataFrame([importance_weight]).T
im_weight_df.columns=['Feature_importance']
im_weight_df['Feature']=im_weight_df.index
im_weight_df_step1 = im_weight_df.sort_values(by='Feature_importance',ascending=False)
step1_eff_features = im_weight_df_step1.loc[:,'Feature'].values.tolist()
#------------------------------------------------------------------------------step2
ori_colnums =X.columns.values.tolist()
list(set(ori_colnums).difference(set(step1_eff_features)))

X = dataset.loc[:, ['gender','Ki.67','Bsym','LN_num','BM','spleen','LN6','SUVmax','SPD','b2mg_LDH','B2mg','LDH','LDH0','HGB','HGB0','Mono','Lym','age_raw','LN_num_6','LDH_300','Lym_Mono','B2mg_3.4']]
y = dataset.loc[:, 'pod_total']
X_COL = X.columns
# X= dataset.iloc[:,1:540]
# y = dataset.loc[:, 'AUC'].values
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=7)

model=XGBClassifier(reg_alpha = 1e-05, 
    colsample_bytree= 0.65, 
    subsample= 0.5,
    gamma = 1.7,
    min_child_weight =3,
    max_depth= 3, 
    objective = "multi:softprob",
    eval_metric= "mlogloss",
    use_label_encoder=False,
    learning_rate = 0.1,
    n_estimators=50, 
    importance_type="weight",
    seed=1024)

model.fit(X_train,y_train)
predictions = model.predict(X_test)
accuracy = accuracy_score(y_test, predictions)
print("Accuracy: %.2f%%" % (accuracy * 100.0))
#73.68%
booster = model.get_booster()
importance_weight = booster.get_score(importance_type="weight")
im_weight_df = pd.DataFrame([importance_weight]).T
im_weight_df.columns=['Feature_importance']
im_weight_df['Feature']=im_weight_df.index
im_weight_df_step2 = im_weight_df.sort_values(by='Feature_importance',ascending=False)
step2_eff_features = im_weight_df_step2.loc[:,'Feature'].values.tolist()

#------------------------------------------------------------------------step3
list(set(step1_eff_features).difference(set(step2_eff_features)))

X = dataset.loc[:, ['Ki.67','Bsym','LN_num','BM','spleen','LN6','SUVmax','SPD','b2mg_LDH','B2mg','LDH','LDH0','HGB','HGB0','Mono','Lym','age_raw','LN_num_6','Lym_Mono','B2mg_3.4']]
y = dataset.loc[:, 'pod_total']
X_COL = X.columns
# X= dataset.iloc[:,1:540]
# y = dataset.loc[:, 'AUC'].values
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=7)

model=XGBClassifier(reg_alpha = 1e-05, 
    colsample_bytree= 0.65, 
    subsample= 0.5,
    gamma = 1.7,
    min_child_weight =3,
    max_depth= 3, 
    objective = "multi:softprob",
    eval_metric= "mlogloss",
    use_label_encoder=False,
    learning_rate = 0.1,
    n_estimators=50, 
    importance_type="weight",
    seed=1024)

model.fit(X_train,y_train)
predictions = model.predict(X_test)
accuracy = accuracy_score(y_test, predictions)
print("Accuracy: %.2f%%" % (accuracy * 100.0))
#76.97%
booster = model.get_booster()
importance_weight = booster.get_score(importance_type="weight")
im_weight_df = pd.DataFrame([importance_weight]).T
im_weight_df.columns=['Feature_importance']
im_weight_df['Feature']=im_weight_df.index
im_weight_df_step3 = im_weight_df.sort_values(by='Feature_importance',ascending=False)
step3_eff_features = im_weight_df_step3.loc[:,'Feature'].values.tolist()
#------------------------------------------------------------------------------step4
list(set(step2_eff_features).difference(set(step3_eff_features)))
# #---------------------------------------------------------feature importance
booster = model.get_booster()
importance_gain = booster.get_score(importance_type="gain")
importance_weight = booster.get_score(importance_type="weight")
# im_gain_df = pd.DataFrame(importance_gain,index=[0]).T
im_gain_df = pd.DataFrame([importance_gain]).T
# im_gain_df = pd.DataFrame.from_dict(importance_gain,orient="index")
im_gain_df.columns=['Feature_importance']
im_gain_df['Feature']=im_gain_df.index
order=['Feature','Feature_importance']
im_gain_df=im_gain_df[order].sort_values(by='Feature_importance',ascending=False)
im_gain_df.to_csv("../output/064_step3_not_fill_weight_feature_importance_gain.txt",sep="\t",index=None)
#--------------------------------------
im_weight_df = pd.DataFrame([importance_weight]).T
im_weight_df.columns=['Feature_importance']
im_weight_df['Feature']=im_weight_df.index
im_weight_df=im_weight_df[order].sort_values(by='Feature_importance',ascending=False)
im_weight_df.to_csv("../output/064_step3_not_fill_weight_feature_importance_weight.txt",sep="\t",index=None)


background = shap.maskers.Independent(X)
def f(x):
    return shap.links.identity(model.predict_proba(x, validate_features=False)[:,1])
explainer = shap.Explainer(f, background, link=shap.links.logit)
shap_values = explainer(X)

# visualize the first prediction's explanation
# shap.plots.bar(shap_values)
# plt.savefig('../output/figure/064_step3_fill_feature_shap_bar.pdf')
# plt.show()


feature_shap =pd.DataFrame(shap_values.values)
feature_shap.columns =X.columns
feature_shap.to_csv("../output/064_step3_weight_not_fill_feature_importance_shap_class.txt",sep="\t",index=None)