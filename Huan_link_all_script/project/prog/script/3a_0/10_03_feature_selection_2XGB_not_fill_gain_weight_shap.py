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
dataset = pd.read_table("./output/09_train_dataset.txt")
# covariates <- c('Ki.67','stage','Bsym','LN_num','LN6','BM','spleen','extend_num','BM_extend','SUVmax','SPD','ECOG','B2MG_re0_train','LDH_re0_train','HGB','age_raw','Lym_Mono')
X = dataset.loc[:, ["LDH","LN_num","HGB","B2mg","ECOG","BM_extend","SPD","Lym_Mono","LN6","BM","SUVmax"]]
y = dataset.loc[:, 'new_pod_total']
X_COL = X.columns
# X= dataset.iloc[:,1:540]
# y = dataset.loc[:, 'AUC'].values
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=7)

#-----------------------Accuracy
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

model.fit(X_train,y_train)
predictions = model.predict(X_test)
accuracy = accuracy_score(y_test, predictions)
print("Accuracy: %.2f%%" % (accuracy * 100.0))
#79.87%

cv = KFold(n_splits=10, shuffle=True, random_state=0)
y_pred = cross_val_predict(model, X, y,cv=cv)
accuracy = accuracy_score(y, y_pred)
print("Accuracy: %.2f%%" % (accuracy * 100.0))
#80.93%
cv = KFold(n_splits=5, shuffle=True, random_state=10)
y_pred = cross_val_predict(model, X, y,cv=cv)
accuracy = accuracy_score(y, y_pred)
print("Accuracy: %.2f%%" % (accuracy * 100.0))
#81.43%
#--------------------------------------------------importance
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


model.fit(X,y)

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
im_gain_df.to_csv("./output/10_03_not_fill_feature_importance_gain.txt",sep="\t",index=None)
#--------------------------------------
im_weight_df = pd.DataFrame([importance_weight]).T
im_weight_df.columns=['Feature_importance']
im_weight_df['Feature']=im_weight_df.index
im_weight_df=im_weight_df[order].sort_values(by='Feature_importance',ascending=False)
im_weight_df.to_csv("./output/10_03_not_fill_feature_importance_weight.txt",sep="\t",index=None)


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
feature_shap.to_csv("./output/10_03_not_fill_feature_importance_shap_class.txt",sep="\t",index=None)