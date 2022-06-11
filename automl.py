# %% [markdown]
#  # Radiomics model building using Auto ML
#  Wookjin Choi <wookjin.choi@jefferson.edu>

# %% [markdown]
#  ## Import modules

# %%
import autosklearn
import autosklearn.classification
import sklearn.model_selection
import sklearn.datasets
import sklearn.metrics
from sklearn import datasets, cluster
from sklearn.tree import DecisionTreeClassifier
from sklearn.feature_selection import GenericUnivariateSelect
from sklearn.feature_selection import RFE
from sklearn.feature_selection import RFECV
from sklearn.feature_selection import mutual_info_classif
from sklearn.feature_selection import SelectFromModel
import pandas as pd
import numpy as np


# %% [markdown]
#  ## 1. Data loading

# %%
clinical_lidc72 = pd.read_csv("metadata/LIDC_72.csv")
clinical_lidc = pd.read_csv("metadata/LIDC.csv")
clinical_lungx = pd.read_csv("metadata/LUNGx.csv")

#features_lidc = pd.read_csv("output/feature-list_nodule-lidc_pylidc.csv")
features_lidc = pd.read_csv("output/feature-list_nodule-lidc.csv")
features_lungx = pd.read_csv("output/feature-list_nodule-lungx.csv")

lidc_tags = features_lidc.Tags.str.split("-").apply(pd.Series)
features_lidc['NID'] = lidc_tags[0].astype(int)
features_lidc['Tags'] = lidc_tags[1]

lungx_tags = features_lungx.Tags.str.split("-").apply(pd.Series)
lungx_tags.iloc[lungx_tags[2].isna(), 1] = "gc"
#lungx_tags.iloc[lungx_tags[2].isna(), 2] = "1mm"

features_lungx['NID'] = lungx_tags[0].astype(int)
features_lungx['Tags'] = lungx_tags[1]

features_lidc = features_lidc.loc[(features_lidc.Tags=='seg') & (features_lidc.NID==1)]
features_lungx = features_lungx.loc[(features_lungx.Tags=='seg')]

#features_lidc = features_lidc.drop(columns=["WeightedElongation","WeightedFlatness","2DWeightedElongation","2DWeightedFlatness"])
#features_lungx = features_lungx.drop(columns=["WeightedElongation","WeightedFlatness","2DWeightedElongation","2DWeightedFlatness"])
features_lidc72 = features_lidc.loc[features_lidc.PID.isin(clinical_lidc72.PID)]
clinical_lidc = clinical_lidc.loc[clinical_lidc.set_index(['PID','NID']).index.isin(features_lidc.set_index(['PID','NID']).index)]
features_lidc = features_lidc.loc[features_lidc.set_index(['PID','NID']).index.isin(clinical_lidc.set_index(['PID','NID']).index)]
features_lidc_train = features_lidc.loc[~features_lidc.PID.isin(clinical_lidc72.PID)]
clinical_lidc_train = clinical_lidc.loc[~clinical_lidc.PID.isin(clinical_lidc72.PID)]

# %%
spiculation_lidc = pd.read_csv("output/nodule-lidc/spiculation/features.csv")
spiculation_lungx = pd.read_csv("output/nodule-lungx/spiculation/features.csv")

spiculation_lidc = spiculation_lidc.loc[spiculation_lidc.NID==1]

spiculation_lidc72 = spiculation_lidc.loc[spiculation_lidc.PID.isin(clinical_lidc72.PID)]
spiculation_lidc = spiculation_lidc.loc[spiculation_lidc.set_index(['PID','NID']).index.isin(clinical_lidc.set_index(['PID','NID']).index)]
spiculation_lidc_train = spiculation_lidc.loc[~spiculation_lidc.PID.isin(clinical_lidc72.PID)]

# %% [markdown]
#  ### 1.1 Data preprocessing

# %%
#X = features_lidc_train.iloc[:, 4:-1].update(spiculation_lungx.iloc[:,3:])
X = pd.concat([features_lidc_train.reset_index().iloc[:, 5:-1], spiculation_lidc_train.reset_index().iloc[:,4:]], axis=1) # with spiculation features
X = X.loc[:,X.isna().sum(axis=0)==0]
y = clinical_lidc_train.Malignancy > 3

#X_72 = features_lidc72.iloc[:, 4:-1]
X_72 = pd.concat([features_lidc72.reset_index().iloc[:, 5:-1], spiculation_lidc72.reset_index().iloc[:,4:]], axis=1) # with spiculation features
X_72 = X_72.loc[:, X.columns]
y_72 = clinical_lidc72.PMalignancy == 2

#X_lungx = features_lungx.iloc[:, 4:-1]
X_lungx = pd.concat([features_lungx.reset_index().iloc[:, 5:-1], spiculation_lungx.reset_index().iloc[:,4:]], axis=1) # with spiculation features
X_lungx = X_lungx.loc[:, X.columns]
y_lungx = clinical_lungx.malignancy > 0


X_lungx_cal = X_lungx.iloc[0:10]
y_lungx_cal = y_lungx.iloc[0:10]
X_lungx_test = X_lungx.iloc[10:]
y_lungx_test = y_lungx.iloc[10:]


# %% [markdown]
#  ### 1.2 Data splitting

# %%
X_train, X_test, y_train, y_test = \
    sklearn.model_selection.train_test_split(X, y, test_size=0.2, random_state=42)


# %% [markdown]
#  ## 2. Feature Selection

# %% [markdown]
#  ### 2.1 Redundant Feature Elimination

# %%
agglo = cluster.FeatureAgglomeration(affinity='precomputed', linkage='complete', n_clusters=None, distance_threshold=0.3)
agglo.fit(1-abs(X_train.corr()))
cv_f = X_train.std()/X_train.mean().abs()

feature_groups = [np.array(X_train.columns[agglo.labels_==i]) for i in np.unique(agglo.labels_)]
feature_reps = [cv_f[f].sort_values().index[-1] for f in feature_groups]
feature_reps


# %%
#X_reduced_tr = pd.DataFrame(agglo.transform(X_train))
#X_reduced_ts = pd.DataFrame(agglo.transform(X_test))
#X_reduced_tr.columns = feature_reps
#X_reduced_ts.columns = feature_reps

X_reduced_tr = X_train[feature_reps]
X_reduced_ts = X_test[feature_reps]

print(X_reduced_tr.shape, X_reduced_ts.shape)
feature_groups[0].shape


# %%
X_train = X_reduced_tr
X_test = X_reduced_ts


# %% [markdown]
#  ### 2.2 Feature Selection

# %%
trans = GenericUnivariateSelect(score_func=lambda X, y: X.mean(axis=0), mode='percentile', param=50)
X_trans = trans.fit_transform(X_train, y_train)

mutual_information = mutual_info_classif(X_train, y_train) 
trans = GenericUnivariateSelect(score_func=mutual_info_classif, mode='k_best', param=30)
X_trans = trans.fit_transform(X_train, y_train)
columns_retained_Select = X_train.columns[trans.get_support()].values
columns_retained_Select


# %%
clf = DecisionTreeClassifier()
clf.fit(X, y)

trans = SelectFromModel(clf, threshold='median', max_features=30)
X_trans = trans.fit_transform(X, y)
columns_retained_FromMode = X.columns[trans.get_support()].values
columns_retained_FromMode


# %%
clf = DecisionTreeClassifier()
trans = RFE(clf, n_features_to_select=30)
X_trans = trans.fit_transform(X_train, y_train)
columns_retained_RFE = X_train.columns[trans.get_support()].values
columns_retained_RFE


# %%
clf = DecisionTreeClassifier()
trans = RFECV(clf)
X_trans = trans.fit_transform(X_train, y_train)
columns_retained_RFECV = X_train.columns[trans.get_support()].values
columns_retained_RFECV


# %%
columns_retained = columns_retained_RFE

# %% [markdown]
#  ## 3. AutoML

# %% [markdown]
#  ### 3.1 TPOT

# %%
from tpot import TPOTClassifier

# X_train, X_test, y_train, y_test = \
#    sklearn.model_selection.train_test_split(X, y, random_state=42)
    
pipeline_optimizer = TPOTClassifier(generations=10, population_size=500, cv=5,
                                    random_state=42, verbosity=2)
# pipeline_optimizer.fit(X_train, y_train)
# print(pipeline_optimizer.score(X_test, y_test))
pipeline_optimizer.fit(X_train[columns_retained], y_train)
#pipeline_optimizer.export('tpot_exported_pipeline.py')


# %% [markdown]
#  #### 3.2.1 TPOT validation

# %%
print("Test AUC:", pipeline_optimizer.score(X_test[columns_retained], y_test))
print("72 AUC:", pipeline_optimizer.score(X_72[columns_retained], y_72))
print("LUNGx AUC:", pipeline_optimizer.score(X_lungx_test[columns_retained], y_lungx_test))


# %% [markdown]
#  ### 3.2 AutoSklearn

# %%
automl = autosklearn.classification.AutoSklearnClassifier(
    time_left_for_this_task = 8000,
    per_run_time_limit = 800,
    seed=42, 
    memory_limit=None,
    metric=autosklearn.metrics.roc_auc,
    resampling_strategy='cv',
    resampling_strategy_arguments=dict(folds=5),
)
automl.fit(X_train[columns_retained], y_train)

print(np.max(automl.cv_results_['mean_test_score']))


# %% [markdown]
#  #### 3.2.1 AutoSklearn validation

# %%
#predict
y_hat = automl.predict(X_test[columns_retained])
y_pred = automl.predict_proba(X_test[columns_retained])

y_hat_72 = automl.predict(X_72[columns_retained])
y_pred_72 = automl.predict_proba(X_72[columns_retained])

y_hat_lungx = automl.predict(X_lungx_test[columns_retained])
y_pred_lungx = automl.predict_proba(X_lungx_test[columns_retained])

# show scores
print("Test AUC: ", sklearn.metrics.roc_auc_score(y_test, y_pred[:,1]))
print("72 AUC: ", sklearn.metrics.roc_auc_score(y_72, y_pred_72[:,1]))
print("LUNGx AUC: ", sklearn.metrics.roc_auc_score(y_lungx_test, y_pred_lungx[:,1]))

print("\nTest Accuracy: ", sklearn.metrics.accuracy_score(y_test, y_hat))
print("72 Accuracy: ", sklearn.metrics.accuracy_score(y_72, y_hat_72))
print("LUNGx Accuracy: ", sklearn.metrics.accuracy_score(y_lungx_test, y_hat_lungx))


# %% [markdown]
#  #### 3.2.2 AutoSklearn models

# %%
# show all models
show_models_str=automl.show_models()
sprint_statistics_str = automl.sprint_statistics()
leaderboard = automl.leaderboard()

print(sprint_statistics_str)
print(leaderboard)

# %%
print(show_models_str[leaderboard.index[0]])

# %%
import matplotlib.pyplot as plt

poT = automl.performance_over_time_
poT.plot(
    x='Timestamp',
    kind='line',
    legend=True,
    title='Auto-sklearn AUC over time',
    grid=True,
)
plt.show()


# %%
def get_metric_result(cv_results):
    results = pd.DataFrame.from_dict(cv_results)
    results = results[results['status'] == "Success"]
    cols = ['rank_test_scores', 'param_classifier:__choice__', 'mean_test_score']
    cols.extend([key for key in cv_results.keys() if key.startswith('metric_')])
    return results[cols]

print("Metric results")
print(get_metric_result(automl.cv_results_).to_string(index=False))
