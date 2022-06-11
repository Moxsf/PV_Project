from pytorch_tabnet.tab_model import TabNetClassifier
from pytorch_tabnet.augmentations import ClassificationSMOTE
import torch
import optuna

from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn import metrics

import pandas as pd
import numpy as np
import pickle5 as pickle

from pathlib import Path
import os
import shutil
import gzip

np.random.seed(0)

## base INFO
ModelF_path = "/beegfs/home/fanxutong/Project_RV/21_nonSigGWAS/02_Pre2Model/04_FModel/05_MoData_NewF/"
ModelO_path = "/beegfs/home/fanxutong/Project_RV/21_nonSigGWAS/06_TabNetC/"
TestF_path = "/beegfs/home/fanxutong/Project_RV/21_nonSigGWAS/03_TestD/03_FTest/05_TData_NewF/"

## Model Out
target = "CLASS"
dataset_name = "DNPGO_model"
#Model_out = ModelO_path + dataset_name

## Model Data IN
filenm = "DNPGO_FModel_SNVs.tsv"
file_in = ModelF_path + filenm

train = pd.read_csv(file_in, header = 0, sep = "\t")
train = train.sample(frac=1).reset_index(drop=True)
n_total = len(train)

train_indices, valid_indices = train_test_split(range(n_total), test_size=0.3, random_state=0)

## test data
RareP_FI = TestF_path + "RareP_FTest_SNVS_order.tsv"
RV300_FI = TestF_path + "RV300_FTest_SNVS_order.tsv"
RV555_FI = TestF_path + "RV555_FTest_SNVS_order.tsv"
RV800_FI = TestF_path + "RV800_FTest_SNVS_order.tsv"

RarePT = pd.read_csv(RareP_FI, header = 0, sep = "\t")
RV300T = pd.read_csv(RV300_FI, header = 0, sep = "\t")
RV555T = pd.read_csv(RV555_FI, header = 0, sep = "\t")
RV800T = pd.read_csv(RV800_FI, header = 0, sep = "\t")

categorical_columns = []
categorical_dims =  {}

for col in train.columns[train.dtypes == object]:
    print(col, train[col].nunique())
    l_enc = LabelEncoder()
    train[col] = train[col].fillna("VV_likely")
    train[col] = l_enc.fit_transform(train[col].values)
    categorical_columns.append(col)
    categorical_dims[col] = len(l_enc.classes_)

for col in RarePT.columns[RarePT.dtypes == object]:
    l_enc = LabelEncoder()
    RarePT[col] = RarePT[col].fillna("VV_likely")
    RarePT[col] = l_enc.fit_transform(RarePT[col].values)
    RV300T[col] = RV300T[col].fillna("VV_likely")
    RV300T[col] = l_enc.fit_transform(RV300T[col].values)
    RV555T[col] = RV555T[col].fillna("VV_likely")
    RV555T[col] = l_enc.fit_transform(RV555T[col].values)
    RV800T[col] = RV800T[col].fillna("VV_likely")
    RV800T[col] = l_enc.fit_transform(RV800T[col].values)

for col in train.columns[train.dtypes == 'float64']:
    train.fillna(train.loc[train_indices, col].mean(), inplace=True)

unused_feat = []
features = [ col for col in train.columns if col not in unused_feat + [target]] 
cat_idxs = [ i for i, f in enumerate(features) if f in categorical_columns]
cat_dims = [ categorical_dims[f] for i, f in enumerate(features) if f in categorical_columns]

if os.getenv("CI", False):
# Take only a subsample to run CI
    X_train = train[features].values[train_indices][:1000,:]
    y_train = train[target].values[train_indices][:1000]
else:
    X_train = train[features].values[train_indices]
    y_train = train[target].values[train_indices]

X_valid = train[features].values[valid_indices]
y_valid = train[target].values[valid_indices]

## test process
X_RarePT = RarePT[features].values
y_RarePT = RarePT[target].values

X_RV300T = RV300T[features].values
y_RV300T = RV300T[target].values

X_RV555T = RV555T[features].values
y_RV555T = RV555T[target].values

X_RV800T = RV800T[features].values
y_RV800T = RV800T[target].values

def metics_cal(y_RarePT, y_prop_RarePT, Datanm):
    AUPR = dict()
    AUROC = dict()
    CLASS = dict()

    AUCi = 1
    Monm = "DNPGO"
    AUROC[AUCi] = metrics.roc_auc_score(y_RarePT, y_prop_RarePT[1])
    AUPR[AUCi] = metrics.average_precision_score(y_RarePT, y_prop_RarePT[1])
    CLASS[AUCi] = Monm + "_1"
    
    if (2 in y_prop_RarePT.columns) or (3 in y_prop_RarePT.columns):

        if "PG" in Monm:
            try:
                y_prop_RarePT['Score_1T'] = y_prop_RarePT[1] + y_prop_RarePT[3]
                y_prop_RarePT['Score_0T'] = y_prop_RarePT[0] + y_prop_RarePT[2] + y_prop_RarePT[4]
            except:
                y_prop_RarePT['Score_1T'] = y_prop_RarePT[1] + y_prop_RarePT[3]
                y_prop_RarePT['Score_0T'] = y_prop_RarePT[0] + y_prop_RarePT[2]

            AUCi = AUCi + 1
            AUROC[AUCi] = metrics.roc_auc_score(y_RarePT, y_prop_RarePT[3])
            AUPR[AUCi] = metrics.average_precision_score(y_RarePT, y_prop_RarePT[3])
            CLASS[AUCi] = Monm + "_3"

            AUCi = AUCi + 1
            AUROC[AUCi] = metrics.roc_auc_score(y_RarePT, y_prop_RarePT['Score_1T'])
            AUPR[AUCi] = metrics.average_precision_score(y_RarePT, y_prop_RarePT['Score_1T'])
            CLASS[AUCi] = Monm + "_1T"

        elif "G" in Monm:
            test_predictions['Score_1T'] = test_predictions[1] + y_prop_RarePT[2]
            test_predictions['Score_0T'] = test_predictions[0]

            AUCi = AUCi + 1
            AUROC[AUCi] = metrics.roc_auc_score(y_RarePT, y_prop_RarePT[2])
            AUPR[AUCi] = metrics.average_precision_score(y_RarePT, y_prop_RarePT[2])
            CLASS[AUCi] = Monm + "_3"

            AUCi = AUCi + 1
            AUROC[AUCi] = metrics.roc_auc_score(y_RarePT, y_prop_RarePT['Score_1T'])
            AUPR[AUCi] = metrics.average_precision_score(y_RarePT, y_prop_RarePT['Score_1T'])
            CLASS[AUCi] = Monm + "_1T"

    Metric_test = pd.DataFrame()
    Metric_test["AUROC_"+Datanm] = pd.Series(AUROC)
    Metric_test["AUPR_"+Datanm] = pd.Series(AUPR)
    Metric_test["CLASS_"+Datanm] = pd.Series(CLASS)
    
    return Metric_test

def objective(trial):

    n_d = trial.suggest_int("n_d", 8, 64)
    n_steps = trial.suggest_int("n_steps", 3, 10)
    gamma = trial.suggest_float("gamma", 1.0, 2.0)
    n_independent = trial.suggest_int("n_independent", 1, 5)
    n_shared = trial.suggest_int("n_shared", 1, 2)
    momentum = trial.suggest_float("momentum", 0.01, 0.4)

    clf = TabNetClassifier(
        n_d=n_d, n_a=n_d, n_steps=n_steps,
        gamma=gamma, n_independent=n_independent, n_shared=n_shared,
        cat_idxs=cat_idxs,
        cat_dims=cat_dims,
        cat_emb_dim=1,
        lambda_sparse=1e-4, momentum=0.3, clip_value=2.,
        optimizer_fn=torch.optim.Adam,
        optimizer_params=dict(lr=2e-2),
        scheduler_params = {"gamma": 0.95, "step_size": 20},
        scheduler_fn=torch.optim.lr_scheduler.StepLR, epsilon=1e-15
    )

    max_epochs = 150 if not os.getenv("CI", False) else 2
    aug = ClassificationSMOTE(p=0.2)

    clf.fit(
        X_train=X_train, y_train=y_train,
        eval_set=[(X_train, y_train), (X_valid, y_valid)],
        eval_name=['train', 'valid'],
        max_epochs=max_epochs, patience=100,
        batch_size=1024 * 2, virtual_batch_size=256, weights=1,
        augmentations=aug
        )

    # To get final results you may need to use a mapping for classes 
    # as you are allowed to use targets like ["yes", "no", "maybe", "I don't know"]
    preds_mapper = { idx : class_name for idx, class_name in enumerate(clf.classes_)}
    print(f"BEST VALID SCORE FOR {dataset_name} : {clf.best_cost}")
    # or you can simply use the predict method
    y_prop_RarePT = pd.DataFrame(clf.predict_proba(X_RarePT))
    RarePT_metric = metics_cal(y_RarePT, y_prop_RarePT, "RareP")

    y_prop_RV300T = pd.DataFrame(clf.predict_proba(X_RV300T))
    RV300T_metric = metics_cal(y_RV300T, y_prop_RV300T, "RV300")

    y_prop_RV555T = pd.DataFrame(clf.predict_proba(X_RV555T))
    RV555T_metric = metics_cal(y_RV555T, y_prop_RV555T, "RV555")

    y_prop_RV800T = pd.DataFrame(clf.predict_proba(X_RV800T))
    RV800T_metric = metics_cal(y_RV800T, y_prop_RV800T, "RV800")

    Metric_combine = pd.concat([RarePT_metric, RV300T_metric, RV555T_metric, RV800T_metric], axis = 1)
    Metric_combine.to_csv("DNPGO_TrainGPU/optuna_metric/{}.tsv".format(trial.number), sep="\t", index=False)

    # Save a trained model to a file.
    with open("DNPGO_TrainGPU/optuna_5C/{}.pickle".format(trial.number), "wb") as fout:
        pickle.dump(clf, fout)
    return 1.0 - metrics.accuracy_score(y_valid, clf.predict(X_valid))

study = optuna.create_study()
study.optimize(objective, n_trials=100)



# # Load the best model.
# with open("{}.pickle".format(study.best_trial.number), "rb") as fin:
#     best_clf = pickle.load(fin)
# print(accuracy_score(y_valid, best_clf.predict(X_valid)))
