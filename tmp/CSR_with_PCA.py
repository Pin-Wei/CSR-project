#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os
import glob
# import math
import pickle
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold
from scipy.linalg import pinv

# from scipy.stats import pearsonr 
# from statsmodels.stats.multitest import fdrcorrection
# import seaborn as sns
# import matplotlib.pyplot as plt
# import matplotlib.ticker as ticker

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

class Config:
    def __init__(self):
        self.n_components = 3
        self.seed = 42 # np.random.randint(0, 10000)

        ## The directory storing the original data:
        self.path_ver = ["local", "server"][1]
        if self.path_ver == "local":
            self.hsk_data_dir = os.path.join("..", "my_Haskins_project", "Data", "Chang_et_al") 
        else: # server
            self.hsk_data_dir = os.path.join("..", "Haskins_project", "Data_single_characters")
        
        self.master_path = os.path.join(self.hsk_data_dir, "Chang_2020_z_CSR.xlsx") 
            ## made through "Modify_data_Chang_for_CSR.R", originated from "AoA_character_naming.rds"
            ## 1,353 characters; columns = ["Char", "subject_id", "RT", "LogCF", "NS", "CON", "PC", "SC", "SAR", "IMG", "AoA"]

        ## The directory to store PCA results and transformed data for CSR function fitting:
        self.csr_data_dir = os.path.join("Data_Linguistic", "Chang_et_al", "Naming")
        self.pca_data_dir = os.path.join(self.csr_data_dir, f"PCA_nc={self.n_components}_seed={self.seed}")
        os.makedirs(self.pca_data_dir, exist_ok=True)

        self.pca_model_path = os.path.join(self.pca_data_dir, "model.pkl")
        self.pc_loadings_path = os.path.join(self.pca_data_dir, "pc_loadings.xlsx")
        self.explained_vars_path = os.path.join(self.pca_data_dir, "pc_explained_vars.xlsx")
        self.char_data_path = os.path.join(self.pca_data_dir, "char_data.xlsx")
        self.subj_data_path_template = os.path.join(self.pca_data_dir, "subj_{}.xlsx")

        ## The directory to store CSR fitting results:
        self.csr_result_dir = os.path.join("Output_Linguistic", "Chang_et_al", "Naming")
        os.makedirs(self.csr_result_dir, exist_ok=True)
        self.csr_result_path = os.path.join(self.csr_result_dir, f"CSR_results_nc={self.n_components}_seed={self.seed}.xlsx")
        self.csr_summary_path = self.csr_result_path.replace(".xlsx", "_summ.xlsx")

def perform_pca(X: pd.DataFrame) -> pd.DataFrame:
    ## Fit PCA model
    pca = PCA(n_components=config.n_components, random_state=config.seed)
    pca.fit(X)

    ## Save PCA model
    with open(config.pca_model_path, "wb") as f:
        pickle.dump(pca, f)

    ## Save PCA loadings
    eignvectors = pca.components_.T
    eignvalues = pca.explained_variance_
    loadings = pd.DataFrame(
        data=eignvectors * np.sqrt(eignvalues), 
        columns=[f"PC{i+1}" for i in range(config.n_components)], 
        index=X.columns
    )
    loadings.to_excel(config.pc_loadings_path)

    ## Save explained variance
    explained = pd.DataFrame(
        {'Explained Variance': pca.explained_variance_, 
        'Explained Variance Ratio': pca.explained_variance_ratio_}, 
        index=[f"PC{i+1}" for i in range(config.n_components)]
    )
    explained.to_excel(config.explained_vars_path)

    ## Reduce original data
    reduced_X = pca.transform(X)
    reduced_X = pd.DataFrame(
        data=reduced_X, 
        columns=[f"PC{i+1}" for i in range(config.n_components)], 
        index=X.index
    )
    reduced_X.reset_index(inplace=True)
    reduced_X.rename(columns={"index": "Char"})
    reduced_X.to_excel(config.char_data_path, index=False)

    return reduced_X

def fit_csr_function(df: pd.DataFrame):
    def _parse_df(df: pd.DataFrame):
        F_vals = df.iloc[:, :-1].values.astype(float) # feature values
        Y = df.iloc[:, -1].values.astype(float)
        nF = F_vals.shape[1] # number of features
        nT = len(Y) # number of trials
        nP = int((nF * 3 + nF ** 2 + 2) / 2) # number of parameters
        nI = sum(range(1, nF)) # number of interaction terms

        return  F_vals, Y, nF, nT, nP, nI

    def _construct_design_matrix(nT: int, nF: int, nI: int, F_vals: np.ndarray) -> np.ndarray:
        X = np.zeros((nT, 1 + nF + nF + nI))
        X[:, 0] = 1 # constant term
        X[:, 1:1 + nF] = F_vals # linear terms
        X[:, 1 + nF:1 + 2 * nF] = F_vals ** 2 # quadratic terms
        col = 1 + nF + nF
        for j in range(nF - 1):
            for k in range(j + 1, nF):
                X[:, col] = F_vals[:, j] * F_vals[:, k] # interaction terms
                col += 1

        return X
    
    def _evalute_model(Coefs: np.ndarray, X: np.ndarray, Y: np.ndarray, nT: int, nP: int) -> tuple:
        Y_pred = X @ Coefs
        residuals = Y - Y_pred
        RSS = np.sum(residuals ** 2)
        TSS = np.sum((Y - np.mean(Y)) ** 2)
        Rsq = 1 - RSS / TSS
        Rsq_adj = 1 - (RSS / (nT - nP - 1)) / (TSS / (nT - 1))
        # re_var = RSS / nT # residual error variance
        # log_lik = -0.5 * nT * (np.log(2 * math.pi * re_var) + 1) # log-likelihood
        # AIC = -2 * log_lik + 2 * nP
        # AICc = AIC + (2 * nP ** 2 + 2 * nP) / (nT - nP - 1)
        # BIC = -2 * log_lik + np.log(nT) * nP
        RMSE = np.sqrt(np.mean(residuals ** 2)) # root mean square error
        NRMSE = RMSE / (Y.max() - Y.min()) # normalized 

        return Y_pred, residuals, [Rsq, Rsq_adj, NRMSE]
    
    def build_coef_names(feature_names: list, nF: int) -> list:
        coef_names = ["X0"] +  feature_names
        coef_names.extend([f"F{i+1}^2" for i in range(nF)])
        for j in range(nF - 1):
            for k in range(j + 1, nF):
                coef_names.append(f"F{j+1}F{k+1}")

        return coef_names

    kf = KFold(n_splits=5, shuffle=True, random_state=config.seed)
    result_list = []

    for fold, (train_index, test_index) in enumerate(kf.split(df)):
        df_train, df_test = df.iloc[train_index], df.iloc[test_index]

        ## Create design matrix and solve regression coefficients:
        F_vals_train, Y_train, nF, nT, nP, nI = _parse_df(df_train)
        X_train = _construct_design_matrix(nT, nF, nI, F_vals_train)
        Coefs = pinv(X_train.T @ X_train) @ (X_train.T @ Y_train)
        Y_pred_train, re_train, model_evals_train = _evalute_model(Coefs, X_train, Y_train, nT, nP)

        ## Evaluate on test set:
        F_vals_test, Y_test, _, nT_, _, _ = _parse_df(df_test)
        X_test = _construct_design_matrix(nT_, nF, nI, F_vals_test)
        Y_pred_test, re_test, model_evals_test = _evalute_model(Coefs, X_test, Y_test, nT_, nP)

        ## Store results:
        results = {"Fold": fold + 1}
        coef_names = build_coef_names(list(df.columns[:-1]), nF)
        results.update(dict(zip(coef_names, Coefs)))
        results.update({
            "Train_Rsq": model_evals_train[0],
            "Train_Rsq_adj": model_evals_train[1], 
            "Train_NRMSE": model_evals_train[-1], 
            "Test_Rsq": model_evals_test[0], 
            "Test_NRMSE": model_evals_test[-1]
        })
        results = pd.DataFrame(results, index=[0])
        result_list.append(results)

    return pd.concat(result_list, ignore_index=True)

## ============================================================================================================
if __name__ == "__main__":
    config = Config()
    overwrite = False

    ## Create dataset with PCA-transformed features for each participant:
    if len(glob.glob(config.subj_data_path_template.format("*"))) == 0 or overwrite:

        master_df = pd.read_excel(config.master_path)
        subj_list = list(master_df["subject_id"].unique())

        ## Perform PCA on single character psycholinguistic dataset
        if not os.path.exists(config.char_data_path) or overwrite:
            char_df = master_df.loc[:, ["Char", "LogCF", "NS", "CON", "PC", "SC", "SAR", "IMG", "AoA"]].drop_duplicates()
            char_df = char_df.set_index("Char")
            reduced_char_df = perform_pca(char_df) # perform PCA and save results
        else:
            reduced_char_df = pd.read_excel(config.char_data_path)

        ## Save PCA-transformed data for each participant
        updated_master = (
            master_df.loc[:, ["Char", "subject_id", "RT"]]
            .merge(reduced_char_df, on="Char", how="left")
        )
        for sid in subj_list:
            subj_df = updated_master[updated_master["subject_id"] == sid]
            df = subj_df.loc[:, [f"PC{i+1}" for i in range(config.n_components)] + ["RT"]]
            df.to_excel(config.subj_data_path_template.format(sid), index=False)

    ## Fit CSR function for each participant
    df_path_list = glob.glob(config.subj_data_path_template.format("*"))
    output_list = []

    for df_path in df_path_list:
        df = pd.read_excel(df_path)
        sid = os.path.basename(df_path).split("_")[1].split(".")[0]

        results_df = fit_csr_function(df) # fit CSR function with 5-fold cross-validation
        results_df.insert(0, "SID", sid)
        output_list.append(results_df)

    output_df = pd.concat(output_list, ignore_index=True)
    output_df.to_excel(config.csr_result_path, index=False)

    data_desc = (
        output_df
        .groupby("SID")
        .mean()
        .iloc[:, 2::]
        .describe()
        .astype("float")
        .map(lambda x: f"{x:.3f}")
        .loc[["min", "mean", "max", "std"], :]
        .T
    )
    data_desc.to_excel(config.csr_summary_path)
