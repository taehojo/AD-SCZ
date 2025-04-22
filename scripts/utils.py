import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer

try:
    import config.config_example as config
except ImportError:
    print("Error: config/config_example.py not found.")
    sys.exit(1)


def load_and_prepare_data(plink_prefix, pipeline_output_dir):
    print(f"Loading data: {plink_prefix}")
    start_time = time.time()
    raw_file = os.path.join(pipeline_output_dir, f"{plink_prefix}_raw.raw")

    if not os.path.exists(raw_file):
        print(f"Error: File not found: {raw_file}")
        return None, None, None, None
    try:
        df_raw = pd.read_csv(raw_file, sep=" ")
        if df_raw.shape[1] < 7:
            df_raw = pd.read_csv(raw_file, sep="\t")
    except Exception as e:
        print(f"Error loading {raw_file}: {e}")
        return None, None, None, None

    if df_raw.shape[1] < 7:
        print(f"Error: Unexpected columns in {raw_file}")
        return None, None, None, None

    y = df_raw["PHENOTYPE"].copy().replace({1: 0, 2: 1})
    valid_idx = y[y.isin([0, 1])].index
    if len(valid_idx) < len(y):
        print(f"Warning: Removing {len(y) - len(valid_idx)} invalid phenotypes.")
        y = y.loc[valid_idx]
        df_raw_filtered = df_raw.loc[valid_idx]
        if len(y) == 0:
            print("Error: No samples remaining.")
            return None, None, None, None
    else:
        df_raw_filtered = df_raw
    y = y.to_numpy()
    print(f"Phenotype counts (0=CN, 1=AD): {pd.Series(y).value_counts()}")

    snp_names = df_raw_filtered.columns[6:].tolist()
    X = df_raw_filtered.iloc[:, 6:].copy()
    print(f"Extracted {X.shape[1]} SNPs for {X.shape[0]} individuals.")
    if X.shape[1] == 0:
        print(f"Warning: No SNP columns found for {plink_prefix}.")
        return None, None, None, None

    print(f"Missing genotypes before imputation: {X.isna().sum().sum()}")
    X_np = X.to_numpy(dtype=float)
    if np.isnan(X_np).any():
        print(f"Imputing using KNN (k={config.KNN_IMPUTE_NEIGHBORS})...")
        impute_start = time.time()
        imputer = KNNImputer(n_neighbors=config.KNN_IMPUTE_NEIGHBORS)
        X_imputed = imputer.fit_transform(X_np)
        print(f"Missing after imputation: {np.isnan(X_imputed).sum()}")
        print(f"Imputation time: {time.time() - impute_start:.2f} sec")
    else:
        print("No missing genotypes found.")
        X_imputed = X_np

    sample_iids = df_raw_filtered["IID"].tolist()
    print(f"Data loading finished. Time: {time.time() - start_time:.2f} sec")
    return X_imputed, y, snp_names, sample_iids
