import os

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# --- Directory Paths (!!! EDIT THESE FOR YOUR ENVIRONMENT !!!) ---
GWAS_DIR = os.path.join(PROJECT_ROOT, "data", "raw_gwas")
ADSP_DIR = "/path/to/your/adsp/plink/files" # !!! EDIT THIS
PIPELINE_OUTPUT_DIR = os.path.join(PROJECT_ROOT, "pipeline_output_no_clumping")
ML_RESULTS_DIR = os.path.join(PROJECT_ROOT, "ml_results_iterative_knn_v1")

# --- Resource Files (!!! EDIT THIS FOR YOUR ENVIRONMENT !!!) ---
CHAIN_FILE = os.path.join(PROJECT_ROOT, "resources", "hg19ToHg38.over.chain.gz") # !!! EDIT THIS

# --- Executable Paths (Used by Python scripts if needed) ---
PYTHON_EXEC = "python3"
# PLINK_EXEC = "/path/to/plink" # Defined in 01_run_pipeline.sh
# LIFTOVER_EXEC = "/path/to/liftOver" # Defined in 01_run_pipeline.sh

# --- GWAS Processing Parameters ---
GWAS_PROC_OUTPUT_DIR = os.path.join(PIPELINE_OUTPUT_DIR, "processed_gwas")
AD_FILENAME = "PGCALZ2sumstatsExcluding23andMe.txt" # !!! EDIT IF DIFFERENT
SCZ_FILENAME = "PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv" # !!! EDIT IF DIFFERENT
CHR_COL_AD = "chr"
POS_COL_AD = "PosGRCh37"
P_COL_AD = "p"
CHR_COL_SCZ = "CHROM"
POS_COL_SCZ = "POS"
P_COL_SCZ = "PVAL"
P_VALUE_THRESHOLD = 5e-8

# --- ML Experiment Parameters ---
AD_SNPS_PREFIX = "adsp_ad_significant_snps_hg38"
SCZ_SNPS_PREFIX = "adsp_scz_significant_snps_hg38"
TEST_SIZE = 0.2
RANDOM_STATE = 42
KNN_IMPUTE_NEIGHBORS = 5
RUN_TABNET = True
N_VALUES = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 500]

# --- Model Hyperparameters ---
LR_PARAMS = {
    "logistic__C": 1.0,
    "logistic__solver": "liblinear",
    "logistic__random_state": RANDOM_STATE,
    "logistic__class_weight": "balanced",
}
RF_PARAMS = {
    "rf__n_estimators": 200,
    "rf__max_depth": 10,
    "rf__min_samples_leaf": 5,
    "rf__random_state": RANDOM_STATE,
    "rf__class_weight": "balanced",
}
XGB_PARAMS = {
    "xgb__n_estimators": 200,
    "xgb__learning_rate": 0.1,
    "xgb__max_depth": 5,
    "xgb__subsample": 0.8,
    "xgb__colsample_bytree": 0.8,
    "xgb__use_label_encoder": False,
    "xgb__eval_metric": "logloss",
    "xgb__random_state": RANDOM_STATE,
}
TABNET_PARAMS = {
    "n_d": 8, "n_a": 8, "n_steps": 3, "gamma": 1.3, "lambda_sparse": 1e-3,
    "optimizer_params": dict(lr=2e-2),
    "scheduler_params": {"step_size":10, "gamma":0.9},
    "mask_type": "sparsemax", "verbose": 0, "seed": RANDOM_STATE
}
TABNET_FIT_PARAMS = {
    "max_epochs": 50, "patience": 10, "batch_size": 1024, "eval_metric": ["auc"]
}

# --- Internal Use File Names ---
AD_SNPS_FILE_PROC = os.path.join(GWAS_PROC_OUTPUT_DIR, "ad_snps.txt")
SCZ_SNPS_FILE_PROC = os.path.join(GWAS_PROC_OUTPUT_DIR, "scz_snps.txt")
