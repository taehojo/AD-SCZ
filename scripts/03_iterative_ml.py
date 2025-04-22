import os
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, auc, f1_score, roc_auc_score,
                             roc_curve)
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

try:
    import config.config_example as config
except ImportError:
    print("Error: config/config_example.py not found.")
    sys.exit(1)

try:
    from utils import load_and_prepare_data
except ImportError:
    print("Error: scripts/utils.py not found.")
    sys.exit(1)

TabNetClassifier = None
torch = None
if config.RUN_TABNET:
    try:
        from pytorch_tabnet.tab_model import TabNetClassifier
        import torch
        config.TABNET_PARAMS["optimizer_fn"] = torch.optim.Adam
        config.TABNET_PARAMS["scheduler_fn"] = torch.optim.lr_scheduler.StepLR
    except ImportError:
        print("Warning: TabNet/PyTorch import failed. Skipping TabNet.")
        config.RUN_TABNET = False
else:
    print("TabNet skipped by configuration.")


def main():
    output_dir = config.PIPELINE_OUTPUT_DIR
    results_dir = config.ML_RESULTS_DIR
    delong_data_dir = os.path.join(results_dir, "delong_test_data")

    lr_pipeline = Pipeline([("scaler", StandardScaler()), ("logistic", LogisticRegression())])
    lr_pipeline.set_params(**config.LR_PARAMS)

    rf_pipeline = Pipeline([("rf", RandomForestClassifier())])
    rf_pipeline.set_params(**config.RF_PARAMS)

    xgb_pipeline = Pipeline([("xgb", XGBClassifier())])
    xgb_pipeline.set_params(**config.XGB_PARAMS)

    OPTIMAL_N_INIT = {"Logistic Regression": 0, "Random Forest": 40, "XGBoost": 70, "TabNet": 50}

    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(delong_data_dir, exist_ok=True)

    print("===== Starting Iterative Feature Addition Experiment =====")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    print("--- Loading Data ---")
    X_ad, y_ad, snp_names_ad, iids_ad = load_and_prepare_data(config.AD_SNPS_PREFIX, output_dir)
    X_scz, y_scz, snp_names_scz, iids_scz = load_and_prepare_data(config.SCZ_SNPS_PREFIX, output_dir)

    if X_ad is None or X_scz is None: sys.exit("Error: Failed to load AD/SCZ data.")
    if not np.array_equal(y_ad, y_scz) or iids_ad != iids_scz: sys.exit("Error: Phenotype/ID mismatch.")

    y_full = y_ad
    n_samples, n_feat_ad = X_ad.shape
    _, n_feat_scz = X_scz.shape

    summary_data = {"Metric": ["Samples", "AD SNPs", "SCZ SNPs"], "Value": [n_samples, n_feat_ad, n_feat_scz]}
    summary_df = pd.DataFrame(summary_data)
    print("--- Dataset Summary ---")
    print(summary_df.to_string(index=False))
    summary_df.to_csv(os.path.join(results_dir, "table1_dataset_summary.csv"), index=False)

    max_n = len(snp_names_scz)
    N_VALUES = sorted(list(set([n for n in config.N_VALUES if n <= max_n])))
    print(f"Evaluating N SCZ SNPs: {N_VALUES}")

    print("--- Train/Test Split ---")
    indices = np.arange(n_samples)
    train_idx, test_idx, y_train, y_test = train_test_split(
        indices, y_full, test_size=config.TEST_SIZE, random_state=config.RANDOM_STATE, stratify=y_full
    )
    X_train_ad, X_test_ad = X_ad[train_idx, :], X_ad[test_idx, :]
    X_train_scz, X_test_scz = X_scz[train_idx, :], X_scz[test_idx, :]
    print(f"Train: {len(y_train)}, Test: {len(y_test)}")

    y_test_path = os.path.join(delong_data_dir, "y_test.npy")
    np.save(y_test_path, y_test)
    print(f"Saved y_test to {y_test_path}")

    print("--- Rank SCZ SNPs ---")
    start_rank = time.time()
    rank_model = XGBClassifier(**{k.replace("xgb__", ""): v for k, v in config.XGB_PARAMS.items()})
    rank_model.fit(X_train_scz, y_train)
    print(f"SCZ ranking model trained. Time: {time.time() - start_rank:.2f} sec")
    importances = rank_model.feature_importances_
    scz_imp_df = pd.DataFrame({"SNP": snp_names_scz, "Importance": importances})
    scz_imp_df.sort_values(by="Importance", ascending=False, inplace=True)
    ranked_scz_names = scz_imp_df["SNP"].tolist()
    scz_imp_df.to_csv(os.path.join(results_dir, "table_scz_snp_ranking.csv"), index=False, float_format="%.6g")
    print("Saved SCZ SNP importance ranking.")

    print("--- Iterative Evaluation ---")
    models = {"Logistic Regression": lr_pipeline, "Random Forest": rf_pipeline, "XGBoost": xgb_pipeline}
    if config.RUN_TABNET and TabNetClassifier is not None: models["TabNet"] = TabNetClassifier

    all_results = []
    roc_data = {name: {} for name in models}
    scz_idx_map = {name: i for i, name in enumerate(snp_names_scz)}
    actual_optimal_n = {}

    for n in N_VALUES:
        print(f"===== Processing N = {n} =====")
        start_iter = time.time()
        if n == 0:
            top_n_names = []
            X_train_scz_n = X_train_scz[:, :0]
            X_test_scz_n = X_test_scz[:, :0]
        else:
            curr_n = min(n, len(ranked_scz_names))
            if curr_n != n: print(f"Adjusted N: {n} -> {curr_n}")
            top_n_names = ranked_scz_names[:curr_n]
            top_n_idx = [scz_idx_map[name] for name in top_n_names if name in scz_idx_map]
            if len(top_n_idx) != curr_n: print(f"Warning: Found {len(top_n_idx)} indices for N={curr_n}")
            X_train_scz_n = X_train_scz[:, top_n_idx]
            X_test_scz_n = X_test_scz[:, top_n_idx]

        n_ad_feat = X_train_ad.shape[1]
        n_scz_feat = X_train_scz_n.shape[1]
        total_feat = n_ad_feat + n_scz_feat
        print(f"Combined features: {n_ad_feat} AD + {n_scz_feat} SCZ = {total_feat} Total")
        X_train_curr = np.hstack((X_train_ad, X_train_scz_n))
        X_test_curr = np.hstack((X_test_ad, X_test_scz_n))
        curr_snp_names = snp_names_ad + top_n_names

        for model_name, model_obj_or_class in models.items():
            start_model = time.time()
            print(f"--- Evaluating {model_name} (N={n}) ---")
            auc_val, acc_val, f1_val = np.nan, np.nan, np.nan
            fpr, tpr = None, None
            model_final = None
            y_prob = None

            try:
                if model_name == "TabNet":
                    tabnet_model = TabNetClassifier(**config.TABNET_PARAMS)
                    tabnet_model.fit(X_train_curr, y_train, eval_set=[(X_test_curr, y_test)], **config.TABNET_FIT_PARAMS)
                    y_prob = tabnet_model.predict_proba(X_test_curr)[:, 1]
                    y_pred = (y_prob > 0.5).astype(int)
                    model_final = tabnet_model
                else:
                    pipeline = model_obj_or_class
                    pipeline.fit(X_train_curr, y_train)
                    y_prob = pipeline.predict_proba(X_test_curr)[:, 1]
                    y_pred = pipeline.predict(X_test_curr)
                    model_final = pipeline

                if y_prob is not None:
                    auc_val, acc_val, f1_val = roc_auc_score(y_test, y_prob), accuracy_score(y_test, y_pred), f1_score(y_test, y_pred)
                    print(f"{model_name}: AUC={auc_val:.4f}, Acc={acc_val:.4f}, F1={f1_val:.4f}")
                    fpr, tpr, _ = roc_curve(y_test, y_prob)
                else: print(f"{model_name}: Probability prediction failed.")

                is_base = n == 0
                is_opt_init = (model_name in OPTIMAL_N_INIT and n == OPTIMAL_N_INIT[model_name])

                if (is_base or is_opt_init) and y_prob is not None:
                    safe_name = model_name.replace(" ", "_").lower()
                    prob_fname = f"proba_{safe_name}_n{n}.npy"
                    prob_path = os.path.join(delong_data_dir, prob_fname)
                    np.save(prob_path, y_prob)
                    print(f"Saved {model_name} (N={n}) probabilities to {prob_path}")

                if n > 0 and model_final is not None and model_name in ["XGBoost", "TabNet"]:
                    imps = None
                    try:
                        if model_name == "XGBoost":
                            model_imp = model_final.named_steps["xgb"]
                            if hasattr(model_imp, "feature_importances_"): imps = model_imp.feature_importances_
                        elif model_name == "TabNet":
                            model_imp = model_final
                            if hasattr(model_imp, "feature_importances_"): imps = model_imp.feature_importances_

                        if imps is not None:
                            temp_imp_df = pd.DataFrame({"SNP": curr_snp_names, "Importance": imps})
                            scz_imp_filt = temp_imp_df[temp_imp_df["SNP"].isin(top_n_names)].copy()
                            scz_imp_filt.sort_values(by="Importance", ascending=False, inplace=True)
                            imp_fname = os.path.join(results_dir, f"importance_scz_{model_name.replace(' ', '_')}_n{n}.csv")
                            scz_imp_filt.to_csv(imp_fname, index=False, float_format="%.6g")
                    except Exception as e_imp: print(f"Importance extraction failed for {model_name} N={n}: {e_imp}")

            except Exception as e:
                print(f"Error during {model_name} N={n}: {e}")
                if 'auc_val' not in locals() or auc_val is np.nan: acc_val, f1_val = np.nan, np.nan
                if y_prob is None: print(f"Could not get y_prob for {model_name} N={n}.")

            all_results.append({"Model": model_name, "N_Added_SCZ_SNPs": n, "Total_SNPs": total_feat, "AUC": auc_val, "Accuracy": acc_val, "F1_Score": f1_val})
            if fpr is not None and tpr is not None: roc_data[model_name][n] = (fpr, tpr)
            print(f"Finished {model_name} N={n}. Time: {time.time() - start_model:.2f} sec")
        print(f"Finished N={n}. Iteration time: {time.time() - start_iter:.2f} sec")

    print("--- Saving Results ---")
    results_df = pd.DataFrame(all_results)
    results_df.sort_values(by=["Model", "N_Added_SCZ_SNPs"], inplace=True)
    results_path = os.path.join(results_dir, "table_all_results_by_N.csv")
    results_df.to_csv(results_path, index=False, float_format="%.4f")
    print(f"Results saved to: {results_path}")

    print("--- Find Optimal N ---")
    best_results = {}
    print("--- Best N (Max AUC) per Model ---")

    for model_name in models:
        model_res = results_df[(results_df["Model"] == model_name) & results_df["AUC"].notna()].copy()
        if model_res.empty:
            print(f"No valid results for {model_name}.")
            best_results[model_name] = {"N": np.nan, "AUC": np.nan, "Baseline_AUC": np.nan, "MaxN_AUC": np.nan}
            actual_optimal_n[model_name] = np.nan
            continue

        best_run = model_res.loc[model_res["AUC"].idxmax()]
        best_n = int(best_run["N_Added_SCZ_SNPs"])
        best_auc = best_run["AUC"]
        actual_optimal_n[model_name] = best_n

        base_run = model_res[model_res["N_Added_SCZ_SNPs"] == 0]
        base_auc = base_run["AUC"].iloc[0] if not base_run.empty else np.nan
        max_n_val = max(N_VALUES) if N_VALUES else 0
        max_n_run = model_res[model_res["N_Added_SCZ_SNPs"] == max_n_val]
        max_n_auc = max_n_run["AUC"].iloc[0] if not max_n_run.empty else np.nan

        best_results[model_name] = {"N": best_n, "AUC": best_auc, "Baseline_AUC": base_auc, "MaxN_AUC": max_n_auc}
        print(f"Model: {model_name}, Best N={best_n}, Best AUC={best_auc:.4f}, Baseline AUC={base_auc:.4f}")

    print("--- Performance Summary Table ---")
    summary_rows = []
    for model_name, res in best_results.items():
        improve = (res["AUC"] - res["Baseline_AUC"]) if not np.isnan(res["AUC"]) and not np.isnan(res["Baseline_AUC"]) else np.nan
        summary_rows.append({"Model": model_name, "Baseline AUC (N=0)": res["Baseline_AUC"], "Best N": res["N"], "Best AUC": res["AUC"], f"AUC at Max N ({max_n_val})": res["MaxN_AUC"], "Improvement": improve})
    summary_df = pd.DataFrame(summary_rows).sort_values(by="Best AUC", ascending=False)
    print(summary_df.round(4).to_string(index=False))
    summary_df.to_csv(os.path.join(results_dir, "table2_3_perf_summary.csv"), index=False, float_format="%.4f")

    overall_best_auc, overall_best_n, overall_best_model = -1.0, None, None
    if best_results:
        valid_models = [m for m, r in best_results.items() if not np.isnan(r["AUC"])]
        if valid_models:
            overall_best_model = max(valid_models, key=lambda m: best_results[m]["AUC"])
            overall_best_n = actual_optimal_n.get(overall_best_model)
            overall_best_auc = best_results[overall_best_model]["AUC"]
            print(f"\nOverall Best: Model={overall_best_model}, N={overall_best_n}, AUC={overall_best_auc:.4f}")
        else: print("\nCannot determine overall best model.")
    else: print("\nCannot determine overall best model.")

    if overall_best_model is not None and not np.isnan(overall_best_n):
        print(f"--- Generating ROC Plot (Best Model: {overall_best_model}) ---")
        plt.figure(figsize=(8, 6))
        roc_fname = f"figure_roc_{overall_best_model}_bestN{overall_best_n}.png"

        fpr_ad, tpr_ad = roc_data[overall_best_model].get(0, (None, None))
        if fpr_ad is not None: plt.plot(fpr_ad, tpr_ad, lw=2, label=f'AD Only (N=0) (AUC = {auc(fpr_ad, tpr_ad):.3f})')

        print(f"Training {overall_best_model} on SCZ-only...")
        y_prob_scz = None
        try:
            scz_model_obj = models[overall_best_model]
            if overall_best_model == "TabNet" and config.RUN_TABNET and TabNetClassifier:
                if not np.isnan(best_results['TabNet']['AUC']):
                    scz_model = TabNetClassifier(**config.TABNET_PARAMS)
                    scz_model.fit(X_train_scz, y_train, eval_set=[(X_test_scz, y_test)], **config.TABNET_FIT_PARAMS)
                    y_prob_scz = scz_model.predict_proba(X_test_scz)[:, 1]
                else: print("Skipping failed SCZ-only TabNet.")
            elif overall_best_model in ["Logistic Regression", "Random Forest", "XGBoost"]:
                scz_pipeline = scz_model_obj
                scz_pipeline.fit(X_train_scz, y_train)
                y_prob_scz = scz_pipeline.predict_proba(X_test_scz)[:, 1]

            if y_prob_scz is not None:
                safe_name = overall_best_model.replace(" ", "_").lower()
                prob_fname = f"proba_scz_only_{safe_name}.npy"
                prob_path = os.path.join(delong_data_dir, prob_fname)
                np.save(prob_path, y_prob_scz)
                print(f"Saved SCZ-only probabilities to {prob_path}")
                fpr_scz, tpr_scz, _ = roc_curve(y_test, y_prob_scz)
                plt.plot(fpr_scz, tpr_scz, lw=2, label=f'SCZ Only (AUC = {auc(fpr_scz, tpr_scz):.3f})')
            else: print(f"Could not get SCZ-only probabilities for {overall_best_model}.")
        except Exception as e_scz: print(f"Error training/evaluating SCZ-only {overall_best_model}: {e_scz}")

        best_n_actual = actual_optimal_n.get(overall_best_model)
        if best_n_actual is not np.nan and best_n_actual is not None:
            fpr_best, tpr_best = roc_data[overall_best_model].get(best_n_actual, (None, None))
            if fpr_best is not None: plt.plot(fpr_best, tpr_best, lw=2, label=f'AD + Top {best_n_actual} SCZ (AUC = {auc(fpr_best, tpr_best):.3f})')
            else: print(f"ROC data for Best N={best_n_actual} not found.")
        else: print(f"Cannot plot Best N curve for {overall_best_model}.")

        plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        plt.xlim([0.0, 1.0]); plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate'); plt.ylabel('True Positive Rate')
        plt.title(f'Comparative ROC Curves ({overall_best_model})')
        plt.legend(loc="lower right", fontsize='small'); plt.grid(True)
        roc_plot_path = os.path.join(results_dir, roc_fname)
        try:
            plt.savefig(roc_plot_path, bbox_inches='tight')
            print(f"ROC plot saved to: {roc_plot_path}")
        except Exception as e_plot: print(f"Error saving ROC plot: {e_plot}")
        plt.close()
    else: print("Cannot generate ROC plot: Best model/N not determined.")

    print(f"===== Script Finished at {time.strftime('%Y-%m-%d %H:%M:%S')} =====")

if __name__ == "__main__":
    main()
