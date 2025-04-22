import gc
import os
import sys
import traceback

import numpy as np
import pandas as pd

try:
    import config.config_example as config
except ImportError:
    print("Error: config/config_example.py not found.")
    sys.exit(1)

CHR_COL = "CHR"
POS_COL = "POS"
P_COL = "P"
SNP_ID_COL = "SNP_ID"


def load_and_preprocess_gwas(file_path, chr_col_name, pos_col_name, p_col_name):
    print(f"Loading GWAS: {file_path}")
    if not os.path.exists(file_path):
        print(f"Error: File not found: {file_path}")
        return None

    df = None
    try:
        df_ws = pd.read_csv(
            file_path, sep=r"\s+", low_memory=False, comment="#", engine="python"
        )
        if p_col_name in df_ws.columns:
            df = df_ws
            print(f"Loaded whitespace separated. Shape: {df.shape}")
    except Exception:
        try:
            df_t = pd.read_csv(
                file_path, sep="\t", low_memory=False, comment="#"
            )
            if p_col_name in df_t.columns:
                df = df_t
                print(f"Loaded tab separated. Shape: {df.shape}")
        except Exception:
            pass

    if df is None:
        print(f"Error: Failed loading {file_path}")
        if config.SCZ_FILENAME in file_path:
            try:
                with open(file_path, "r") as f:
                    header = None
                    for line in f:
                        if line.startswith("#CHROM"):
                            header = line.strip().split("\t")
                            break
                    if header:
                        df = pd.read_csv(file_path, sep="\t", comment="#")
                        df.columns = header
                        print(f"Loaded with manual header. Shape: {df.shape}")

                    if df is None:
                        print("Manual header failed.")
                        return None

            except Exception as e_manual:
                print(f"Manual header read failed: {e_manual}")
                return None
        else:
            return None

    try:
        df.columns = df.columns.str.strip()

        required_cols_original = {chr_col_name, pos_col_name, p_col_name}
        if not required_cols_original.issubset(df.columns):
            missing = required_cols_original - set(df.columns)
            print(f"Error: Missing columns {missing} in {file_path}.")
            return None

        cols_to_select = {
            chr_col_name: CHR_COL,
            pos_col_name: POS_COL,
            p_col_name: P_COL,
        }
        df_selected = df[list(cols_to_select.keys())].copy()
        df_selected.rename(columns=cols_to_select, inplace=True)

        initial_rows = len(df_selected)
        print("Preprocessing...")

        df_selected[POS_COL] = pd.to_numeric(df_selected[POS_COL], errors="coerce")
        df_selected.dropna(subset=[POS_COL], inplace=True)
        df_selected[POS_COL] = df_selected[POS_COL].astype(int)

        df_selected[CHR_COL] = (
            df_selected[CHR_COL]
            .astype(str)
            .str.strip()
            .str.replace("chr", "", regex=False)
        )
        valid_chroms = [str(i) for i in range(1, 23)]
        df_selected = df_selected[df_selected[CHR_COL].isin(valid_chroms)].copy()
        if df_selected.empty:
            print("Error: No valid autosomal chromosomes found.")
            return None
        df_selected[CHR_COL] = df_selected[CHR_COL].astype(int)

        df_selected[P_COL] = pd.to_numeric(df_selected[P_COL], errors="coerce")
        df_selected.dropna(subset=[P_COL], inplace=True)
        df_selected = df_selected[
            (df_selected[P_COL] > 0) & (df_selected[P_COL] <= 1)
        ]

        initial_rows_before_dedup = len(df_selected)
        df_selected.drop_duplicates(
            subset=[CHR_COL, POS_COL], keep="first", inplace=True
        )
        print(f"Removed {initial_rows_before_dedup - len(df_selected)} duplicates.")

        df_selected[SNP_ID_COL] = (
            df_selected[CHR_COL].astype(str)
            + ":"
            + df_selected[POS_COL].astype(str)
        )

        df_selected.sort_values(by=[CHR_COL, POS_COL], inplace=True)

        rows_after = len(df_selected)
        print(f"Preprocessing finished. Shape: {df_selected.shape}")
        print(f"Rows removed: {initial_rows - rows_after}")

        return df_selected[[SNP_ID_COL, CHR_COL, POS_COL, P_COL]]

    except Exception as e:
        traceback.print_exc()
        print(f"Unexpected error during processing {file_path}: {e}")
        return None


if __name__ == "__main__":
    print("--- GWAS Preprocessing ---")

    os.makedirs(config.GWAS_PROC_OUTPUT_DIR, exist_ok=True)

    ad_file = os.path.join(config.GWAS_DIR, config.AD_FILENAME)
    ad_df = load_and_preprocess_gwas(
        ad_file, config.CHR_COL_AD, config.POS_COL_AD, config.P_COL_AD
    )

    if ad_df is None or ad_df.empty:
        sys.exit("Error processing AD data.")

    ad_sig = ad_df[ad_df[P_COL] <= config.P_VALUE_THRESHOLD].copy()
    print(f"Found {len(ad_sig):,} significant AD SNPs.")

    ad_out = config.AD_SNPS_FILE_PROC
    ad_input = ad_sig[[SNP_ID_COL, P_COL]].rename(columns={SNP_ID_COL: "SNP", P_COL: "P"})
    ad_input.to_csv(ad_out, sep="\t", index=False, header=True)
    print(f"Saved AD SNPs to: {ad_out}")

    del ad_df, ad_sig, ad_input
    gc.collect()

    scz_file = os.path.join(config.GWAS_DIR, config.SCZ_FILENAME)
    scz_df = load_and_preprocess_gwas(
        scz_file, config.CHR_COL_SCZ, config.POS_COL_SCZ, config.P_COL_SCZ
    )

    if scz_df is None or scz_df.empty:
        sys.exit("Error processing SCZ data.")

    scz_sig = scz_df[scz_df[P_COL] <= config.P_VALUE_THRESHOLD].copy()
    print(f"Found {len(scz_sig):,} significant SCZ SNPs.")

    scz_out = config.SCZ_SNPS_FILE_PROC
    scz_input = scz_sig[[SNP_ID_COL, P_COL]].rename(columns={SNP_ID_COL: "SNP", P_COL: "P"})
    scz_input.to_csv(scz_out, sep="\t", index=False, header=True)
    print(f"Saved SCZ SNPs to: {scz_out}")

    del scz_df, scz_sig, scz_input
    gc.collect()

    print("--- GWAS Preprocessing Finished ---")
