# AD/SCZ Cross-Disorder ML Prediction

## Overview

This project uses Machine Learning (ML) to test if adding Schizophrenia (SCZ) associated genetic variants (SNPs) improves Alzheimer's Disease (AD) risk prediction compared to using AD SNPs alone. It utilizes PGC AD/SCZ GWAS summary statistics and ADSP WGS reference genotypes. 

## Prerequisites

* Linux (Bash environment)
* PLINK (v1.9+) & UCSC liftOver tool
* Python 3.8+ (see `environment.yml` for packages)
* AD & SCZ GWAS summary statistics (GRCh37/hg19)
* Reference genotype panel (e.g., ADSP WGS - GRCh38)
* `hg19ToHg38.over.chain.gz` liftOver chain file

## Setup

1.  **Clone:**
    ```bash
    git clone https://github.com/taehojo/AD-SCZ.git
    cd your-project-name
    ```
2.  **Create Environment:**
    ```bash
    conda env create -f environment.yml
    conda activate ad_scz_env # Use environment name from environment.yml
    ```
3.  **Configure:**
    * Copy `config/config_example.py` to `config/config.py`.
    * **Edit `config/config.py`**: Set all required paths (GWAS, ADSP, CHAIN_FILE) and review parameters.
    * **Edit `scripts/01_run_pipeline.sh`**: Set paths for `PLINK_EXEC` and `LIFTOVER_EXEC`.
    * Ensure input data and chain file are placed according to the configured paths.

## Usage

1.  **Run Data Pipeline:**
    ```bash
    bash scripts/01_run_pipeline.sh
    ```
2.  **Run ML Experiment:**
    ```bash
    python scripts/03_iterative_ml.py
    ```

## Outputs

Results (performance tables, SNP rankings, importance scores, ROC plots) are saved in the directory specified by `ML_RESULTS_DIR` in `config.py`. Intermediate data pipeline files are saved in `PIPELINE_OUTPUT_DIR`.

## License

Copyright 2025, Taeho Jo and Carol S. Lim. All rights reserved.

## Citation

Jo, T., & Lim, C. S. (2025). Cross-disorder machine learning uncovers schizophrenia risk variants predictive of Alzheimer's disease.
