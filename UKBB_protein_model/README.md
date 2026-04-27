# Protein-based %GZMB+ prediction model

This repository contains an R workflow for training a protein-based lasso model to predict a single-cell phenotype (`%GZMB+`) across multiple cohorts. The script can optionally apply the trained model to UK Biobank (UKBB) proteomics data and export UKBB predictions.

The main script is:

```bash
protein_model_public.R
```

## Overview

The workflow performs the following steps:

1. Loads phenotype and proteomics data for the SoundLife, ABF300, KCL, and ALTRA cohorts.
2. Harmonizes shared protein features across cohorts.
3. Performs cohort-level batch correction using ComBat.
4. Generates PCA plots for the corrected proteomics matrix.
5. Trains a lasso model using `glmnet`.
6. Generates model performance plots and cross-validation summaries.
7. Optionally loads UKBB proteomics data and writes predicted `%GZMB+` values for UKBB participants.

## Input files

The core analysis requires four non-UKBB cohorts.

### Required input files

| Argument | Description |
| --- | --- |
| `--seattle-pct` | Sound Life single-cell phenotype file. Expected columns include `donor_id`, `minor_celltypes`, and `pct`. |
| `--seattle-olink` | Sound Life Olink proteomics file. Expected as an Excel file with columns such as `Subject`, `Baseline Age`, `sex`, `Visit`, `NPX_bridged`, and `Assay`. |
| `--stl-annotations` | ABF300 phenotype annotation file. The script expects columns corresponding to `ID`, `pheb`, `phek`, `temra`, `baseline age`, `sex`, and `ratio`. |
| `--stl-olink` | ABF300 Olink proteomics file in long format with `SampleID`, `SampleType`, `Assay`, and `NPX`. |
| `--kcl-annotations` | KCL annotation file with cell-type percentages and demographics. |
| `--kcl-olink` | KCL Olink proteomics file in long format with `SampleID`, `SampleType`, `Assay`, and `NPX`. |
| `--ra-annotations` | ALTRA phenotype annotation file. |
| `--ra-olink` | ALTRA Olink proteomics file with assay-level NPX values. |

### Optional UKBB input files

UKBB inputs are optional. If no UKBB input is provided, the script still trains the model, generates plots, and writes model outputs. UKBB prediction output is skipped.

If UKBB is provided, the script requires:

| Argument | Description |
| --- | --- |
| `--ukbb-metadata` | UKBB metadata file with `eid`, `SEX`, and `BL_AGE`. |
| `--ukbb-excluded` | File containing UKBB participant IDs to exclude. |
| `--ukbb-olink-file` | A UKBB proteomics file. 

## Example run with UKBB proteomics file

```bash
Rscript protein_model_public.R \
  --seattle-pct pct_Seattle.csv \
  --seattle-olink BRI_Olink_with_demographics.xlsx \
  --stl-annotations ABF300_pct_of_cells_ratios_noMAIT.csv \
  --stl-olink P24-006_Artyomov_3k_Extended_NPX.csv \
  --kcl-annotations all_percentages_all_cohort.csv \
  --kcl-olink P26-010_Artyomov_Extended_NPX.csv \
  --ra-annotations Donor_infor_TemK_B.csv \
  --ra-olink RA_Olink.csv \
  --ukbb-metadata UKBB_pheno.tsv \
  --ukbb-excluded w88907_20250818.csv \
  --ukbb-olink-file merged_ukbb_olink.tsv.gz
```

## Example run without UKBB

```bash
Rscript protein_model_public.R \
  --seattle-pct pct_Seattle.csv \
  --seattle-olink BRI_Olink_with_demographics.xlsx \
  --stl-annotations ABF300_pct_of_cells_ratios_noMAIT.csv \
  --stl-olink P24-006_Artyomov_3k_Extended_NPX.csv \
  --kcl-annotations all_percentages_all_cohort.csv \
  --kcl-olink P26-010_Artyomov_Extended_NPX.csv \
  --ra-annotations Donor_infor_TemK_B.csv \
  --ra-olink RA_Olink.csv
```

## Output files

By default, outputs are written to the `outputs/` directory.

| Output argument | Default path | Description |
| --- | --- | --- |
| `--out-union-cor1-f` | `outputs/UNION_cor1_f.txt` | Harmonized and ComBat-corrected analysis table used for model fitting. |
| `--out-ukbb-predictions` | `outputs/ukbb_pheb_predictions.txt` | UKBB predicted %GZMB+ values. This file is written only when UKBB inputs are provided. |
| `--out-r2-metrics` | `outputs/model_r2_metrics.txt` | Model R^2 metrics across cohorts and subsets. |
| `--out-coefficients` | `outputs/model_coefficients.txt` |  Model coefficients. |
| `--out-cross-validation` | `outputs/cross_validation_r2.txt` | Cross-validation R^2 results. |
| `--out-lasso-model` | `outputs/protein_lasso_model.rds` | Saved R object containing model outputs and metadata. |
| `--out-pca-plot` | `outputs/pca_plot.pdf` | PCA plot after ComBat correction. |
| `--out-r2-plot` | `outputs/r2_plot.pdf` | Observed versus predicted %GZMB+ plots. |
| `--out-crossval-plot` | `outputs/cross_validation_r2_plot.pdf` | Cross-validation R^2 distribution plot. |

Any output path can be overridden from the command line. For example:

```bash
Rscript protein_model_public.R \
  [required input arguments] \
  --out-ukbb-predictions results/ukbb_predictions.tsv \
  --out-lasso-model results/protein_model.rds \
  --out-pca-plot figures/pca_plot.pdf
```
