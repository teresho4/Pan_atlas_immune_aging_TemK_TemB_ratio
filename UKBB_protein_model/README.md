# Protein-based %GZMB+ prediction model

This repository contains an R workflow for training a protein-based lasso model to predict a single-cell phenotype (`%GZMB+`) across multiple cohorts. The script can optionally apply the trained model to UK Biobank (UKBB) proteomics data and export UKBB predictions.

The main script is:

```bash
protein_model.R
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

The core analysis requires four non-UKBB cohorts. Olink files are available here: https://www.synapse.org/Synapse:syn69762629

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
| `--ukbb-olink-file` | One pre-merged UKBB proteomics file in wide format. |
| `--ukbb-olink-instance-1` … `--ukbb-olink-instance-6` | Alternative mode: all six wide-format UKBB Olink files, merged by `eid`. |

### Exact input schemas and column meanings

Column names and category labels are case-sensitive. Phenotypes are percentages on the **0–100 scale**. Protein values are numeric Olink NPX measurements (normalized protein expression on a log2 scale).

#### 1. SoundLife phenotypes: `--seattle-pct` (`pct_Seattle.csv`)

Comma-separated, long-format CSV with a header; one row per donor and cell type.

| Column | Meaning / expected values | Script use |
| --- | --- | --- |
| `donor_id` | Donor identifier, e.g. `BR1004`. | Matches `Subject` in the SoundLife Olink file. |
| `minor_celltypes` | Cell-type label; must include exactly `Tem GZMB+` and `Tem GZMK+`. | Pivoted into the two phenotype columns. |
| `pct` | Percentage of cells assigned to that type. | `Tem GZMB+` becomes `pheb`; `Tem GZMK+` becomes `phek`. |

Demographics are taken from the Olink file instead.

#### 2. SoundLife proteomics: `--seattle-olink`

Excel workbook, **first sheet**, long format.

| Column | Meaning / expected values |
| --- | --- |
| `Subject` | Donor ID matching phenotype `donor_id`. |
| `Baseline Age` | Numeric baseline age in years; becomes `BL_AGE`. |
| `sex` | Sex label; exactly `male` maps to 1 and other non-missing labels map to 0. Use lowercase `male` / `female`. |
| `Visit` | Visit name. Main analysis uses exactly `Flu Year 1 Day 0`. |
| `Assay` | Protein/assay name; becomes a protein column. |
| `NPX_bridged` | Numeric bridged NPX measurement. Duplicate Subject/Assay values within a visit are averaged. |

The script also constructs tables for other visit labels, but the final SoundLife input uses baseline plus hard-coded exceptions: `BR1029`, `BR1035`, and `BR2004` use `Other - Non-Flu`; `BR2027` uses `Flu Year 1 Stand-Alone`.

#### 3. ABF300 phenotypes: `--stl-annotations` (`ABF300_pct_of_cells_ratios_noMAIT.csv`)

Comma-separated, wide-format CSV with a header. **Exactly eight columns in the following order are required**. Do not add an exported row-index column.

| Position | Supplied column | Meaning | Internal name |
| --- | --- | --- | --- |
| 1 | `Tube_id` | Sample/tube ID matching Olink `SampleID`. | `ID` |
| 2 | `%TemB` | Percentage of Tem GZMB+ cells. | `pheb` |
| 3 | `%TemK` | Percentage of Tem GZMK+ cells. | `phek` |
| 4 | `%Temra` | Percentage of Temra cells. | `temra` |
| 5 | `Age` | Age in years. | `BL_AGE` |
| 6 | `Sex` | `Male` / `Female`; `Male` maps to 1, other non-missing labels to 0. | `SEX` |
| 7 | `emK_emB` | Tem GZMK+ / Tem GZMB+ percentage ratio. | `rat` |
| 8 | `emB_emK` | Reciprocal ratio: Tem GZMB+ / Tem GZMK+. | `ratio` |


#### 4. ABF300 proteomics: `--stl-olink`

Comma-separated, long-format CSV with a header.

| Column | Meaning / script use |
| --- | --- |
| `SampleID` | Sample ID matching annotation `Tube_id`. |
| `SampleType` | Olink sample category; retained as a pivot identifier. Unlike KCL, there is no explicit `SAMPLE` filter here. |
| `Assay` | Protein/assay name, pivoted into columns; listed Olink control assays are removed. |
| `NPX` | Numeric NPX value; repeated measurements within SampleID/SampleType/Assay are averaged. |


#### 5. KCL phenotypes: `--kcl-annotations` (`all_percentages_all_cohort.csv`)

Comma-separated, long-format CSV with a header. Other cohorts may be present, but **only `Dataset == "KCL"` is used**.

| Column | Meaning / expected values | Script use |
| --- | --- | --- |
| `Dataset` | Cohort label; exactly `KCL`. | Selects KCL records. |
| `donor_id` | Donor ID matching KCL Olink `SampleID`. | Becomes `ID`. |
| `Annotation` | Cell-type label; must include `CD8 Tem GZMB+` and `CD8 Tem GZMK+`. | Pivoted into phenotype columns. |
| `percent` | Cell-type percentage. | GZMB+ becomes `pheb`, GZMK+ becomes `phek`; duplicates are averaged within pivot identifiers. |
| `Sex` | Exactly `Male` or `Female`. | Maps to 1 or 0; unrecognized labels become NA. |
| `Ethnicity` | Donor ethnicity annotation. | Required by the initial selection and retained during pivoting, but dropped before modeling. |
| `Age` | Numeric age in years. | Becomes `BL_AGE`. |


#### 6. KCL proteomics: `--kcl-olink`

Comma-separated, long-format CSV with a header. Requires **`SampleID`, `SampleType`, `Assay`, and `NPX`**, with the meanings described for ABF300. `SampleID` must match KCL `donor_id`. Only records with exactly `SampleType == "SAMPLE"` are retained after pivoting. 

#### 7. ALTRA phenotypes: `--ra-annotations` (`Donor_infor_TemK_B.csv`)

Comma-separated CSV with a header. Preserve **these eight columns in this order**. 

| Position | Column | Meaning / script use |
| --- | --- | --- |
| 1 | `sample.sampleKitGuid` | Sample-kit ID matching the ALTRA Olink file; becomes `ID`. |
| 2 | `Tem_GZMK` | Tem GZMK+ percentage; becomes `phek`. |
| 3 | `Tem_GZMB` | Tem GZMB+ percentage; becomes `pheb`. |
| 4 | `RatioKB` | Tem GZMK+ / Tem GZMB+ ratio; becomes `rat`. |
| 5 | `subject.diseaseGroup` | Disease-group label; becomes `disease`. Exactly `Control (HC1)` identifies healthy controls used in training. Other groups form the unhealthy subset. |
| 6 | `subject.progressionGroup` | Progression-group annotation; not selected as a final model covariate. |
| 7 | `sample.visitName` | Visit annotation; this loader does not filter by visit. |
| 8 | `subject.subjectGuid` | Subject identifier; the merge uses sample-kit ID instead. |


#### 8. ALTRA proteomics: `--ra-olink`

Comma-separated, long-format CSV with a header.

| Column | Meaning / script use |
| --- | --- |
| `sample.sampleKitGuid` | Sample-kit ID matching the phenotype table. |
| `subject.biologicalSex` | Sex label; exactly `Male` maps to 1, other non-missing labels to 0. |
| `sample.subjectAgeAtDraw` | Numeric age at blood draw in years; becomes `BL_AGE`. |
| `olink.assay` | Protein/assay name, pivoted into protein columns. |
| `olink.NPX_norm` | Numeric normalized NPX; repeated assay values are averaged within sample/sex/age groups. |

#### 9. Optional UKBB metadata: `--ukbb-metadata`

| Column | Meaning / expected values |
| --- | --- |
| `eid` | UKBB participant ID matching the proteomics file. |
| `SEX` | The script assumes **1 = male, 2 = female**, then recodes 2 to 0 and other non-missing values to 1. A 0/1-coded input would therefore incorrectly map 0 to male. |
| `BL_AGE` | Numeric baseline age in years. |

#### 10. Optional UKBB exclusions: `--ukbb-excluded`

**Tab-separated text without a header**; participant IDs must be in the **first column**. The script renames this column to `eid`; other columns, if present, are ignored. Despite the `.csv` filename in the example, this reader explicitly uses a tab delimiter. An empty header-only file is not the expected format.

#### 11. Optional UKBB proteomics: one file or six files

Use either `--ukbb-olink-file` **or all six** `--ukbb-olink-instance-1` through `--ukbb-olink-instance-6`

| Mode | Expected columns and layout |
| --- | --- |
| One merged file | Headered wide-format table: one ID column (`eid` or `ID`) and numeric protein columns named by assay/protein. If neither recognized name is present after uppercasing, the first column is treated as `ID`. |
| Six instance files | Each is a headered wide-format table with an exact lowercase `eid` column plus numeric protein columns. The six tables are merged by `eid` using inner joins; only participants present in all six remain. Avoid overlapping protein column names that would receive merge suffixes. |


### Internal phenotype names and matching rules

| Internal name | Meaning |
| --- | --- |
| `ID` | Donor/sample identifier used for cohort-specific phenotype–proteomics joins. |
| `pheb` | Observed Tem GZMB+ percentage; the lasso prediction target. |
| `phek` | Observed Tem GZMK+ percentage. |
| `rat` | Tem GZMK+ / Tem GZMB+ ratio. |
| `BL_AGE` | Age in years, taken from the cohort-specific source described above. |
| `SEX` | Harmonized sex: 1 = male, 0 = female, subject to the input coding rules above. |
| `COHORT` | Internal cohort label: `SEA`, `STL`, `KCL`, `RA`, or `UKBB`; `RAH` is used as a training-selection alias for healthy ALTRA controls. |


## Example run with UKBB proteomics file

```bash
Rscript protein_model.R \
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
Rscript protein_model.R \
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

