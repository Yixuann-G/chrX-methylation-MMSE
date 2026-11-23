# X-Chromosome DNA Methylation and Cognitive Variation in Centenarians
## 1. Project overview

This study includes centenarians from the Dutch 100-plus Study.
It uses whole-blood PacBio DNA methylation data.
It examines the association between chrX CpG methylation and MMSE scores.

This repository contains three steps and seven Python scripts.

## 2. Repository structure

| Directory | File | Purpose |
|---|---|---|
| Root directory | `README.md` | Project description and execution order |
| Root directory | `environment.yml` | Records the actual `xci_0` environment |
| Root directory | `.gitignore` | Excludes local data and output files |
| `Step_1/` | `make_centenarian_sample_list.py` | Build sample list |
| `Step_1/` | `make_chrX_methylation_matrix_centenarian.py` | Build methylation matrix |
| `Step_1/` | `make_chrX_coverage_matrix_centenarian.py` | Build coverage matrix |
| `Step_2/` | `Coverage_and_CpG_completeness_exploration.py` | 比较 QC thresholds |
| `Step_2/` | `make_chrX_usable_matrix.py` | Apply QC，Export usable matrix |
| `Step_3/` | `run_primary_association.py` | Run Primary analysis |
| `Step_3/` | `sex_stratified_chrX.py` | Run Female、Male analyses |

`Input/` and `Output/` are local directories used during execution.

Actual sample data and generated outputs are not pushed to GitHub.

`environment.yml`

`.gitignore`

## 3. Software environment

Use `xci_0` Environment。

Main dependencies:

| Dependency | Purpose |
|---|---|
| Python | Run scripts |
| pandas | Read, join, and export tables |
| NumPy | Perform matrix operations |
| SciPy | Perform regression and statistical tests |
| Matplotlib | Generate QC plots |

Dependency versions are not recorded in the scripts.

The actual versions must be exported from the `xci_0` environment.

Run preparation:

1. Activate the `xci_0` environment.
2. Export the actual environment configuration.
3. Update the paths in each script.
4. Create the local working directories.
5. Run the scripts in order: Step_1 → Step_2 → Step_3.

```bash
conda activate xci_0
conda env export --no-builds > environment.yml

mkdir -p Step_1/Input Step_1/Output
mkdir -p Step_2/Output
mkdir -p Step_3/Input Step_3/Output
```

Run the commands above from the repository root directory.

The current scripts use hardcoded paths and do not accept CLI parameters.

Default working directory:

```text
/project/holstegelab/Share/yixuan/workflow
```

All script paths must point to the same working directory.

## 4. Input data

### 4.1 Raw Inputs

| Input | Required Columns / Content |
|---|---|
| `Step_1/Input/pb_samples_phenotypes_20240425.txt` | TSV; `pacbio_id`, `ID_GWAS`, `ID_100plus`, `PHENOTYPE`, `SEX` |
| `100plus_mmse.csv` | `studysubjectid`, `harmonized_contact_nr`, `mmse_score_imputed` |
| Per-sample `*.cpg.combined.bed` | Columns 1, 2, 4, and 6: chromosome, start, methylation, coverage |
| `Step_3/Input/final_model_samples_correct_age_cellprops.tsv` | `pacbio_id`, `SEX`, `mmse_score_imputed`, `PC1`–`PC10`, education, age, and six blood cell proportions |

Default MMSE path:

```text
/project/holstegelab/Share/yaran/sample/100plus_mmse.csv
```

Default BED path:

```text
/project/holstegelab/Share/yaran/epigenetics/phased/denovo/{pacbio_id}/{pacbio_id}.cpg.combined.bed
```

### 4.2 ID Mapping and Format Requirements

1. Use `pacbio_id` to align samples with matrices.
2. Join MMSE data using `ID_100plus ↔ studysubjectid`.
3. Standardize `SEX` values as `F` and `M`.
4. Ensure that `pacbio_id` is unique in each sample table.
5. Ensure that chrX `start` positions are unique within each sample.
6. Ensure that `studysubjectid` is unique among MMSE records with contact 0.

### 4.3 Step_3 covariates

| Model Variable | Input Column / Rule |
|---|---|
| MMSE | `mmse_score_imputed` |
| Education | Prefer `education`; use `verhage_education_level_numeric` if the `education` column is absent |
| Age | Prefer `age`; use `age_entrance_interview` if the `age` column is absent |
| Cell proportions | `CD8T`, `CD4T`, `NK`, `Bcell`, `Mono`, `Neu` |
| B-cell alias | Use `Bcel` if the `Bcell` column is absent |

If the input contains an `age` column, its values must be derived from `age_entrance_interview`.

The seven scripts do not build the Step_3 covariate table.

This file must be prepared in advance.

### 4.4 Matrix Format

Both raw matrices use the **sample × CpG** format.

The first three columns are:

```text
pacbio_id    ID_GWAS    ID_100plus
```

The remaining columns follow the format:

```text
chrX:{start}
```

CpG coordinates are taken directly from the BED `start` column.

The scripts do not apply a `+1` coordinate conversion.

**Methylation units must be consistent.**

Step_1 reads column 4 of the BED files without conversion.

Step_3 divides the methylation values by 100.

Therefore, the input must contain methylation percentages on a scale of 0–100.

## 5. Step_1 — Sample Selection and Matrix Construction

Run the scripts in the following order.

All outputs in this section are saved in `Step_1/Output/`.

### 5.1 `make_centenarian_sample_list.py`

**Purpose:** Build the centenarian sample list.

**Input:** Phenotype TSV, MMSE CSV, and BED file paths.

Processing steps:

1. Retain samples with `PHENOTYPE == "Centenarian"`.
2. Exclude samples whose `pacbio_id` contains `brain`.
3. Exclude samples with missing BED files.
4. Join MMSE records with `harmonized_contact_nr == "0"`.
5. Export the sample list.

**Output：** `centenarian_sample_list.tsv`。

Output columns：

```text
pacbio_id
ID_GWAS
ID_100plus
SEX
mmse_score_imputed
bed_path
```

This step retains samples with missing MMSE values.

### 5.2 `make_chrX_methylation_matrix_centenarian.py`

**Purpose:** Build the raw chrX methylation matrix.

**Input:** Sample list and the BED files specified by `bed_path`.

Processing steps:

1. Collect the union of chrX CpG sites across all samples.
2. Sort CpG sites by genomic coordinate.
3. Build the matrix using column 4 of the BED files.
4. Retain unobserved values as NaN.
5. Export the methylation matrix.

**Output：** `chrX_methylation_matrix.tsv`。

### 5.3 `make_chrX_coverage_matrix_centenarian.py`

**Purpose:** Build a coverage matrix aligned with the methylation matrix.

**Input:** Sample list, methylation matrix, and BED files.

Processing steps:

1. Preserve the sample and CpG order of the methylation matrix.
2. Build the coverage matrix using column 6 of the BED files.
3. Validate that coverage values are non-negative, finite integers.
4. Export the coverage matrix.

**Output:** `chrX_coverage_matrix.tsv`.

The rows and columns of the two matrices correspond exactly.

### 5.4 Execution Order

```bash
python Step_1/make_centenarian_sample_list.py
python Step_1/make_chrX_methylation_matrix_centenarian.py
python Step_1/make_chrX_coverage_matrix_centenarian.py
```

## 6. Step_2 — Coverage Quality Control

Both scripts first exclude sample `8783353`.

All outputs in this section are saved in `Step_2/Output/`.

### 6.1 `Coverage_and_CpG_completeness_exploration.py`

**Purpose:** Compare coverage and completeness thresholds.

**Input:** Coverage matrix and sample list.

Processing steps:

1. Exclude sample `8783353`.
2. Verify that 247 samples remain.
3. Compare two sets of coverage thresholds.
4. Compare completeness thresholds of 50%, 80%, and 90%.
5. Export QC plots and summary statistics.

Coverage settings:

| Mode | Female | Male |
|---|---|---|
| Same absolute threshold | 5, 10, 15, 20× | 5, 10, 15, 20× |
| Sex-adjusted threshold | 5, 10, 15, 20× | 2.5, 5, 7.5, 10× |

Completeness is calculated as the number of samples meeting the coverage threshold divided by the total of 247 samples.

NaN coverage values are treated as failing the threshold.

**Output:**

| File | Content |
|---|---|
| `QC_coverage_thresholds_247samples.png` | Number of CpGs meeting the coverage threshold per sample |
| `QC_completeness_thresholds_247samples.png` | Number of retained CpGs for each threshold combination |
| `QC_threshold_summary_247samples.tsv` | Summary statistics for 24 threshold combinations |

### 6.2 `make_chrX_usable_matrix.py`

**Purpose:** Apply the final QC criteria and export the matrix for association analysis.

**Input:** Methylation matrix, coverage matrix, and sample list.

Processing steps:

1. Exclude sample `8783353`.
2. Apply coverage thresholds of Female ≥10× and Male ≥5×.
3. Retain CpGs with coverage completeness ≥80%.
4. Mask methylation values as NaN where coverage is below the required threshold.
5. Mask methylation values as NaN where coverage is NaN.
6. Preserve existing NaN methylation values.
7. Export the final NaN matrix.

**Final Output:**

```text
chrX_usable_matrix_10X_80pct_NaN.tsv
```

With 247 samples, at least 198 samples must meet the coverage threshold.

Completeness is calculated based on coverage.

It is not equivalent to the proportion of non-missing methylation values.

### 6.3 Execution Order

```bash
python Step_2/Coverage_and_CpG_completeness_exploration.py
python Step_2/make_chrX_usable_matrix.py
```

## 7. Step_3 — Association analysis

Both scripts directly read:

```text
Step_2/Output/chrX_usable_matrix_10X_80pct_NaN.tsv
Step_3/Input/final_model_samples_correct_age_cellprops.tsv
```

Common processing steps:

1. Align the two inputs using `pacbio_id`.
2. Exclude samples with missing MMSE values or covariates.
3. Divide methylation percentages by 100.
4. Use CpG-specific complete cases for each CpG.
5. Fit the regression model.
6. Export the full association results and candidate results.

No imputation is performed on methylation values.

All outputs in this section are saved in `Step_3/Output/`.

### 7.1 `run_primary_association.py`

合并 Female 与 Male，Run Primary model：

```text
MMSE ~ CpG methylation + sex + age + education
       + PC1–PC10
       + CD8T + CD4T + NK + Bcell + Mono + Neu
```

The model includes an intercept.

Sex is encoded as Female = 0 and Male = 1.

The code requires 235 samples with complete MMSE and covariate data.

**Output:**

| File | Content |
|---|---|
| `primary_all_CpGs.tsv` | All successfully estimated CpGs |
| `Table_1_primary_suggestive_CpGs.tsv` | Candidate CpGs with `P < 1e-5`, sorted by P-value |

### 7.2 `sex_stratified_chrX.py`

分别在 Female 与 Male 中 Run：

```text
MMSE ~ CpG methylation + age + education
       + PC1–PC10
       + CD8T + CD4T + NK + Bcell + Mono + Neu
```

The model includes an intercept but excludes sex.

The code requires the input to contain 1,159,396 CpG columns.

**Output:**

| File | Content |
|---|---|
| `female_all_CpGs.tsv` | All CpGs in the female analysis |
| `male_all_CpGs.tsv` | All CpGs in the male analysis |
| `Table_female_suggestive_CpGs.tsv` | Female candidate CpGs with `P < 1e-5` |
| `Table_male_suggestive_CpGs.tsv` | Male candidate CpGs with `P < 1e-5` |

The full stratified result tables retain CpGs for which model estimates could not be obtained.

The statistical values for these CpGs are recorded as NaN.

The full primary result table excludes CpGs for which model estimates could not be obtained.

### 7.3 Result columns

| Column | Description |
|---|---|
| `CpG position` | CpG identifier |
| `beta` | Methylation regression coefficient adjusted for covariates; uses the 0–1 fraction scale |
| `SE` | Standard error of beta |
| `P-value` | Unadjusted P-value from a two-sided t-test |
| `N` | Actual number of complete-case samples for the CpG |

N varies across CpGs and is not equal to a fixed cohort size.

The values 235 and 1,159,396 are validation checks hardcoded in the scripts.

### 7.4 Execution Order

```bash
python Step_3/run_primary_association.py
python Step_3/sex_stratified_chrX.py
```

## 8. Analysis Settings

| Item | Final Analysis Setting |
|---|---|
| Sample QC | Exclude sample `8783353` |
| Coverage | Female ≥10×; Male ≥5× |
| Completeness | ≥80%; denominator of 247 samples |
| Methylation missing values | Use the NaN matrix; no imputation |
| Association | Use CpG-specific complete cases |
| Age | `age_entrance_interview` |
| Primary covariates | Sex, age, education, PC1–PC10, and six blood cell proportions |
| Stratified covariates | Same as the primary model, excluding sex |
| Candidate threshold | Unadjusted `P < 1 × 10⁻⁵` |
| Execution method | Load all CpGs at once and fit each CpG separately; no chunking |

The regression skips CpGs with insufficient samples or for which model estimates cannot be obtained.

This is not equivalent to predefined variance filtering.
