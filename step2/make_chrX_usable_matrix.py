#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gc
import os
import sys
import numpy as np
import pandas as pd

sys.stdout.reconfigure(line_buffering=True)


# =========================
# Input / output paths
# =========================

BASE = "/project/holstegelab/Share/yixuan/workflow"

METHYLATION_MATRIX = f"{BASE}/Step_1/Output/chrX_methylation_matrix.tsv"
COVERAGE_MATRIX = f"{BASE}/Step_1/Output/chrX_coverage_matrix.tsv"

PHENOTYPE_FILE = f"{BASE}/Step_1/Output/centenarian_sample_list.tsv"

OUT_NAN_MATRIX = f"{BASE}/Step_2/Output/chrX_usable_matrix_10X_80pct_NaN.tsv"
OUT_IMPUTED_MATRIX = f"{BASE}/Step_2/Output/chrX_usable_matrix_10X_80pct_imputed.tsv"


# =========================
# Fixed parameters
# =========================

META_COLS = ["pacbio_id", "ID_GWAS", "ID_100plus"]

OUTLIER_PACBIO_ID = "8783353"

FEMALE_THRESHOLD = 10.0
MALE_THRESHOLD = 5.0

PHENO_ID_COL = "pacbio_id"
PHENO_SEX_COL = "SEX"

MALE_VALUES = {"M", "Male", "male", "MALE"}
FEMALE_VALUES = {"F", "Female", "female", "FEMALE"}

FLOAT_FORMAT = None


# =========================
# Helper functions
# =========================

def select_keep_cpgs_from_coverage(coverage_path, phenotype_path):
    """Retain CpGs passing sex-specific coverage in >=80% of QC samples."""
    print(f"[READ] coverage matrix: {coverage_path}")

    # Read sample IDs in the same order as the coverage matrix.
    ids = pd.read_csv(
        coverage_path,
        sep="\t",
        usecols=["pacbio_id"],
        dtype=str,
    )["pacbio_id"]

    if ids.isna().any() or ids.duplicated().any():
        raise ValueError("Missing or duplicated pacbio_id in coverage matrix")

    # Exclude the predefined coverage outlier before CpG completeness QC.
    keep_samples = ids.ne(OUTLIER_PACBIO_ID).to_numpy()
    if np.count_nonzero(~keep_samples) != 1:
        raise ValueError("Expected exactly one predefined coverage outlier")

    filtered_ids = ids.loc[keep_samples].tolist()
    sex_map = load_sex_map(phenotype_path)
    thresholds = make_threshold_vector(filtered_ids, sex_map)

    # Only load CpG coverage columns. Keep their original column order.
    cov_df = pd.read_csv(
        coverage_path,
        sep="\t",
        usecols=lambda col: col.startswith("chrX:"),
        dtype=np.float32,
        na_values=["NA", "NaN", "nan", ""],
        low_memory=False,
    )

    if cov_df.shape[0] != len(ids):
        raise ValueError("Coverage matrix and sample ID row counts differ")
    if cov_df.shape[1] == 0:
        raise ValueError("Coverage matrix contains no chrX CpGs")

    cov = cov_df.to_numpy(dtype=np.float32, copy=False)[keep_samples]

    # Female >=10X, male >=5X; NaN never passes coverage QC.
    passing_counts = np.count_nonzero(cov >= thresholds[:, None], axis=0)
    minimum_samples = int(np.ceil(0.80 * len(filtered_ids)))
    keep_mask = passing_counts >= minimum_samples
    keep_cpgs = cov_df.columns[keep_mask].tolist()

    print(f"[QC] Samples after outlier removal: {len(filtered_ids)}")
    print(f"[QC] Minimum samples passing coverage: {minimum_samples}")
    print(f"[QC] Retained CpGs: {len(keep_cpgs):,}")

    if not keep_cpgs:
        raise ValueError("No CpGs pass the coverage and completeness thresholds")

    del cov_df, cov, passing_counts
    gc.collect()
    return keep_cpgs


def read_header(path):
    with open(path, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
    return header


def read_selected_matrix(path, keep_cpgs, label):
    print(f"[READ] selected {label} matrix:")
    print(f"       {path}")

    header = read_header(path)

    if header[:len(META_COLS)] != META_COLS:
        raise ValueError(
            f"{label} matrix metadata columns are wrong. "
            f"Expected {META_COLS}, found {header[:len(META_COLS)]}"
        )

    header_set = set(header)
    missing_cpgs = [cpg for cpg in keep_cpgs if cpg not in header_set]

    if missing_cpgs:
        raise ValueError(
            f"{label} matrix is missing {len(missing_cpgs)} kept CpGs. "
            f"Examples: {missing_cpgs[:5]}"
        )

    usecols = META_COLS + keep_cpgs

    dtype_map = {col: str for col in META_COLS}
    dtype_map.update({cpg: "float32" for cpg in keep_cpgs})

    df = pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        dtype=dtype_map,
        na_values=["NA", "NaN", "nan", ""],
        keep_default_na=True,
        low_memory=False,
    )

    df = df[usecols]

    print(f"[INFO] {label} shape: {df.shape}")
    print(f"[INFO] {label} samples: {df.shape[0]}")
    print(f"[INFO] {label} CpGs: {df.shape[1] - len(META_COLS)}")

    return df


def check_alignment(meth_df, cov_df, label):
    print(f"[CHECK] alignment: {label}")

    if meth_df.shape != cov_df.shape:
        raise ValueError(
            f"Shape mismatch after {label}: "
            f"methylation {meth_df.shape}, coverage {cov_df.shape}"
        )

    if meth_df.columns.tolist() != cov_df.columns.tolist():
        raise ValueError(f"Column mismatch after {label}.")

    for col in META_COLS:
        meth_values = meth_df[col].astype(str).tolist()
        cov_values = cov_df[col].astype(str).tolist()

        if meth_values != cov_values:
            raise ValueError(f"Sample order mismatch in column {col} after {label}.")

    print("[CHECK] alignment OK")


def remove_outlier(meth_df, cov_df):
    print(f"[FILTER] remove outlier pacbio_id == {OUTLIER_PACBIO_ID}")

    mask_keep = meth_df["pacbio_id"].astype(str) != OUTLIER_PACBIO_ID
    n_removed = int((~mask_keep).sum())

    print(f"[INFO] removed rows: {n_removed}")

    meth_df = meth_df.loc[mask_keep].reset_index(drop=True)
    cov_df = cov_df.loc[mask_keep].reset_index(drop=True)

    print(f"[INFO] samples after removal: {meth_df.shape[0]}")

    return meth_df, cov_df


def load_sex_map(path):
    print(f"[READ] phenotype file:")
    print(f"       {path}")

    pheno = pd.read_csv(path, sep="\t", dtype={PHENO_ID_COL: str})

    if PHENO_ID_COL not in pheno.columns:
        raise KeyError(f"Phenotype file missing column: {PHENO_ID_COL}")

    if PHENO_SEX_COL not in pheno.columns:
        raise KeyError(f"Phenotype file missing column: {PHENO_SEX_COL}")

    sex_map = {}

    for pid, sex in zip(pheno[PHENO_ID_COL].astype(str), pheno[PHENO_SEX_COL].astype(str)):
        sex = sex.strip()

        if sex in MALE_VALUES:
            sex_map[pid] = "male"
        elif sex in FEMALE_VALUES:
            sex_map[pid] = "female"
        else:
            raise ValueError(f"Unknown SEX value for pacbio_id {pid}: {sex}")

    print("[INFO] phenotype SEX counts:")
    print(pheno[PHENO_SEX_COL].astype(str).str.strip().value_counts(dropna=False))

    return sex_map


def make_threshold_vector(pacbio_ids, sex_map):
    thresholds = []

    missing = []

    for pid in pacbio_ids:
        pid = str(pid)

        if pid not in sex_map:
            missing.append(pid)
            continue

        if sex_map[pid] == "male":
            thresholds.append(MALE_THRESHOLD)
        elif sex_map[pid] == "female":
            thresholds.append(FEMALE_THRESHOLD)
        else:
            raise ValueError(f"Invalid sex label for pacbio_id {pid}: {sex_map[pid]}")

    if missing:
        raise ValueError(
            f"{len(missing)} samples are missing SEX in phenotype file. "
            f"Examples: {missing[:5]}"
        )

    thresholds = np.array(thresholds, dtype=np.float32)

    print(f"[INFO] threshold vector length: {len(thresholds)}")
    print(f"[INFO] female threshold: {FEMALE_THRESHOLD}")
    print(f"[INFO] male threshold: {MALE_THRESHOLD}")
    print(f"[INFO] threshold counts:")
    print(pd.Series(thresholds).value_counts().sort_index())

    return thresholds


def write_matrix(path, meta_df, arr, cpg_cols):
    print(f"[WRITE] {path}")

    out_df = pd.concat(
        [
            meta_df.reset_index(drop=True),
            pd.DataFrame(arr, columns=cpg_cols),
        ],
        axis=1,
    )

    out_df.to_csv(
        path,
        sep="\t",
        index=False,
        na_rep="NaN",
        float_format=FLOAT_FORMAT,
    )

    del out_df
    gc.collect()


def impute_column_means_inplace(arr):
    print("[IMPUTE] column-mean imputation")

    col_means = np.nanmean(arr, axis=0).astype(np.float32)

    bad_cols = np.isnan(col_means)
    n_bad = int(bad_cols.sum())

    if n_bad > 0:
        bad_idx = np.where(bad_cols)[0][:10].tolist()
        raise ValueError(
            f"{n_bad} CpG columns have all NaN values. "
            f"Example column indexes: {bad_idx}"
        )

    n_imputed = 0

    for j in range(arr.shape[1]):
        missing = np.isnan(arr[:, j])

        if missing.any():
            n_imputed += int(missing.sum())
            arr[missing, j] = col_means[j]

    print(f"[IMPUTE] imputed cells: {n_imputed}")

    if np.isnan(arr).any():
        raise ValueError("NaN still exists after imputation.")

    return arr


# =========================
# Main
# =========================

def main():
    print("=" * 80)
    print("Build chrX 10X/80pct NaN matrix and imputed matrix")
    print("=" * 80)

    # 1. Determine retained CpGs directly from the Step_1 coverage matrix
    keep_cpgs = select_keep_cpgs_from_coverage(
        COVERAGE_MATRIX, PHENOTYPE_FILE
    )

    # 2. Read META_COLS + kept CpGs from methylation matrix
    meth_df = read_selected_matrix(
        METHYLATION_MATRIX,
        keep_cpgs,
        label="methylation",
    )

    # 3. Read META_COLS + kept CpGs from coverage matrix
    cov_df = read_selected_matrix(
        COVERAGE_MATRIX,
        keep_cpgs,
        label="coverage",
    )

    # 4. Check alignment before outlier removal
    check_alignment(meth_df, cov_df, label="before outlier removal")

    # 5. Remove outlier sample
    meth_df, cov_df = remove_outlier(meth_df, cov_df)

    # 6. Check alignment again after outlier removal
    check_alignment(meth_df, cov_df, label="after outlier removal")

    # 7. Build threshold vector from SEX
    sex_map = load_sex_map(PHENOTYPE_FILE)
    threshold_vector = make_threshold_vector(
        meth_df["pacbio_id"].astype(str).tolist(),
        sex_map,
    )

    # Keep metadata before converting blocks
    meta_df = meth_df[META_COLS].copy()
    cpg_cols = keep_cpgs

    # 8. Convert selected methylation block to numpy array
    print("[ARRAY] convert methylation block to numpy")
    meth = meth_df[cpg_cols].to_numpy(dtype=np.float32, copy=True)

    # 9. Convert selected coverage block to numpy array
    print("[ARRAY] convert coverage block to numpy")
    cov = cov_df[cpg_cols].to_numpy(dtype=np.float32, copy=True)

    print(f"[INFO] methylation array shape: {meth.shape}")
    print(f"[INFO] coverage array shape: {cov.shape}")

    del meth_df, cov_df
    gc.collect()

    # 10. Mask methylation by coverage
    print("[MASK] coverage below sex-specific threshold -> methylation NaN")

    bad_coverage = np.isnan(cov) | (cov < threshold_vector[:, None])
    n_bad = int(bad_coverage.sum())

    meth[bad_coverage] = np.nan

    del cov, bad_coverage
    gc.collect()

    nan_count = int(np.isnan(meth).sum())
    total_count = meth.size
    nan_fraction = nan_count / total_count

    print(f"[INFO] NaN cells after masking: {nan_count}")
    print(f"[INFO] total methylation cells: {total_count}")
    print(f"[INFO] NaN fraction: {nan_fraction:.6f}")

    # 11. Write NaN matrix first
    write_matrix(
        OUT_NAN_MATRIX,
        meta_df,
        meth,
        cpg_cols,
    )

    print("[DONE] NaN matrix written successfully")

    # 12. Use NaN matrix to generate imputed matrix
    impute_column_means_inplace(meth)

    # 13. Write imputed matrix
    write_matrix(
        OUT_IMPUTED_MATRIX,
        meta_df,
        meth,
        cpg_cols,
    )

    print("[DONE] imputed matrix written successfully")

    print("=" * 80)
    print("Finished")
    print(f"NaN matrix:     {OUT_NAN_MATRIX}")
    print(f"Imputed matrix: {OUT_IMPUTED_MATRIX}")
    print("=" * 80)


if __name__ == "__main__":
    main()
