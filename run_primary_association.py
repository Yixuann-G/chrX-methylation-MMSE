from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import lstsq
from scipy.stats import t

# Paths
BASE = Path("/project/holstegelab/Share/yixuan/workflow")
METH_FILE = BASE / "Step_2/Output/chrX_usable_matrix_10X_80pct_NaN.tsv"
COV_FILE = BASE / "Step_3/Input/final_model_samples_correct_age_cellprops.tsv"
OUT = BASE / "Step_3/Output"
OUT.mkdir(parents=True, exist_ok=True)

P_THRESHOLD = 1e-5

# 1. Load and align samples
meth = pd.read_csv(METH_FILE, sep="\t")
cov = pd.read_csv(COV_FILE, sep="\t", dtype={"pacbio_id": str})

meth["pacbio_id"] = meth["pacbio_id"].astype(str)
cov = cov.set_index("pacbio_id")

if not meth["pacbio_id"].is_unique or not cov.index.is_unique:
    raise ValueError("Duplicate pacbio_id detected")

common = meth["pacbio_id"].isin(cov.index)
meth = meth.loc[common].reset_index(drop=True)
cov = cov.loc[meth["pacbio_id"]].copy()

if "education" not in cov:
    cov["education"] = cov["verhage_education_level_numeric"]

if "age" not in cov:
    cov["age"] = cov["age_entrance_interview"]

if "Bcell" not in cov and "Bcel" in cov:
    cov["Bcell"] = cov["Bcel"]

sex = cov["SEX"].astype(str).str.upper()
cov["sex_numeric"] = sex.map({
    "F": 0, "FEMALE": 0, "0": 0,
    "M": 1, "MALE": 1, "1": 1
})

variables = (
    ["sex_numeric"]
    + [f"PC{i}" for i in range(1, 11)]
    + ["education", "age"]
    + ["CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"]
)

covariates = cov[variables].apply(pd.to_numeric, errors="coerce")
y = pd.to_numeric(cov["mmse_score_imputed"], errors="coerce")

valid = y.notna() & covariates.notna().all(axis=1)
meth = meth.loc[valid.to_numpy()].reset_index(drop=True)
y = y.loc[valid].to_numpy(dtype=float)
C = covariates.loc[valid].to_numpy(dtype=float)
C = np.column_stack([np.ones(len(y)), C])

if len(y) != 235:
    raise ValueError(f"Expected 235 samples, got {len(y)}")

# 2. Load CpG methylation matrix
cpgs = [c for c in meth.columns if c.startswith("chrX:")]
X = np.asfortranarray(meth[cpgs].to_numpy(dtype=np.float32))
X /= 100.0
del meth

# 3. CpG-wise regression with complete cases
results = []

for j, site in enumerate(cpgs):
    x = X[:, j].astype(float)
    keep = np.isfinite(x)
    n = int(keep.sum())

    if n <= C.shape[1] + 1:
        continue

    A = C[keep]
    yy = y[keep]
    xx = x[keep]

    fitted, _, rank, _ = lstsq(
        A, np.column_stack([yy, xx]),
        lapack_driver="gelsy"
    )

    ry = yy - A @ fitted[:, 0]
    rx = xx - A @ fitted[:, 1]

    ssx = rx @ rx
    df = n - rank - 1

    if df <= 0 or ssx <= 1e-12:
        continue

    beta = (rx @ ry) / ssx
    residual = ry - beta * rx
    se = np.sqrt((residual @ residual) / df / ssx)

    if se <= 0 or not np.isfinite(se):
        continue

    p = 2 * t.sf(abs(beta / se), df)
    results.append((site, beta, se, p, n))

del X

# 4. Save all CpG results
res = pd.DataFrame(
    results,
    columns=["CpG position", "beta", "SE", "P-value", "N"]
)

res.to_csv(
    OUT / "primary_all_CpGs.tsv",
    sep="\t",
    index=False
)

# 5. Save suggestive CpGs
hits = res.loc[res["P-value"] < P_THRESHOLD].copy()
hits = hits.sort_values("P-value")

hits.to_csv(
    OUT / "Table_1_primary_suggestive_CpGs.tsv",
    sep="\t",
    index=False
)

print(f"Samples: {len(y)}")
print(f"CpGs tested: {len(res)}")
print(f"Suggestive CpGs: {len(hits)}")
