"""
Sex-stratified chrX methylation–MMSE association.

Female / Male:
MMSE ~ CpG + PC1-PC10 + education + age
       + CD8T + CD4T + NK + Bcell + Mono + Neu

Per-CpG complete-case regression; no methylation imputation.
Outputs: four TSV files.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import lstsq
from scipy.stats import t

# 1. Paths
BASE = Path("/project/holstegelab/Share/yixuan/workflow")
METH_FILE = BASE / "Step_2/Output/chrX_usable_matrix_10X_80pct_NaN.tsv"
COV_FILE = BASE / "Step_3/Input/final_model_samples_correct_age_cellprops.tsv"
OUT = BASE / "Step_3/Output"
OUT.mkdir(parents=True, exist_ok=True)

P_LIMIT = 1e-5
CELLS = ["CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"]
COVARS = [f"PC{i}" for i in range(1, 11)] + [
    "education", "age"
] + CELLS

# 2. Load and align samples
cov = pd.read_csv(
    COV_FILE, sep="\t", dtype={"pacbio_id": str}
).set_index("pacbio_id")

if "education" not in cov:
    cov["education"] = cov["verhage_education_level_numeric"]
if "age" not in cov:
    cov["age"] = cov["age_entrance_interview"]
if "Bcell" not in cov and "Bcel" in cov:
    cov["Bcell"] = cov["Bcel"]

cov["SEX"] = cov["SEX"].astype(str).str.strip().str.upper().str[0]

needed = ["mmse_score_imputed"] + COVARS
cov[needed] = cov[needed].apply(pd.to_numeric, errors="coerce")
cov = cov.dropna(subset=needed)
cov = cov[cov["SEX"].isin(["F", "M"])]

meth = pd.read_csv(
    METH_FILE, sep="\t", dtype={"pacbio_id": str}
).set_index("pacbio_id")

cov = cov.loc[cov.index.isin(meth.index)]
cpgs = meth.columns[meth.columns.str.startswith("chrX:")]
M = meth.loc[cov.index, cpgs].to_numpy(dtype=np.float64)
del meth

assert len(cpgs) == 1_159_396, "Unexpected CpG count"

y_all = cov["mmse_score_imputed"].to_numpy(float)
Z_all = np.column_stack([
    np.ones(len(cov)), cov[COVARS].to_numpy(float)
])

# 3. Female and Male regression
for name, sex in [("female", "F"), ("male", "M")]:
    selected = (cov["SEX"] == sex).to_numpy()
    y = y_all[selected]
    Z = Z_all[selected]
    X = M[selected]

    beta = np.full(len(cpgs), np.nan)
    se = np.full(len(cpgs), np.nan)
    pval = np.full(len(cpgs), np.nan)
    nval = np.zeros(len(cpgs), dtype=int)

    for i in range(len(cpgs)):
        x = X[:, i] / 100.0
        valid = np.isfinite(x)
        n = int(valid.sum())
        nval[i] = n

        if n <= Z.shape[1] + 1:
            continue

        Zv = Z[valid]
        data = np.column_stack((x[valid], y[valid]))

        coef, _, rank, _ = lstsq(
            Zv, data, check_finite=False
        )
        rx, ry = (data - Zv @ coef).T

        df = n - rank - 1
        denom = rx @ rx

        if df <= 0 or denom <= 1e-14:
            continue

        b = (rx @ ry) / denom
        resid = ry - b * rx
        s = np.sqrt((resid @ resid) / df / denom)

        if not np.isfinite(s) or s <= 0:
            continue

        beta[i] = b
        se[i] = s
        pval[i] = 2 * t.sf(abs(b / s), df)

    del X

    # 4. Save all CpG results
    results = pd.DataFrame({
        "CpG position": cpgs,
        "beta": beta,
        "SE": se,
        "P-value": pval,
        "N": nval
    })

    results.to_csv(
        OUT / f"{name}_all_CpGs.tsv",
        sep="\t",
        index=False
    )

    # 5. Save suggestive CpGs
    hits = results.loc[
        results["P-value"] < P_LIMIT
    ].sort_values("P-value")

    hits.to_csv(
        OUT / f"Table_{name}_suggestive_CpGs.tsv",
        sep="\t",
        index=False
    )

    print(
        f"{name}: samples={len(y)}, "
        f"suggestive CpGs={len(hits)}",
        flush=True
    )

print("Done.")
