from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("/project/holstegelab/Share/yixuan/workflow/Step_1/Output")
SAMPLE_LIST = BASE / "centenarian_sample_list.tsv"
METH_MATRIX = BASE / "chrX_methylation_matrix.tsv"
OUT_MATRIX = BASE / "chrX_coverage_matrix.tsv"

ID_COLS = ["pacbio_id", "ID_GWAS", "ID_100plus"]

# 1. Read the methylation matrix structure
columns = pd.read_csv(METH_MATRIX, sep="\t", nrows=0).columns.tolist()
metadata = pd.read_csv(
    METH_MATRIX, sep="\t", usecols=ID_COLS, dtype=str
)[ID_COLS]
cpg_cols = [c for c in columns if c not in ID_COLS]

if not cpg_cols or any(
    not c.startswith("chrX:") or not c[5:].isdigit() for c in cpg_cols
):
    raise ValueError("Missing or invalid CpG columns")

coords = pd.Index([int(c.split(":")[1]) for c in cpg_cols])
if coords.has_duplicates:
    raise ValueError("Duplicated CpG coordinates")

# 2. Map sample IDs to BED paths
samples = pd.read_csv(SAMPLE_LIST, sep="\t", dtype=str)

for frame in (metadata, samples):
    if frame["pacbio_id"].isna().any() or frame["pacbio_id"].duplicated().any():
        raise ValueError("Missing or duplicated pacbio_id")

if samples["bed_path"].isna().any():
    raise ValueError("Missing BED path")

bed_paths = samples.set_index("pacbio_id")["bed_path"]
missing = metadata.loc[
    ~metadata["pacbio_id"].isin(bed_paths.index), "pacbio_id"
]
if len(missing):
    raise ValueError(f"Samples absent from sample list: {missing.tolist()}")

# 3. Fill coverage values from BED column 6
matrix = np.full((len(metadata), len(coords)), np.nan, dtype=np.float32)

for i, pid in enumerate(metadata["pacbio_id"]):
    path = bed_paths.loc[pid]
    bed = pd.read_csv(
        path, sep="\t", header=None,
        usecols=[0, 1, 5], names=["chrom", "start", "coverage"],
        dtype={"chrom": str, "start": np.int64, "coverage": np.float32},
    )

    if bed.empty or bed.isna().any().any() or (bed["coverage"] < 0).any():
        raise ValueError(f"{pid}: empty BED, missing fields or negative coverage")

    chrX = bed.loc[bed["chrom"].eq("chrX")]
    if chrX.empty or chrX["start"].duplicated().any():
        raise ValueError(f"{pid}: no chrX rows or duplicated chrX starts")

    indices = coords.get_indexer(chrX["start"])
    matched = indices >= 0
    values = chrX.loc[matched, "coverage"].to_numpy()

    if not len(values):
        raise ValueError(f"{pid}: no matching CpGs")
    if not np.isfinite(values).all() or not np.allclose(values, np.round(values)):
        raise ValueError(f"{pid}: coverage must contain finite integer counts")

    matrix[i, indices[matched]] = values
    fill_rate = len(values) / len(coords)

    if fill_rate < 0.95:
        print(f"WARNING: {pid}: fill rate {fill_rate:.2%}", flush=True)

    print(f"Coverage: {i + 1}/{len(metadata)}", flush=True)

# 4. Export in exactly the methylation matrix column order
out = pd.concat(
    [metadata, pd.DataFrame(matrix, columns=cpg_cols)],
    axis=1,
)[columns]

out.to_csv(OUT_MATRIX, sep="\t", index=False)

print(f"Matrix shape: {matrix.shape}")
print(f"Saved: {OUT_MATRIX}")