from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("/project/holstegelab/Share/yixuan/workflow/Step_1/Output")
SAMPLE_LIST = BASE / "centenarian_sample_list.tsv"
OUT_MATRIX = BASE / "chrX_methylation_matrix.tsv"

ID_COLS = ["pacbio_id", "ID_GWAS", "ID_100plus"]
samples = pd.read_csv(SAMPLE_LIST, sep="\t", dtype=str)
samples[ID_COLS + ["bed_path"]]  # Check required columns

# 1. Collect the union of chrX coordinates
coords = set()

for path in samples["bed_path"]:
    bed = pd.read_csv(
        path, sep="\t", header=None,
        usecols=[0, 1], names=["chrom", "start"],
        dtype={"chrom": str, "start": np.int32},
    )
    coords.update(bed.loc[bed["chrom"].eq("chrX"), "start"])

coords = np.array(sorted(coords), dtype=np.int32)
coord_index = pd.Index(coords)


# 2. Fill methylation values from BED column 4
matrix = np.full((len(samples), len(coords)), np.nan, dtype=np.float32)

for i, row in samples.iterrows():
    bed = pd.read_csv(
        row["bed_path"], sep="\t", header=None,
        usecols=[0, 1, 3], names=["chrom", "start", "meth_score"],
        dtype={"chrom": str, "start": np.int32, "meth_score": np.float32},
    )
    chrX = bed.loc[bed["chrom"].eq("chrX"), ["start", "meth_score"]]

    if chrX["start"].duplicated().any():
        print(f"WARNING: {row['pacbio_id']}: averaging duplicate starts", flush=True)
        chrX = chrX.groupby("start", as_index=False)["meth_score"].mean()

    indices = coord_index.get_indexer(chrX["start"])
    valid = indices >= 0

    if not valid.all():
        print(f"WARNING: {row['pacbio_id']}: skipping new coordinates", flush=True)

    matrix[i, indices[valid]] = chrX["meth_score"].to_numpy()[valid]
    print(f"Methylation: {i + 1}/{len(samples)}", flush=True)

# 3. Export metadata and matrix
out = pd.concat(
    [
        samples[ID_COLS],
        pd.DataFrame(matrix, columns=[f"chrX:{c}" for c in coords]),
    ],
    axis=1,
)
out.to_csv(OUT_MATRIX, sep="\t", index=False)

print(f"Matrix shape: {matrix.shape}")
print(f"NaN count: {np.isnan(matrix).sum():,}")
print(f"Saved: {OUT_MATRIX}")
