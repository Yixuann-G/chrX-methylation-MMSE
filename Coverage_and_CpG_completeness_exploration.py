from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import StrMethodFormatter

# ============================================================
# 1. Paths and parameters
# ============================================================

BASE = Path("/project/holstegelab/Share/yixuan/workflow")

COVERAGE = BASE / "Step_1/Output/chrX_coverage_matrix.tsv"
SAMPLES = BASE / "Step_1/Output/centenarian_sample_list.tsv"

OUT_DIR = BASE / "Step_2/Output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUTLIER_ID = "8783353"

THRESHOLDS = [5, 10, 15, 20]
COMPLETENESS = [50, 80, 90]

COLORS = {
    "F": "#B65C85",
    "M": "#487EAD"
}

MODES = [
    "Same absolute threshold",
    "Sex-adjusted threshold"
]

EXPECTED_SAMPLES = 247
EXPECTED_CPGS = 1159396


# ============================================================
# 2. Load sample IDs and sex
# ============================================================

print("Loading sample metadata...", flush=True)

ids = pd.read_csv(
    COVERAGE,
    sep="\t",
    usecols=["pacbio_id"],
    dtype=str
)["pacbio_id"]

samples = pd.read_csv(
    SAMPLES,
    sep="\t",
    usecols=["pacbio_id", "SEX"],
    dtype=str
)

for values in (ids, samples["pacbio_id"]):
    if values.isna().any() or values.duplicated().any():
        raise ValueError("Missing or duplicated pacbio_id")

sex_series = samples.set_index("pacbio_id")["SEX"]
sex_all = ids.map(sex_series).to_numpy()

if not np.isin(sex_all, ["F", "M"]).all():
    raise ValueError("Missing or invalid SEX")

# IMPORTANT: Remove coverage outlier BEFORE QC calculations

keep_mask = ids.to_numpy() != OUTLIER_ID

if np.count_nonzero(~keep_mask) != 1:
    raise ValueError(
        f"Expected exactly one outlier: {OUTLIER_ID}"
    )

sex = sex_all[keep_mask]

print(
    f"Raw samples: {len(ids)}; "
    f"After outlier removal: {len(sex)}",
    flush=True
)


# ============================================================
# 3. Load coverage matrix
# ============================================================

print("Reading coverage matrix...", flush=True)

cov = pd.read_csv(
    COVERAGE,
    sep="\t",
    usecols=lambda c: c.startswith("chrX:"),
    dtype=np.float32
).to_numpy()

if cov.shape[0] != len(ids):
    raise ValueError("Coverage and metadata rows do not match")

# Remove outlier row
cov = cov[keep_mask, :]

n_samples, n_cpg = cov.shape

if n_samples != EXPECTED_SAMPLES:
    raise ValueError(
        f"Expected {EXPECTED_SAMPLES} samples, got {n_samples}"
    )

if n_cpg == 0:
    raise ValueError("No CpGs found")

if np.isinf(cov).any() or np.any(cov < 0):
    raise ValueError("Invalid coverage values")

print(
    f"Final samples: {n_samples}; "
    f"Raw CpGs: {n_cpg:,}",
    flush=True
)


# ============================================================
# 4. Calculate QC statistics
# ============================================================

sample_counts = {}
retained_counts = {}

for mode in MODES:

    sample_counts[mode] = []
    retained_counts[mode] = []

    for t in THRESHOLDS:

        if mode == "Sex-adjusted threshold":
            thresholds = np.where(
                sex == "M",
                t / 2,
                t
            )
        else:
            thresholds = np.full(n_samples, t)

        # Missing coverage does not pass
        passed = cov >= thresholds[:, None]

        # Number of passing CpGs per sample
        sample_counts[mode].append(
            np.count_nonzero(passed, axis=1)
        )

        # Number of passing samples per CpG
        passing_samples = np.count_nonzero(
            passed,
            axis=0
        )

        retained_counts[mode].append([
            np.count_nonzero(
                passing_samples * 100 >= p * n_samples
            )
            for p in COMPLETENESS
        ])

        del passed, passing_samples

    retained_counts[mode] = np.array(
        retained_counts[mode]
    )

del cov


# ============================================================
# 5. Verify final QC count
# ============================================================

mode = "Sex-adjusted threshold"

threshold_index = THRESHOLDS.index(10)
completeness_index = COMPLETENESS.index(80)

final_count = retained_counts[mode][
    threshold_index,
    completeness_index
]

minimum_samples = int(np.ceil(
    n_samples * 0.80
))

print("\n========== FINAL QC ==========", flush=True)
print("Female coverage: >=10X")
print("Male coverage: >=5X")
print("Completeness: >=80%")
print("Samples:", n_samples)
print("Minimum passing samples:", minimum_samples)
print("Retained CpGs:", f"{final_count:,}")
print("Expected CpGs:", f"{EXPECTED_CPGS:,}")

if final_count != EXPECTED_CPGS:
    print(
        "WARNING: QC count does not match expected result.",
        flush=True
    )
else:
    print("QC count matches expected result.", flush=True)


# ============================================================
# 6. Figure 1: Passing CpGs per sample
# ============================================================

fig, axes = plt.subplots(
    1, 2,
    figsize=(14, 6),
    sharey=True
)

x = np.arange(len(THRESHOLDS))

jitter = np.random.default_rng(42).uniform(
    -0.06,
    0.06,
    n_samples
)

for ax, mode in zip(axes, MODES):

    for j, t in enumerate(THRESHOLDS):

        counts = sample_counts[mode][j]

        for label, offset in [
            ("F", -0.18),
            ("M", 0.18)
        ]:

            group = sex == label

            if not group.any():
                continue

            position = j + offset

            ax.boxplot(
                counts[group],
                positions=[position],
                widths=0.26,
                patch_artist=True,
                showfliers=False,
                manage_ticks=False,
                boxprops=dict(
                    facecolor=COLORS[label],
                    alpha=0.25
                ),
                medianprops=dict(color="black")
            )

            ax.scatter(
                position + jitter[group],
                counts[group],
                color=COLORS[label],
                s=12,
                alpha=0.6,
                zorder=3
            )

    labels = [
        f"F≥{t} / M≥{t if mode == MODES[0] else t / 2:g}"
        for t in THRESHOLDS
    ]

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_title(mode)
    ax.set_xlabel("Coverage threshold (×)")
    ax.set_ylim(bottom=0)

    ax.yaxis.set_major_formatter(
        StrMethodFormatter("{x:,.0f}")
    )

    ax.grid(axis="y", alpha=0.2)

axes[0].set_ylabel(
    "Number of passing CpGs per sample"
)

axes[1].legend(handles=[
    Patch(
        facecolor=COLORS["F"],
        label="Female"
    ),
    Patch(
        facecolor=COLORS["M"],
        label="Male"
    )
])

fig.suptitle(
    "Chromosome X coverage QC (247 samples)",
    fontsize=15
)

fig.tight_layout()

fig.savefig(
    OUT_DIR / "QC_coverage_thresholds_247samples.png",
    dpi=300,
    facecolor="white"
)

plt.close(fig)


# ============================================================
# 7. Figure 2: Retained CpGs by completeness
# ============================================================

fig, axes = plt.subplots(
    1, 2,
    figsize=(15, 6),
    sharey=True
)

width = 0.24

BAR_COLORS = {
    50: "#1f77b4",
    80: "#ff7f0e",
    90: "#2ca02c"
}

for ax, mode in zip(axes, MODES):

    for k, percentage in enumerate(COMPLETENESS):

        values = retained_counts[mode][:, k]

        bars = ax.bar(
            x + (k - 1) * width,
            values,
            width,
            color=BAR_COLORS[percentage],
            label=f"≥{percentage}% samples"
        )

        for j, (bar, value) in enumerate(zip(bars, values)):

            # Highlight selected QC combination
            is_selected = (
                mode == "Sex-adjusted threshold"
                and THRESHOLDS[j] == 10
                and percentage == 80
            )

            if is_selected:
                bar.set_edgecolor("black")
                bar.set_linewidth(2.5)

            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value,
                f"{value:,}",
                ha="center",
                va="bottom",
                rotation=90,
                fontsize=8,
                fontweight="bold" if is_selected else "normal"
            )

    labels = [
        f"F≥{t} / M≥{t if mode == MODES[0] else t / 2:g}"
        for t in THRESHOLDS
    ]

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_title(mode)
    ax.set_xlabel("Coverage threshold (×)")
    ax.set_ylim(0, n_cpg * 1.25)

    ax.yaxis.set_major_formatter(
        StrMethodFormatter("{x:,.0f}")
    )

    ax.grid(axis="y", alpha=0.2)
    ax.legend()

axes[0].set_ylabel("Number of CpGs retained")

fig.suptitle(
    "ChrX CpG retention after outlier removal (N=247)",
    fontsize=15
)

fig.tight_layout()

fig.savefig(
    OUT_DIR / "QC_completeness_thresholds_247samples.png",
    dpi=300,
    facecolor="white"
)

plt.close(fig)


# ============================================================
# 8. Export QC summary
# ============================================================

summary_rows = []

for mode in MODES:
    for j, t in enumerate(THRESHOLDS):
        for k, percentage in enumerate(COMPLETENESS):

            summary_rows.append({
                "Mode": mode,
                "Female_threshold": t,
                "Male_threshold": (
                    t if mode == MODES[0] else t / 2
                ),
                "Completeness_percent": percentage,
                "Samples": n_samples,
                "Retained_CpGs": int(
                    retained_counts[mode][j, k]
                )
            })

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    OUT_DIR / "QC_threshold_summary_247samples.tsv",
    sep="\t",
    index=False
)

print("\nSaved figures and QC summary to:", OUT_DIR)
print("DONE", flush=True)