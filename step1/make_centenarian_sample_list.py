#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from pathlib import Path

import pandas as pd


sys.stdout.reconfigure(line_buffering=True)


PHENO_FILE = Path(
    "/project/holstegelab/Share/yixuan/workflow/Step_1/Input/"
    "pb_samples_phenotypes_20240425.txt"
)

MMSE_FILE = Path(
    "/project/holstegelab/Share/yaran/sample/100plus_mmse.csv"
)

BED_ROOT = Path(
    "/project/holstegelab/Share/yaran/epigenetics/phased/denovo"
)

OUT_FILE = Path(
    "/project/holstegelab/Share/yixuan/workflow/Step_1/Output/"
    "centenarian_sample_list.tsv"
)


def main():
    print("Step 1: Reading phenotype file...")
    pheno = pd.read_csv(PHENO_FILE, sep="\t", dtype=str)

    required_pheno_cols = [
        "pacbio_id",
        "ID_GWAS",
        "ID_100plus",
        "PHENOTYPE",
        "SEX",
    ]

    missing_cols = [c for c in required_pheno_cols if c not in pheno.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in phenotype file: {missing_cols}")

    print(f"Phenotype rows before filtering: {len(pheno)}")

    print("Step 2: Keeping PHENOTYPE == 'Centenarian'...")
    samples = pheno.loc[
        pheno["PHENOTYPE"].eq("Centenarian")
    ].copy()

    print(f"Rows after Centenarian filter: {len(samples)}")

    print("Step 3: Removing pacbio_id containing 'brain'...")
    brain_mask = samples["pacbio_id"].str.contains(
        "brain",
        case=False,
        na=False,
    )

    if brain_mask.any():
        print("Samples removed because pacbio_id contains 'brain':")
        print(samples.loc[brain_mask, ["pacbio_id", "ID_100plus", "SEX"]].to_string(index=False))

    samples = samples.loc[~brain_mask].copy()

    print(f"Rows after removing brain samples: {len(samples)}")
    print("SEX counts before BED existence check:")
    print(samples["SEX"].value_counts(dropna=False).to_string())

    print("Step 4: Creating bed_path...")
    samples["bed_path"] = samples["pacbio_id"].apply(
        lambda x: str(BED_ROOT / x / f"{x}.cpg.combined.bed")
    )

    print("Step 5: Checking BED file existence...")
    samples["bed_exists"] = samples["bed_path"].apply(lambda x: Path(x).exists())

    missing_bed = samples.loc[~samples["bed_exists"]].copy()

    if len(missing_bed) > 0:
        print(f"BED files not found for {len(missing_bed)} samples:")
        print(
            missing_bed[
                ["pacbio_id", "ID_GWAS", "ID_100plus", "SEX", "bed_path"]
            ].to_string(index=False)
        )
    else:
        print("All BED files exist.")

    samples = samples.loc[samples["bed_exists"]].copy()
    samples = samples.drop(columns=["bed_exists"])

    print(f"Rows after removing samples without BED: {len(samples)}")
    print("SEX counts after BED existence check:")
    print(samples["SEX"].value_counts(dropna=False).to_string())

    print("Step 6: Reading MMSE file...")
    mmse = pd.read_csv(MMSE_FILE, dtype=str)

    required_mmse_cols = [
        "studysubjectid",
        "harmonized_contact_nr",
        "mmse_score_imputed",
    ]

    missing_mmse_cols = [c for c in required_mmse_cols if c not in mmse.columns]
    if missing_mmse_cols:
        raise ValueError(f"Missing columns in MMSE file: {missing_mmse_cols}")

    print(f"MMSE rows before contact filtering: {len(mmse)}")

    print("Step 7: Keeping harmonized_contact_nr == 0...")
    mmse0 = mmse.loc[
        mmse["harmonized_contact_nr"].astype(str).eq("0"),
        ["studysubjectid", "mmse_score_imputed"],
    ].copy()

    print(f"MMSE rows after contact filtering: {len(mmse0)}")

    print("Step 8: Left joining MMSE by ID_100plus <-> studysubjectid...")
    out = samples.merge(
        mmse0,
        how="left",
        left_on="ID_100plus",
        right_on="studysubjectid",
    )

    out = out[
        [
            "pacbio_id",
            "ID_GWAS",
            "ID_100plus",
            "SEX",
            "mmse_score_imputed",
            "bed_path",
        ]
    ].copy()

    print("Step 9: Writing output sample list...")
    out.to_csv(OUT_FILE, sep="\t", index=False)

    print(f"Output written to: {OUT_FILE}")
    print(f"Final sample count: {len(out)}")

    print("Final SEX value_counts:")
    print(out["SEX"].value_counts(dropna=False).to_string())

    print("MMSE missing count:")
    print(out["mmse_score_imputed"].isna().sum())


if __name__ == "__main__":
    main()
