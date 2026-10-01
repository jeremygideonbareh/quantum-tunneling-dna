"""Baseline QM features per context from DNAkmerQM (Masuda & Sahakyan 2024).

DNAkmerQM gives B/A/Z-DNA features for every 7-mer. We keep B-DNA energy and
energy-difference features, take the pyrimidine-centred 7-mers (each duplex
once), average over flanks to get 3-mer and 5-mer values, standardise, and
add the first three principal components. These are the "generic QM"
baseline for the nested-model test (pre-registration H5).

    python scripts/dnakmerqm_baseline.py
"""

import argparse
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


def load(zip_path):
    with zipfile.ZipFile(zip_path) as z:
        (name,) = [n for n in z.namelist() if n.endswith(".txt")]
        with z.open(name) as fh:
            return pd.read_csv(fh).set_index("seq")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="data/raw/DNAkmerQM/7-mer")
    ap.add_argument("--out", default="data/processed")
    args = ap.parse_args()
    src = Path(args.src)

    feats = load(src / "energy.zip").join(load(src / "denergy.zip"), rsuffix="_d")
    feats = feats[[c for c in feats.columns if c.startswith("B_")]]
    feats = feats.loc[[s for s in feats.index if s[3] in "CT"]]
    feats = feats.loc[:, feats.std() > 0].dropna(axis=1)

    for k in (3, 5):
        flank = k // 2
        key = feats.index.str[3 - flank:4 + flank]
        agg = feats.groupby(key).mean()
        z = (agg - agg.mean()) / agg.std()
        u, s, vt = np.linalg.svd(z.values, full_matrices=False)
        pcs = pd.DataFrame(u[:, :3] * s[:3], index=agg.index, columns=["PC1", "PC2", "PC3"])
        explained = (s ** 2 / (s ** 2).sum())[:3]
        out = pcs.join(agg)
        out.index.name = "context"
        out.to_csv(Path(args.out) / f"dnakmerqm_B_k{k}.csv")
        print(f"k={k}: {len(out)} contexts, {agg.shape[1]} features, "
              f"PC1-3 explain {np.round(explained * 100, 1).tolist()} %")


if __name__ == "__main__":
    main()
