"""Turn COSMIC SBS signatures into per-site transition rates by trinucleotide.

A COSMIC signature gives P(mutation type | signature) over the 96 classes,
so it mixes context mutability with how common each trinucleotide is in the
genome. Dividing by GRCh38 trinucleotide counts gives relative per-site
mutability, which is what a physical fingerprint should predict.

Inputs (data/raw/cosmic/), both shipped inside public PyPI packages:
  COSMIC_v3.4_SBS_GRCh38.txt   (SigProfilerAssignment)
  context_counts_GRCh38_96.csv (SigProfilerMatrixGenerator)

    python scripts/cosmic_sitewise_rates.py
"""

import argparse
from pathlib import Path

import pandas as pd

from qtdna.contexts import TRANSITION, is_cpg

SIGNATURES = {
    # mismatch-repair deficiency
    "SBS6": "MMR", "SBS14": "POLE+MMR", "SBS15": "MMR", "SBS20": "POLD1+MMR",
    "SBS21": "MMR", "SBS26": "MMR", "SBS44": "MMR",
    # polymerase proofreading
    "SBS10a": "POLE", "SBS10b": "POLE", "SBS10c": "POLD1", "SBS10d": "POLD1",
    # controls
    "SBS1": "CpG deamination (control)", "SBS5": "clock-like (control)",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="data/raw/cosmic")
    ap.add_argument("--out", default="data/processed/cosmic_sitewise_k3.csv")
    args = ap.parse_args()
    src = Path(args.src)

    sig = pd.read_csv(src / "COSMIC_v3.4_SBS_GRCh38.txt", sep="\t", index_col=0)
    counts = pd.read_csv(src / "context_counts_GRCh38_96.csv", index_col=0).sum(axis=1)

    rows = []
    for label, p in sig[list(SIGNATURES)].iterrows():
        ctx = label[0] + label[2] + label[6]          # A[C>T]G -> ACG
        ref, alt = label[2], label[4]
        if TRANSITION[ref] != alt:
            continue
        for name, val in p.items():
            rows.append({"signature": name, "process": SIGNATURES[name], "context": ctx,
                         "centre": ref, "prob": val, "n_sites": int(counts[ctx]),
                         "cpg": is_cpg(ctx)})
    df = pd.DataFrame(rows)
    df["rate"] = df["prob"] / df["n_sites"]
    df["rate_norm"] = df["rate"] / df.groupby(["signature", "centre"])["rate"].transform("mean")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"wrote {len(df)} rows for {df.signature.nunique()} signatures -> {args.out}")


if __name__ == "__main__":
    main()
