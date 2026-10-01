"""Per-site transition rates by sequence context for E. coli MA data.

Input: the Jago et al. PNAS 2026 compilation (Lagator group), which collects
120,208 base substitutions from Lee 2012, Foster 2015, Long 2016,
Tincher 2017, Foster 2018 and Niccum 2018, plus the MG1655 reference.
Get it with: git clone https://github.com/Lagator-Group/extended-sequence-context data/raw/extended-sequence-context

Output (data/processed/): ecoli_rates_k{3,5}.csv with one row per
(repair group, context): n_mut, n_sites, rate, rate_norm (rate / mean rate of
that centre within the group), cpg flag.

    python scripts/ecoli_context_rates.py
"""

import argparse
from pathlib import Path

import pandas as pd

from qtdna.contexts import is_cpg
from qtdna.mutations import count_transitions, read_fasta

GROUPS = {
    "proofreading(+) MMR(-)": "MMR_def",          # primary: replication errors MMR would fix
    "proofreading(-) MMR(-)": "MMR_def_pol_exo",  # mutD5 + mutL: raw polymerase errors
    "proofreading(-) MMR(+)": "pol_exo",          # mutD5 alone
    "proofreading(+) MMR(+)": "WT_like",          # too few mutations for context analysis
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="data/raw/extended-sequence-context")
    ap.add_argument("--out", default="data/processed")
    args = ap.parse_args()
    src, out = Path(args.src), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    muts = pd.read_csv(src / "supplementary_data_1.csv", skiprows=1)  # row 0 = descriptions
    genome = read_fasta(src / "mg1655_NC_000913.3.fasta")
    (chrom,) = genome
    muts["chrom"] = chrom

    for k in (3, 5):
        tables = []
        for group, label in GROUPS.items():
            df = count_transitions(muts[muts["group"] == group], genome, k)
            df.insert(0, "group", label)
            df["centre"] = df["context"].str[k // 2]
            df["cpg"] = df["context"].map(is_cpg)
            df["rate_norm"] = df["rate"] / df.groupby("centre")["rate"].transform("mean")
            tables.append(df)
        res = pd.concat(tables, ignore_index=True)
        res.to_csv(out / f"ecoli_rates_k{k}.csv", index=False)
        print(f"k={k}: wrote {len(res)} rows; transitions per group:",
              res.groupby("group")["n_mut"].sum().to_dict())


if __name__ == "__main__":
    main()
