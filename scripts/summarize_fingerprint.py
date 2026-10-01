"""Descriptive summary of a fingerprint run (no mutation data involved).

Rebuilds results/fingerprint_k{k}_{env}.csv from the per-context JSONs and
writes results/fingerprint_k{k}_{env}.md with ranges per centre. Deliberately
does NOT touch mutation data or the pre-registered hypotheses (H6 included).

    python scripts/summarize_fingerprint.py --k 3 --env gas
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from qtdna.contexts import contexts

KT = 0.0019872 * 310.0
COLS = ["reaction_energy_kcal", "barrier_kcal", "reverse_barrier_kcal",
        "E_half_kcal", "kappa_adiabatic", "kappa_sudden", "F"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--env", default="gas")
    args = ap.parse_args()
    d = Path("results") / f"fingerprint_k{args.k}_{args.env}"
    rows = [json.loads(p.read_text()) for p in sorted(d.glob("*.json"))]
    df = pd.DataFrame(rows).set_index("context")
    missing = sorted(set(contexts(args.k)) - set(df.index))
    df["well_deeper_than_RT"] = df["tautomer_is_minimum"] & (df["reverse_barrier_kcal"] > KT)
    kappa = df["kappa_sudden"].astype(float).fillna(1.0).where(df["well_deeper_than_RT"], 1.0)
    df["F"] = np.log(kappa) - df["barrier_kcal"] / KT
    # energy at the half-transfer point X = 0 (single-proton-transfer region)
    df["E_half_kcal"] = [float(np.interp(0.0, r["X"], r["profile_kcal"])) for r in rows]
    keep = ["centre", "env", "tautomer_is_minimum", "well_deeper_than_RT"] + COLS
    df[keep].sort_index().to_csv(d.parent / f"fingerprint_k{args.k}_{args.env}.csv")

    lines = [f"# Fingerprint summary: k={args.k}, {args.env}, GFN2-xTB",
             "", f"Contexts computed: {len(df)} / {len(contexts(args.k))}"
             + (f" (missing: {', '.join(missing)})" if missing else ""),
             f"Tautomer is a minimum in {int(df.tautomer_is_minimum.sum())} / {len(df)} contexts; "
             f"well deeper than RT ({KT:.2f} kcal/mol) in {int(df.well_deeper_than_RT.sum())}.", ""]
    for centre, name in (("C", "G:C → G*:C* (C-centred)"), ("T", "A:T → A*:T* (T-centred)")):
        g = df[df.centre == centre]
        if g.empty:
            continue
        lines += [f"## {name}, n = {len(g)}", "", "| quantity | min | median | max | range |",
                  "|---|---|---|---|---|"]
        for c in COLS:
            v = g[c].astype(float).dropna()
            if len(v):
                lines.append(f"| {c} | {v.min():.2f} | {v.median():.2f} | {v.max():.2f} | {v.max() - v.min():.2f} |")
        lines.append(f"\nTautomer well deeper than RT: {int(g.well_deeper_than_RT.sum())} / {len(g)}")
        lo, hi = g["F"].idxmin(), g["F"].idxmax()
        lines += ["", f"Lowest / highest F: {lo} ({g.F[lo]:.2f}) / {hi} ({g.F[hi]:.2f}); "
                  f"spread corresponds to a {np.exp(g.F[hi] - g.F[lo]):.1f}-fold rate difference.", ""]
    (d.parent / f"fingerprint_k{args.k}_{args.env}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
