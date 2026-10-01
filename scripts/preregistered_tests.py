"""Pre-registered fingerprint-vs-mutation tests (docs/preregistration.md).

This script refuses to run without the OSF registration URL: the whole point
of pre-registration is that nobody looks at these correlations first.

    python scripts/preregistered_tests.py --osf https://osf.io/xxxxx --env gas

Implements H1/H2 (E. coli), H3 (COSMIC part; yeast/human added when the data
arrive), H4 (proofreading groups). Writes results/preregistered_<env>.csv.
"""

import argparse
import itertools
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr

KT = 0.0019872 * 310.0  # kcal/mol at 310 K
MMR_SIGS = ["SBS6", "SBS15", "SBS21", "SBS26", "SBS44"]
PROOF_SIGS = ["SBS10a", "SBS10b", "SBS10c", "SBS10d"]


def primary_feature(fp: pd.DataFrame) -> pd.Series:
    """F = ln(kappa_sudden) - dE_barrier / RT (pre-registered primary).

    Where there is no tautomer well deeper than RT, kappa is undefined; use
    kappa = 1 (pre-registered rule) and keep the context.
    """
    kappa = fp["kappa_sudden"].astype(float).fillna(1.0)
    kappa[fp["reverse_barrier_kcal"] < KT] = 1.0
    return np.log(kappa) - fp["barrier_kcal"] / KT


def exact_spearman(x, y, n_perm=200_000, seed=0):
    """One-sided (positive) Spearman with permutation p (exact if n <= 9)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    rho = spearmanr(x, y).statistic
    rx = pd.Series(x).rank().values
    ry = pd.Series(y).rank().values
    if len(x) <= 9:
        perms = (np.array(p) for p in itertools.permutations(ry))
        null = np.array([np.corrcoef(rx, p)[0, 1] for p in perms])
    else:
        rng = np.random.default_rng(seed)
        null = np.array([np.corrcoef(rx, rng.permutation(ry))[0, 1] for _ in range(n_perm)])
    p = (np.sum(null >= rho) + 1) / (len(null) + 1)
    return rho, p


def holm(pvals):
    order = np.argsort(pvals)
    adj = np.empty(len(pvals))
    running = 0.0
    for i, idx in enumerate(order):
        running = max(running, min(1.0, (len(pvals) - i) * pvals[idx]))
        adj[idx] = running
    return adj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--osf", required=True, help="OSF pre-registration URL (must exist first)")
    ap.add_argument("--env", default="gas")
    ap.add_argument("--h2", choices=["A", "B"], required=True,
                    help="A: T-centred predictor = -E_half (single-PT energy); B: F as for G:C")
    args = ap.parse_args()
    if not args.osf.startswith("https://osf.io/"):
        raise SystemExit("Post docs/preregistration.md on OSF first and pass its URL.")

    fp = pd.read_csv(f"results/fingerprint_k3_{args.env}.csv").set_index("context")
    fp["F"] = primary_feature(fp)
    if args.h2 == "A":  # lower half-transfer energy -> more errors, so use -E_half
        t = fp["centre"] == "T"
        fp.loc[t, "F"] = -fp.loc[t, "E_half_kcal"]
    eco = pd.read_csv("data/processed/ecoli_rates_k3.csv")
    cos = pd.read_csv("data/processed/cosmic_sitewise_k3.csv")

    rows = []
    for centre, hyp in (("C", "H1"), ("T", "H2")):
        ctx = fp.index[fp["centre"] == centre]
        for group, h in (("MMR_def", hyp), ("MMR_def_pol_exo", "H4"), ("pol_exo", "H4")):
            r = eco[eco.group == group].set_index("context").loc[ctx, "rate"]
            rho, p = exact_spearman(fp.loc[ctx, "F"], r)
            rows.append(dict(hyp=h, dataset=f"E. coli {group}", centre=centre, n=len(ctx), rho=rho, p=p))
        for sig in MMR_SIGS + PROOF_SIGS:
            s = cos[(cos.signature == sig) & (cos.centre == centre) & ~cos.cpg].set_index("context")
            c = ctx.intersection(s.index)
            rho, p = exact_spearman(fp.loc[c, "F"], s.loc[c, "rate"])
            rows.append(dict(hyp="H3" if sig in MMR_SIGS else "H4", dataset=f"COSMIC {sig}",
                             centre=centre, n=len(c), rho=rho, p=p))
    res = pd.DataFrame(rows)
    main_tests = res.hyp.isin(["H1", "H2", "H3"])
    res.loc[main_tests, "p_holm"] = holm(res.loc[main_tests, "p"].values)
    for centre in "CT":
        sub = res[(res.hyp == "H3") & (res.centre == centre)]
        z = norm.isf(sub["p"]).sum() / math.sqrt(len(sub))
        rows_meta = dict(hyp="H3-meta", dataset="Stouffer COSMIC MMR", centre=centre,
                         n=int(sub.n.iloc[0]), rho=np.nan, p=float(norm.sf(z)))
        res = pd.concat([res, pd.DataFrame([rows_meta])], ignore_index=True)
    out = Path("results") / f"preregistered_{args.env}.csv"
    res.to_csv(out, index=False)
    print(res.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print(f"\nOSF registration: {args.osf} (H2 option {args.h2})\nwrote {out}")


if __name__ == "__main__":
    main()
