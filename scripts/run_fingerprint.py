"""Compute the tunnelling fingerprint for every context (resumable).

    python scripts/run_fingerprint.py --k 3 --env gas
    python scripts/run_fingerprint.py --k 3 --env water      # ALPB implicit water
    python scripts/run_fingerprint.py --k 3 --env gas ACG GTA  # just these

Reads structures/k{k}/<ctx>.pdb (scripts/build_bdna_clusters.py), writes
results/fingerprint_k{k}_{env}/<ctx>.json and <ctx>_scan.xyz, then a combined
fingerprint_k{k}_{env}.csv.
"""

import argparse
import json
import time
from pathlib import Path

import pandas as pd
from ase.io import write
from tblite.ase import TBLite

from qtdna.contexts import contexts
from qtdna.scan import dpt_scan, frozen_indices, proton_triples, read_pdb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--env", choices=["gas", "water"], default="gas")
    ap.add_argument("--points", type=int, default=15)
    ap.add_argument("--fmax", type=float, default=0.02)
    ap.add_argument("--structures", default="structures")
    ap.add_argument("ctx", nargs="*")
    args = ap.parse_args()

    out = Path("results") / f"fingerprint_k{args.k}_{args.env}"
    out.mkdir(parents=True, exist_ok=True)
    solv = {"solvation": ("alpb", "water")} if args.env == "water" else {}

    def make_calc():
        return TBLite(method="GFN2-xTB", verbosity=0, **solv)

    for ctx in args.ctx or contexts(args.k):
        f = out / f"{ctx}.json"
        if f.exists():
            continue
        t0 = time.time()
        cl = read_pdb(Path(args.structures) / f"k{args.k}" / f"{ctx}.pdb")
        centre = ctx[args.k // 2]
        triples = proton_triples(cl, args.k, centre)
        res, images = dpt_scan(cl.atoms, triples, frozen_indices(cl, args.k), make_calc,
                               points=args.points, fmax=args.fmax)
        res.update(context=ctx, centre=centre, env=args.env, method="GFN2-xTB",
                   seconds=round(time.time() - t0))
        f.write_text(json.dumps(res, indent=2))
        frames = [cl.atoms.copy() for _ in images]
        for fr, p in zip(frames, images):
            fr.positions = p
        write(out / f"{ctx}_scan.xyz", frames)
        print(f"{ctx}: dE={res['reaction_energy_kcal']:.2f} barrier={res['barrier_kcal']:.2f} "
              f"reverse={res['reverse_barrier_kcal']:.2f} ({res['seconds']} s)", flush=True)

    rows = [json.loads(p.read_text()) for p in sorted(out.glob("*.json"))]
    keep = ["context", "centre", "env", "reaction_energy_kcal", "barrier_kcal",
            "reverse_barrier_kcal", "tautomer_is_minimum", "kappa_adiabatic", "kappa_sudden"]
    pd.DataFrame(rows)[keep].to_csv(out.parent / f"fingerprint_k{args.k}_{args.env}.csv", index=False)


if __name__ == "__main__":
    main()
