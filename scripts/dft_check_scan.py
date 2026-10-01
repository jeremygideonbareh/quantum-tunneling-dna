"""DFT single points on an xTB scan (PySCF, free alternative to ORCA).

Re-evaluates every image of an xTB scan (scan.xyz) with DFT and writes the
DFT profile next to the xTB one, so the barrier / reaction energy / ranking
can be compared directly.

    python scripts/dft_check_scan.py results/smoke_gc --xc b3lyp --basis def2-svp
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
from ase.io import read
from pyscf import dft, gto, lib

HARTREE_TO_KCAL = 627.509474


def single_point(atoms, xc, basis, disp):
    mol = gto.M(atom=[(s, p) for s, p in zip(atoms.get_chemical_symbols(), atoms.positions)],
                basis=basis, charge=0, spin=0, verbose=0)
    mf = dft.RKS(mol, xc=xc).density_fit()
    if disp:
        mf.disp = disp
    mf.conv_tol = 1e-9
    mf.grids.level = 3
    e = mf.kernel()
    if not mf.converged:
        raise RuntimeError("SCF did not converge")
    return e


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scan_dir")
    ap.add_argument("--xc", default="b3lyp")
    ap.add_argument("--basis", default="def2-svp")
    ap.add_argument("--disp", default="d3bj", help="'' to disable")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    lib.num_threads(args.threads)

    d = Path(args.scan_dir)
    images = read(d / "scan.xyz", index=":")
    xtb = json.loads((d / "summary.json").read_text())

    disp = args.disp or None  # D3 needs `pip install pyscf-dispersion`

    energies = []
    for i, atoms in enumerate(images):
        t0 = time.time()
        energies.append(single_point(atoms, args.xc, args.basis, disp))
        rel = (energies[-1] - energies[0]) * HARTREE_TO_KCAL
        print(f"image {i:2d}  E_rel = {rel:7.2f} kcal/mol  ({time.time() - t0:.0f} s)", flush=True)

    e = (np.array(energies) - energies[0]) * HARTREE_TO_KCAL
    e_xtb = np.array([p["E_kcal"] for p in xtb["scan"]])
    i_ts = int(np.argmax(e))
    out = {
        "level": f"{args.xc}{'-' + disp if disp else ''}/{args.basis} // GFN2-xTB geometries",
        "barrier_kcal": float(e[i_ts]),
        "reaction_energy_kcal": float(e[-1]),
        "reverse_barrier_kcal": float(e[i_ts] - e[-1]),
        "tautomer_is_minimum_on_xtb_path": bool(i_ts < len(e) - 1),
        "xtb_barrier_kcal": xtb["barrier_kcal"],
        "xtb_reaction_energy_kcal": xtb["reaction_energy_kcal"],
        "profile": [{"X_A": p["X_A"], "E_xtb": float(a), "E_dft": float(b)}
                    for p, a, b in zip(xtb["scan"], e_xtb, e)],
    }
    tag = f"{args.xc}_{args.basis}".replace("-", "").lower()
    (d / f"dft_{tag}.json").write_text(json.dumps(out, indent=2))
    print(json.dumps({k: v for k, v in out.items() if k != "profile"}, indent=2))


if __name__ == "__main__":
    main()
