"""Phase 1 smoke test: G:C double proton transfer scan with GFN2-xTB.

Single methyl-capped G:C pair, gas phase.

1. Optimise canonical G:C and the G*:C* tautomer (both protons moved).
2. Relaxed scan along X = xi1 + xi2, where xi = d(donor-H) - d(H...acceptor)
   for the two transferring protons (G N1-H -> C N3, C N4-H -> G O6).
   Constraining only the sum lets the path be stepwise or concerted.
3. Report forward/reverse barrier, reaction energy and a WKB tunnelling
   factor at 310 K along the mass-weighted path.

This is a pipeline test, not a result: one pair, no stacking, no backbone,
no solvent, and xTB barriers for proton transfer must be checked against DFT.

    python scripts/smoke_test_gc_scan.py --out results/smoke_gc
"""

import argparse
import json
from pathlib import Path

import numpy as np
from ase.constraints import FixInternals
from ase.io import write
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from tblite.ase import TBLite

from qtdna.basepair import build_gc_pair
from qtdna.tunneling import mass_weighted_path, wkb_kappa

EV_TO_KCAL = 23.0605


def xi(atoms, d, h, a):
    return atoms.get_distance(d, h) - atoms.get_distance(h, a)


def place_proton(atoms, d, h, a, target):
    """Put h on the d->a line so that xi = target."""
    pd, pa = atoms.positions[d], atoms.positions[a]
    L = np.linalg.norm(pa - pd)
    atoms.positions[h] = pd + (pa - pd) / L * (L + target) / 2.0


def zpe(atoms, name, workdir):
    """Harmonic zero-point energy (eV) and number of imaginary modes."""
    vib = Vibrations(atoms, name=str(workdir / f"vib_{name}"))
    vib.clean()
    vib.run()
    e = vib.get_energies()
    vib.clean()
    n_imag = int(np.sum(np.abs(e.imag) > 1e-3))
    return float(0.5 * np.sum(e.real[6:])), n_imag


def relax(atoms, fmax, steps=400):
    BFGS(atoms, logfile=None).run(fmax=fmax, steps=steps)
    return atoms.get_potential_energy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/smoke_gc")
    ap.add_argument("--method", default="GFN2-xTB")
    ap.add_argument("--points", type=int, default=21)
    ap.add_argument("--fmax", type=float, default=0.01)
    ap.add_argument("--no-zpe", action="store_true")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    atoms, idx = build_gc_pair()
    atoms.calc = TBLite(method=args.method, verbosity=0)
    e_free = relax(atoms, args.fmax)
    write(out / "gc_opt.xyz", atoms)

    hb = {k: float(atoms.get_distance(idx[a], idx[b])) for k, (a, b) in {
        "O6...N4": ("G_O6", "C_N4"), "N1...N3": ("G_N1", "C_N3"),
        "N2...O2": ("G_N2", "C_O2")}.items()}
    print("H-bond heavy-atom distances (A):", {k: round(v, 3) for k, v in hb.items()})

    triples = [(idx["G_N1"], idx["G_H1"], idx["C_N3"]),
               (idx["C_N4"], idx["C_H4"], idx["G_O6"])]
    x_r = [xi(atoms, *t) for t in triples]

    taut = atoms.copy()
    taut.calc = TBLite(method=args.method, verbosity=0)
    for t in triples:
        place_proton(taut, *t, -xi(atoms, *t))
    e_taut = relax(taut, args.fmax)
    x_p = [xi(taut, *t) for t in triples]
    taut_is_minimum = sum(x_p) > 0
    write(out / "gcstar_opt.xyz", taut)

    reactant = atoms.positions.copy()
    combo = [[d, h, 1.0] for d, h, _ in triples] + [[h, a, -1.0] for _, h, a in triples]
    fracs = np.linspace(0.0, 1.0, args.points)
    energies, images = [e_free], [reactant]
    for frac in fracs[1:-1]:
        targets = [(1 - frac) * a + frac * b for a, b in zip(x_r, x_p)]
        atoms.set_constraint()
        for t, x in zip(triples, targets):
            place_proton(atoms, *t, x)
        atoms.set_constraint(FixInternals(bondcombos=[[sum(targets), combo]]))
        energies.append(relax(atoms, args.fmax))
        images.append(atoms.positions.copy())
    atoms.set_constraint()
    energies.append(e_taut)
    images.append(taut.positions.copy())
    grid = (1 - fracs) * sum(x_r) + fracs * sum(x_p)
    write(out / "scan.xyz", [atoms.__class__(atoms.symbols, positions=p) for p in images])

    if not args.no_zpe:
        z_r, imag_r = zpe(atoms, "gc", out)
        z_p, imag_p = zpe(taut, "gcstar", out)

    e = np.array(energies) - e_free
    i_ts = int(np.argmax(e))
    product_is_minimum = taut_is_minimum and i_ts < len(e) - 1
    # Two limits for the tunnelling path length:
    #   adiabatic - every atom follows the relaxed scan (heavy atoms lengthen s)
    #   sudden    - only the two protons move, frame frozen (Slocombe-style regime)
    s = mass_weighted_path(images, atoms.get_masses())
    proton_mass = np.zeros(len(atoms))
    for _, h, _ in triples:
        proton_mass[h] = atoms.get_masses()[h]
    s_sudden = mass_weighted_path(images, proton_mass)
    summary = {
        "method": args.method,
        "model": "methyl-capped G:C, gas phase, relaxed scan of xi(N1-H) + xi(N4-H)",
        "hbond_distances_A": hb,
        "barrier_kcal": float(e[i_ts] * EV_TO_KCAL),
        "reaction_energy_kcal": float(e[-1] * EV_TO_KCAL),
        "reverse_barrier_kcal": float((e[i_ts] - e[-1]) * EV_TO_KCAL),
        "tautomer_is_minimum": bool(product_is_minimum),
        "wkb_kappa_310K_adiabatic": wkb_kappa(s, e, T=310.0) if product_is_minimum else None,
        "wkb_kappa_310K_sudden": wkb_kappa(s_sudden, e, T=310.0) if product_is_minimum else None,
        "zpe_corrected_reaction_energy_kcal": None if args.no_zpe else
            float((e_taut + z_p - e_free - z_r) * EV_TO_KCAL),
        "imaginary_modes_gc_gcstar": None if args.no_zpe else [imag_r, imag_p],
        "scan": [{"X_A": float(d), "E_kcal": float(x * EV_TO_KCAL), "s": float(si)}
                 for d, x, si in zip(grid, e, s)],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    for p in summary["scan"]:
        print(f"  X = {p['X_A']:+.3f} A   E = {p['E_kcal']:7.2f} kcal/mol")
    print(json.dumps({k: v for k, v in summary.items() if k != "scan"}, indent=2))


if __name__ == "__main__":
    main()
