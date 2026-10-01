"""Double-proton-transfer scans on a central base pair inside a B-DNA cluster.

Reaction coordinate: X = xi1 + xi2 with xi = d(D-H) - d(H...A) for the two
transferring protons. Constraining only the sum lets the path be stepwise or
concerted. Flanking pairs and every C1' methyl cap are frozen at the ideal
B-DNA positions so the stack cannot fall apart.

    G*:C* : G N1-H -> C N3   and  C N4-H -> G O6
    A*:T* : T N3-H -> A N1   and  A N6-H -> T O4
"""

from dataclasses import dataclass

import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms, FixInternals
from ase.optimize import BFGS

from .tunneling import mass_weighted_path, wkb_kappa

EV_TO_KCAL = 23.0605

# (donor residue, donor atom, acceptor residue, acceptor atom); residues "pyr"/"pur"
PROTONS = {
    "C": [("pur", "N1", "pyr", "N3"), ("pyr", "N4", "pur", "O6")],
    "T": [("pyr", "N3", "pur", "N1"), ("pur", "N6", "pyr", "O4")],
}


@dataclass
class Cluster:
    atoms: Atoms
    names: list          # PDB atom names
    chain: list
    resi: list


def read_pdb(path) -> Cluster:
    sym, pos, names, chain, resi = [], [], [], [], []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            names.append(line[12:16].strip())
            chain.append(line[21])
            resi.append(int(line[22:26]))
            pos.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            sym.append(line[76:78].strip() or names[-1][0])
    return Cluster(Atoms(sym, positions=pos), names, chain, resi)


def central_pair(cl: Cluster, k: int):
    """Indices of the central pyrimidine (strand A) and purine (strand B) residues."""
    i = k // 2 + 1
    pyr = [n for n in range(len(cl.names)) if cl.chain[n] == "A" and cl.resi[n] == i]
    pur = [n for n in range(len(cl.names)) if cl.chain[n] == "B" and cl.resi[n] == -i]
    return pyr, pur


def _find(cl, idx, name):
    (j,) = [n for n in idx if cl.names[n] == name]
    return j


def proton_triples(cl: Cluster, k: int, centre: str):
    pyr, pur = central_pair(cl, k)
    res = {"pyr": pyr, "pur": pur}
    pos = cl.atoms.positions
    triples = []
    for dres, dname, ares, aname in PROTONS[centre]:
        d, a = _find(cl, res[dres], dname), _find(cl, res[ares], aname)
        hs = [n for n in res[dres] if cl.atoms[n].symbol == "H"
              and np.linalg.norm(pos[n] - pos[d]) < 1.2]
        h = min(hs, key=lambda n: np.linalg.norm(pos[n] - pos[a]))
        triples.append((d, h, a))
    return triples


def frozen_indices(cl: Cluster, k: int):
    pyr, pur = central_pair(cl, k)
    central = set(pyr) | set(pur)
    return [n for n in range(len(cl.names))
            if (n not in central and cl.atoms[n].symbol != "H") or cl.names[n] == "C1'"]


def xi(atoms, d, h, a):
    return atoms.get_distance(d, h) - atoms.get_distance(h, a)


def place_proton(atoms, d, h, a, target):
    pd, pa = atoms.positions[d], atoms.positions[a]
    L = np.linalg.norm(pa - pd)
    atoms.positions[h] = pd + (pa - pd) / L * (L + target) / 2.0


def relax(atoms, fmax, steps=600):
    BFGS(atoms, logfile=None).run(fmax=fmax, steps=steps)
    return atoms.get_potential_energy()


def dpt_scan(atoms, triples, frozen, make_calc, points=15, fmax=0.02):
    """Return a dict of fingerprint features plus the raw profile."""
    atoms = atoms.copy()
    atoms.calc = make_calc()
    atoms.set_constraint(FixAtoms(indices=frozen))
    e_r = relax(atoms, fmax)
    x_r = [xi(atoms, *t) for t in triples]
    reactant = atoms.positions.copy()

    taut = atoms.copy()
    taut.calc = make_calc()
    taut.set_constraint(FixAtoms(indices=frozen))
    for t in triples:
        place_proton(taut, *t, -xi(atoms, *t))
    e_p = relax(taut, fmax)
    x_p = [xi(taut, *t) for t in triples]

    combo = [[d, h, 1.0] for d, h, _ in triples] + [[h, a, -1.0] for _, h, a in triples]
    fracs = np.linspace(0.0, 1.0, points)
    energies, images = [e_r], [reactant]
    for f in fracs[1:-1]:
        targets = [(1 - f) * a + f * b for a, b in zip(x_r, x_p)]
        atoms.set_constraint()
        for t, x in zip(triples, targets):
            place_proton(atoms, *t, x)
        atoms.set_constraint([FixAtoms(indices=frozen),
                              FixInternals(bondcombos=[[sum(targets), combo]])])
        energies.append(relax(atoms, fmax))
        images.append(atoms.positions.copy())
    energies.append(e_p)
    images.append(taut.positions.copy())

    e = (np.array(energies) - e_r)
    i_ts = int(np.argmax(e))
    tautomer_is_minimum = sum(x_p) > 0 and 0 < i_ts < len(e) - 1
    masses = atoms.get_masses()
    s_ad = mass_weighted_path(images, masses)
    m_h = np.zeros(len(atoms))
    for _, h, _ in triples:
        m_h[h] = masses[h]
    s_sd = mass_weighted_path(images, m_h)
    return {
        "barrier_kcal": float(e[i_ts] * EV_TO_KCAL),
        "reaction_energy_kcal": float(e[-1] * EV_TO_KCAL),
        "reverse_barrier_kcal": float((e[i_ts] - e[-1]) * EV_TO_KCAL),
        "tautomer_is_minimum": bool(tautomer_is_minimum),
        "kappa_adiabatic": wkb_kappa(s_ad, e) if tautomer_is_minimum else None,
        "kappa_sudden": wkb_kappa(s_sd, e) if tautomer_is_minimum else None,
        "xi_reactant": [float(x) for x in x_r],
        "xi_tautomer": [float(x) for x in x_p],
        "profile_kcal": [float(v * EV_TO_KCAL) for v in e],
        "X": [float((1 - f) * sum(x_r) + f * sum(x_p)) for f in fracs],
    }, images
