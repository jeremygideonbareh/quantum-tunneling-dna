"""Build a methyl-capped Watson-Crick base pair without external DNA tools.

The N9/N1 methyl caps stand in for the sugars. Geometry is a rough guess
(RDKit embed + three parallel H-bonds) that xTB then relaxes. This is only for
the single-pair smoke test; stacked 3-bp / 5-bp models come from 3DNA `fiber`
or DNAkmerQM geometries (see docs/research_plan.md, Phase 1).
"""

import numpy as np
from ase import Atoms
from rdkit import Chem
from rdkit.Chem import AllChem

# Each entry: SMILES with atom-map numbers on the atoms we need to address.
# Map numbers label the heavy atoms of the three H-bonds.
G_SMILES = "Cn1cnc2c1nc([NH2:4])[nH:2]c2=[O:1]"   # 9-methylguanine
C_SMILES = "Cn1ccc([NH2:5])[n:6]c1=[O:7]"         # 1-methylcytosine
# Pairing (G atom map -> C atom map): O6...H-N4, N1-H...N3, N2-H...O2
PAIRS = [(1, 5), (2, 6), (4, 7)]


def _embed(smiles: str, seed: int = 7):
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(mol)
    pos = mol.GetConformer().GetPositions()
    sym = [a.GetSymbol() for a in mol.GetAtoms()]
    amap = {a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum()}
    hyd = {m: [n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetSymbol() == "H"]
           for m, i in amap.items()}
    ring = pos[list(mol.GetRingInfo().AtomRings()[0])].mean(0)
    return sym, pos, amap, hyd, ring


def _closest(pos, candidates, ref):
    return min(candidates, key=lambda i: np.linalg.norm(pos[i] - ref))


def _kabsch(p, q):
    """Rotation R, translation t minimising |R p + t - q|."""
    pc, qc = p.mean(0), q.mean(0)
    h = (p - pc).T @ (q - qc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1, 1, d]) @ u.T
    return r, qc - r @ pc


def build_gc_pair():
    """Return (atoms, idx) where idx names the atoms used in the DPT scan."""
    gs, gp, gm, gh, _ = _embed(G_SMILES)
    cs, cp, cm, ch, cring = _embed(C_SMILES)

    h1 = gh[2][0]
    n1 = gp[gm[2]]
    d = (gp[h1] - n1) / np.linalg.norm(gp[h1] - n1)
    h2 = _closest(gp, gh[4], gp[h1])  # N2 hydrogen on the Watson-Crick side
    targets = np.array([gp[gm[1]] + 2.85 * d,     # C N4 opposite G O6
                        gp[gm[2]] + 2.95 * d,     # C N3 opposite G N1
                        gp[gm[4]] + 2.90 * d])    # C O2 opposite G N2
    src = np.array([cp[cm[5]], cp[cm[6]], cp[cm[7]]])
    # The three edge atoms are nearly collinear, so also pin the cytosine ring
    # centroid: in-plane with guanine, on the far side of the H-bond edge.
    depth = np.linalg.norm(cring - src.mean(0))
    targets = np.vstack([targets, targets.mean(0) + depth * d])
    src = np.vstack([src, cring])
    r, t = _kabsch(src, targets)
    cp = cp @ r.T + t

    # N4 hydrogen that points at G O6
    h4 = _closest(cp, ch[5], gp[gm[1]])

    ng = len(gs)
    atoms = Atoms(gs + cs, positions=np.vstack([gp, cp]))
    idx = {
        "G_O6": gm[1], "G_N1": gm[2], "G_H1": h1, "G_N2": gm[4], "G_H2": h2,
        "C_N4": ng + cm[5], "C_H4": ng + h4, "C_N3": ng + cm[6], "C_O2": ng + cm[7],
    }
    return atoms, idx
