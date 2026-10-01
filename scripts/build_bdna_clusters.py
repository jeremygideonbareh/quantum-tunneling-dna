"""Build methyl-capped B-DNA base-pair stacks for every context (PyMOL fnab).

For each pyrimidine-centred k-mer, build ideal double-stranded B-DNA, keep
only the bases plus C1' (which becomes a methyl cap), add hydrogens and
write structures/k{k}/<context>.pdb. Strand A carries the context as written;
the central pair is (pyrimidine on A) : (purine on B).

PyMOL's open-source wheel needs numpy<2, so run this in its own venv:
    python -m venv .pymolenv && .pymolenv/bin/pip install pymol-open-source-whl "numpy<2"
    .pymolenv/bin/python scripts/build_bdna_clusters.py --k 3
"""

import argparse
import sys
from itertools import product
from pathlib import Path

from pymol import cmd

BACKBONE = "P+OP1+OP2+O1P+O2P+O5'+C5'+C4'+O4'+C3'+O3'+C2'+O2'"


def contexts(k):
    flank = (k - 1) // 2
    out = []
    for centre in "CT":
        for left in product("ACGT", repeat=flank):
            for right in product("ACGT", repeat=flank):
                out.append("".join(left) + centre + "".join(right))
    return out


def build(seq, path):
    cmd.reinitialize()
    cmd.fnab(seq, name="dna", mode="DNA", form="B", dbl_helix=1)
    cmd.remove(f"dna and name {BACKBONE}")
    cmd.h_add("dna")
    cmd.alter("dna", "type='ATOM'")
    cmd.save(str(path), "dna")
    return cmd.count_atoms("dna")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--out", default="structures")
    ap.add_argument("contexts", nargs="*")
    args = ap.parse_args()
    out = Path(args.out) / f"k{args.k}"
    out.mkdir(parents=True, exist_ok=True)
    for ctx in args.contexts or contexts(args.k):
        n = build(ctx, out / f"{ctx}.pdb")
        print(f"{ctx}: {n} atoms", file=sys.stderr)


if __name__ == "__main__":
    main()
