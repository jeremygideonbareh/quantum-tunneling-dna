"""Sequence-context bookkeeping.

Convention (same as COSMIC SBS-96): every context is written on the strand
whose centre base is a pyrimidine (C or T). So a G:C pair is labelled by its
C, an A:T pair by its T, and the 5'/3' neighbours are read on that strand.

    16 neighbour pairs x 2 centres = 32 trinucleotide contexts
    256 neighbour pairs x 2 centres = 512 pentanucleotide contexts
"""

from itertools import product

BASES = "ACGT"
PYRIMIDINES = "CT"
_COMP = str.maketrans("ACGT", "TGCA")

# The two transition classes the tunnelling hypothesis predicts.
# G*:C* tautomer pair -> G:T / A:C mispairs -> C>T (G:C>A:T)
# A*:T* tautomer pair -> A:C / G:T mispairs -> T>C (A:T>G:C)
TRANSITION = {"C": "T", "T": "C"}


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def canonical(kmer: str) -> str:
    """Return the pyrimidine-centred orientation of an odd-length k-mer."""
    kmer = kmer.upper()
    if len(kmer) % 2 == 0:
        raise ValueError("k-mer must have odd length")
    centre = kmer[len(kmer) // 2]
    return kmer if centre in PYRIMIDINES else revcomp(kmer)


def contexts(k: int = 3) -> list[str]:
    """All pyrimidine-centred k-mers (k odd), sorted like COSMIC."""
    if k % 2 == 0 or k < 1:
        raise ValueError("k must be a positive odd integer")
    flank = (k - 1) // 2
    out = []
    for centre in PYRIMIDINES:
        for left in product(BASES, repeat=flank):
            for right in product(BASES, repeat=flank):
                out.append("".join(left) + centre + "".join(right))
    return sorted(out, key=lambda s: (s[flank], s))


def is_cpg(kmer: str) -> bool:
    """True if the centre C is followed by G (or centre G preceded by C).

    CpG C>T is dominated by methyl-C deamination, not replication error,
    so these contexts need separate handling (see docs/research_plan.md).
    """
    kmer = canonical(kmer)
    mid = len(kmer) // 2
    return kmer[mid] == "C" and kmer[mid + 1] == "G"


def sbs96_label(kmer3: str) -> str:
    """COSMIC-style label for the transition at a trinucleotide, e.g. A[C>T]G."""
    kmer3 = canonical(kmer3)
    if len(kmer3) != 3:
        raise ValueError("need a trinucleotide")
    c = kmer3[1]
    return f"{kmer3[0]}[{c}>{TRANSITION[c]}]{kmer3[2]}"
