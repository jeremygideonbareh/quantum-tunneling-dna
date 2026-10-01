"""Duplex stability baseline: SantaLucia (1998) unified nearest-neighbour dG37.

Replication-error hotspots are enriched in low-stability DNA (Hasenauer et al.
NAR 2025), so every fingerprint test must control for local stability.
"""

from .contexts import revcomp

# kcal/mol at 37 C, 1 M NaCl; SantaLucia, PNAS 95:1460 (1998), Table 2
_NN = {"AA": -1.00, "AT": -0.88, "TA": -0.58, "CA": -1.45, "GT": -1.44,
       "CT": -1.28, "GA": -1.30, "CG": -2.17, "GC": -2.24, "GG": -1.84}
NN_DG37 = {**_NN, **{revcomp(k): v for k, v in _NN.items()}}


def stacking_dg(seq: str) -> float:
    """Sum of nearest-neighbour dG37 over all steps in seq (more negative = more stable)."""
    seq = seq.upper()
    return sum(NN_DG37[seq[i:i + 2]] for i in range(len(seq) - 1))
