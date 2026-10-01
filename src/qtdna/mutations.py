"""Count transitions per sequence context and normalise by opportunity.

Input: a table of single-base substitutions (chrom, pos 1-based, ref, alt)
and the reference genome. Output, per pyrimidine-centred k-mer:

    n_mut      transitions observed (C>T for C-centred, T>C for T-centred)
    n_sites    times that k-mer occurs in the genome (both strands collapsed)
    rate       n_mut / n_sites  (per-site mutability, the thing to correlate)

Normalising by n_sites matters: raw counts mostly track how common a
k-mer is, not how mutable it is.
"""

from collections import Counter
from pathlib import Path

import pandas as pd

from .contexts import PYRIMIDINES, TRANSITION, canonical, contexts, revcomp


def read_fasta(path: str | Path) -> dict[str, str]:
    seqs, name, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def kmer_opportunities(genome: dict[str, str], k: int = 3) -> Counter:
    flank = k // 2
    counts = Counter()
    for seq in genome.values():
        for i in range(flank, len(seq) - flank):
            kmer = seq[i - flank:i + flank + 1]
            if "N" not in kmer:
                counts[canonical(kmer)] += 1
    return counts


def count_transitions(muts: pd.DataFrame, genome: dict[str, str], k: int = 3) -> pd.DataFrame:
    """muts needs columns chrom, pos (1-based), ref, alt."""
    flank = k // 2
    observed, mismatched_ref = Counter(), 0
    for chrom, pos, ref, alt in muts[["chrom", "pos", "ref", "alt"]].itertuples(index=False):
        seq = genome[str(chrom)]
        i = int(pos) - 1
        if i - flank < 0 or i + flank >= len(seq):
            continue
        kmer = seq[i - flank:i + flank + 1]
        if kmer[flank] != ref.upper():
            mismatched_ref += 1
            continue
        if ref.upper() not in PYRIMIDINES:
            kmer, alt = revcomp(kmer), revcomp(alt.upper())
        if TRANSITION[kmer[flank]] == alt.upper():
            observed[kmer] += 1
    if mismatched_ref:
        raise ValueError(f"{mismatched_ref} mutations disagree with the reference "
                         "(wrong genome build or 0- vs 1-based positions?)")

    sites = kmer_opportunities(genome, k)
    rows = [{"context": c, "n_mut": observed[c], "n_sites": sites[c]} for c in contexts(k)]
    df = pd.DataFrame(rows)
    df["rate"] = df["n_mut"] / df["n_sites"].where(df["n_sites"] > 0)
    return df
