import pandas as pd
import pytest

from qtdna.mutations import count_transitions, kmer_opportunities


GENOME = {"chr": "AACGTTACAGG"}
#               1234567890a   (1-based)


def test_opportunities_collapse_strands():
    sites = kmer_opportunities(GENOME, 3)
    assert sum(sites.values()) == len(GENOME["chr"]) - 2
    # ACG at pos 2-4, and CGT (revcomp ACG) at pos 3-5
    assert sites["ACG"] == 2


def test_transitions_counted_on_pyrimidine_strand():
    muts = pd.DataFrame({
        "chrom": ["chr", "chr", "chr"],
        "pos": [3, 4, 8],          # C>T in ACG; G>A in CGT (= C>T in ACG); C>A (not a transition)
        "ref": ["C", "G", "C"],
        "alt": ["T", "A", "A"],
    })
    df = count_transitions(muts, GENOME, 3).set_index("context")
    assert df.loc["ACG", "n_mut"] == 2
    assert df["n_mut"].sum() == 2
    assert df.loc["ACG", "rate"] == 1.0


def test_reference_mismatch_raises():
    muts = pd.DataFrame({"chrom": ["chr"], "pos": [5], "ref": ["G"], "alt": ["A"]})
    with pytest.raises(ValueError):
        count_transitions(muts, GENOME, 3)
