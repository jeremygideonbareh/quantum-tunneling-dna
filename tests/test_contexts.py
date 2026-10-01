from qtdna.contexts import canonical, contexts, is_cpg, revcomp, sbs96_label


def test_counts():
    assert len(contexts(3)) == 32
    assert len(contexts(5)) == 512
    assert len(set(contexts(5))) == 512


def test_all_pyrimidine_centred():
    assert all(c[1] in "CT" for c in contexts(3))


def test_canonical_flips_purine_centre():
    assert canonical("AGT") == "ACT"
    assert canonical("ACT") == "ACT"
    assert canonical(revcomp("TTCAG")) == "TTCAG"


def test_cpg():
    assert is_cpg("ACG")
    assert is_cpg("CGT")  # revcomp ACG, centre G preceded by C
    assert not is_cpg("ACA")
    assert sum(is_cpg(c) for c in contexts(3)) == 4


def test_sbs96_label():
    assert sbs96_label("ACG") == "A[C>T]G"
    assert sbs96_label("TAC") == "G[T>C]A"
