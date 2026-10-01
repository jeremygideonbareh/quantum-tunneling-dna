from pathlib import Path

import pytest

from qtdna.scan import central_pair, frozen_indices, proton_triples, read_pdb

STRUCT = Path(__file__).resolve().parents[1] / "structures" / "k3"


@pytest.mark.parametrize("ctx, expect", [
    ("ACG", [("N1", "N3"), ("N4", "O6")]),   # G N1-H -> C N3, C N4-H -> G O6
    ("GTA", [("N3", "N1"), ("N6", "O4")]),   # T N3-H -> A N1, A N6-H -> T O4
])
def test_proton_triples(ctx, expect):
    cl = read_pdb(STRUCT / f"{ctx}.pdb")
    triples = proton_triples(cl, 3, ctx[1])
    got = [(cl.names[d], cl.names[a]) for d, _, a in triples]
    assert got == expect
    for d, h, a in triples:
        assert cl.atoms[h].symbol == "H"
        assert cl.atoms.get_distance(d, h) < 1.15            # covalent N-H
        assert 1.6 < cl.atoms.get_distance(h, a) < 2.3       # Watson-Crick H-bond


def test_flanks_and_caps_frozen_central_free():
    cl = read_pdb(STRUCT / "ACG.pdb")
    pyr, pur = central_pair(cl, 3)
    frozen = set(frozen_indices(cl, 3))
    caps = [n for n, name in enumerate(cl.names) if name == "C1'"]
    assert set(caps) <= frozen
    assert not ({n for n in pyr + pur if cl.names[n] != "C1'"} & frozen)
    assert {cl.names[n] for n in pyr if cl.names[n] != "C1'"} >= {"N3", "N4", "O2"}
