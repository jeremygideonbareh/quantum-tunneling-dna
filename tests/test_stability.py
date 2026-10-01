from qtdna.contexts import revcomp
from qtdna.stability import NN_DG37, stacking_dg


def test_all_16_steps_and_strand_symmetry():
    assert len(NN_DG37) == 16
    for step, dg in NN_DG37.items():
        assert NN_DG37[revcomp(step)] == dg


def test_gc_more_stable_than_at():
    assert stacking_dg("GCG") < stacking_dg("ATA")
    assert stacking_dg("ACG") == stacking_dg(revcomp("ACG"))
