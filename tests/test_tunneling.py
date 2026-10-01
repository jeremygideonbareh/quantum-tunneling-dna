import numpy as np

from qtdna.tunneling import HBAR_MW, skodje_truhlar, wigner, wkb_kappa, wkb_transmission


def test_wkb_exact_for_parabolic_barrier():
    v0, k = 0.5, 2.0  # eV, eV/(amu*A^2)
    s = np.linspace(-1.0, 1.0, 2001)
    v = v0 - 0.5 * k * s**2
    omega = np.sqrt(k)
    e = np.array([0.40, 0.45, 0.49])  # well inside the truncated region
    exact = 1.0 / (1.0 + np.exp(2 * np.pi * (v0 - e) / (HBAR_MW * omega)))
    assert np.allclose(wkb_transmission(s, v, e), exact, rtol=2e-2)


def test_corrections_agree_when_tunnelling_is_weak():
    w = wigner(300.0)
    st = skodje_truhlar(300.0, v_forward=0.5, delta_e=0.0)
    assert 1.0 < w < 1.1
    assert abs(st - w) < 0.02


def test_wkb_kappa_classical_limit_and_tunnelling():
    s = np.linspace(-3.0, 3.0, 601)
    # broad, heavy barrier: tunnelling negligible -> kappa ~ 1
    broad = 0.5 * np.exp(-(s / 2.0) ** 2)
    assert 0.9 < wkb_kappa(s, broad) < 1.3
    # thin barrier: strong tunnelling -> kappa >> 1
    thin = 0.5 * np.exp(-(s / 0.1) ** 2)
    assert wkb_kappa(s, thin) > 10.0
