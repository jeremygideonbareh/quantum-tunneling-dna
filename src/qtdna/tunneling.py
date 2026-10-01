"""Tunnelling corrections from a 1D reaction profile.

Three levels, cheapest first. All energies in eV, temperatures in K.

* Wigner      - needs only the imaginary barrier frequency; valid for weak tunnelling.
* Skodje-Truhlar (parabolic, truncated) - frequency + barrier height.
* Semiclassical WKB on the scanned profile - uses the actual barrier shape
  along the mass-weighted path. This is the number to put in the fingerprint;
  the other two are sanity checks.

None of these include environment-assisted (open quantum system) effects
treated by Slocombe et al. 2022; we state that as a limitation.
"""

import numpy as np

KB = 8.617333262e-5          # eV/K
HBAR_MW = 0.064654           # hbar in sqrt(eV*amu)*Angstrom
CM1_TO_EV = 1.239841984e-4


def wigner(imag_freq_cm: float, T: float = 310.0) -> float:
    u = imag_freq_cm * CM1_TO_EV / (KB * T)
    return 1.0 + u * u / 24.0


def skodje_truhlar(imag_freq_cm: float, v_forward: float, delta_e: float,
                   T: float = 310.0) -> float:
    """Skodje & Truhlar, J. Phys. Chem. 85, 624 (1981)."""
    alpha = 2.0 * np.pi / (imag_freq_cm * CM1_TO_EV)
    beta = 1.0 / (KB * T)
    v = v_forward if delta_e <= 0 else v_forward - delta_e
    if beta < alpha:
        x = beta * np.pi / alpha
        return x / np.sin(x) - beta / (alpha - beta) * np.exp((beta - alpha) * v)
    return beta / (beta - alpha) * (np.exp((beta - alpha) * v) - 1.0)


def mass_weighted_path(images, masses) -> np.ndarray:
    """Cumulative mass-weighted arc length s (sqrt(amu)*Angstrom) along images."""
    s = [0.0]
    m = np.asarray(masses)[:, None]
    for a, b in zip(images[:-1], images[1:]):
        s.append(s[-1] + float(np.sqrt((m * (b - a) ** 2).sum())))
    return np.asarray(s)


def wkb_transmission(s: np.ndarray, v: np.ndarray, energies: np.ndarray) -> np.ndarray:
    """P(E) = 1 / (1 + exp(2*theta)), theta = int sqrt(2(V-E)) ds / hbar.

    Above the barrier top, theta is continued as -theta(2*Vmax - E) (the
    standard reflection used with the Kemble formula), so P -> 1 rather
    than sticking at 1/2.
    """
    s_fine = np.linspace(s[0], s[-1], 4000)
    v_fine = np.interp(s_fine, s, v)
    vmax = v_fine.max()

    def theta(e):
        integrand = np.sqrt(np.clip(2.0 * (v_fine - e), 0.0, None))
        return np.trapezoid(integrand, s_fine) / HBAR_MW

    out = np.empty_like(energies, dtype=float)
    for i, e in enumerate(energies):
        th = theta(e) if e <= vmax else -theta(2.0 * vmax - e)
        out[i] = 1.0 / (1.0 + np.exp(np.clip(2.0 * th, -700.0, 700.0)))
    return out


def wkb_kappa(s: np.ndarray, v: np.ndarray, T: float = 310.0, n: int = 600) -> float:
    """Thermal tunnelling factor: quantum rate / classical over-barrier rate.

    v is referenced to the reactant minimum (v[0] = 0). Energies start at the
    endpoint that lies higher, so both sides are classically allowed.
    """
    v = np.asarray(v) - v[0]
    vmax = v.max()
    e0 = max(v[0], v[-1])
    beta = 1.0 / (KB * T)
    energies = np.linspace(e0, vmax + 20.0 / beta, n)
    p = wkb_transmission(s, v, energies)
    quantum = np.trapezoid(p * np.exp(-beta * energies), energies)
    classical = np.exp(-beta * vmax) / beta
    return float(quantum / classical)
