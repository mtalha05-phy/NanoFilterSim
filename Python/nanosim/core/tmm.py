"""Transfer Matrix Method core (Python reference implementation).

CONVENTIONS (identical in MATLAB and C++):
  * Time dependence e^{+i w t} (Macleod). Layer matrix exactly as in the project spec:
        M_j = [[cos d, i sin d / q], [i q sin d, cos d]]
  * Database index is n + ik (k >= 0). The solver uses N = n - ik internally, which is the
    passive-medium form of the e^{+i w t} convention. (R, T, A are identical to the
    e^{-i w t}, n+ik formulation with the conjugate matrix.)
  * Tangential wavevector conserved: n0 sin(th0) = N_j sin(th_j).
        eta_j = N_j cos(th_j) = sqrt(N_j^2 - (n0 sin th0)^2),
    branch chosen so Im(eta_j) <= 0 (decaying / outgoing wave in passive media;
    evanescent lossless waves get eta = -i|eta|).
  * q_TE = eta,  q_TM = N^2/eta (= N / cos th).   delta_j = 2 pi eta_j d_j / lambda.
  * Layer 1 is adjacent to the incident medium; M = M_1 M_2 ... M_N; substrate behind layer N.
  * Incident medium must be lossless (real index); angle th0 in [0, 90) degrees.
  * No clipping of R, T, A is ever performed.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class TMMResult:
    wavelength_nm: np.ndarray
    R: np.ndarray
    T: np.ndarray
    A: np.ndarray
    r: np.ndarray
    t: np.ndarray
    theta_deg: float
    polarization: str


def _wavelengths(wl):
    wl = np.atleast_1d(np.asarray(wl, dtype=float))
    if wl.ndim != 1 or wl.size == 0:
        raise ValueError("wavelength array must be 1-D and non-empty")
    if not np.all(np.isfinite(wl)) or np.any(wl <= 0):
        raise ValueError("wavelengths must be finite and > 0 (nm)")
    return wl


def internal_index(nk):
    """Database n+ik -> solver N = n-ik."""
    return np.conjugate(np.asarray(nk, dtype=complex))


def normal_eta(N, s):
    """eta = N cos(theta) with the passive-medium branch (Im eta <= 0)."""
    eta = np.sqrt(np.asarray(N * N - s * s, dtype=complex))
    return np.where(eta.imag > 0, -eta, eta)


def admittance(N, eta, pol):
    pol = pol.upper()
    if np.any(np.abs(eta) < 1e-14):
        raise ValueError("zero normal wavevector component (grazing/critical): admittance undefined")
    if pol == "TE":
        return eta
    if pol == "TM":
        return N * N / eta
    raise ValueError("polarization must be 'TE' or 'TM'")


def phase_thickness(wl_nm, eta, d_nm):
    return 2.0 * np.pi * eta * d_nm / wl_nm


def layer_matrix(delta, q):
    c, s = np.cos(delta), np.sin(delta)
    M = np.empty(np.shape(delta) + (2, 2), dtype=complex)
    M[..., 0, 0] = c
    M[..., 0, 1] = 1j * s / q
    M[..., 1, 0] = 1j * q * s
    M[..., 1, 1] = c
    return M


def solve(wavelength_nm, nk_inc, nk_layers, d_nm, nk_sub, theta_deg=0.0, pol="TE"):
    """Low-level solver. nk_* are complex (n+ik) arrays on the wavelength grid
    (nk_layers is a list of arrays, ordered incident side -> substrate side)."""
    wl = _wavelengths(wavelength_nm)
    pol = pol.upper()
    if pol not in ("TE", "TM"):
        raise ValueError("polarization must be 'TE' or 'TM'")
    if not (np.isfinite(theta_deg) and 0.0 <= theta_deg < 90.0):
        raise ValueError("incidence angle must satisfy 0 <= theta < 90 deg")
    d_nm = list(d_nm)
    if len(d_nm) != len(nk_layers):
        raise ValueError("number of thicknesses must equal number of layers")
    for d in d_nm:
        if not np.isfinite(d) or d < 0:
            raise ValueError("layer thickness must be finite and >= 0 (nm)")

    n0c = np.asarray(nk_inc, dtype=complex)
    if np.any(np.abs(n0c.imag) > 1e-12):
        raise ValueError("incident medium must be lossless (real index)")
    n0 = n0c.real
    s = n0 * np.sin(np.deg2rad(theta_deg))
    eta0 = normal_eta(n0.astype(complex), s)
    q0 = admittance(n0.astype(complex), eta0, pol)

    M = np.zeros(wl.shape + (2, 2), dtype=complex)
    M[..., 0, 0] = 1.0
    M[..., 1, 1] = 1.0
    for nk, d in zip(nk_layers, d_nm):          # M = M_1 M_2 ... M_N
        N = internal_index(nk)
        eta = normal_eta(N, s)
        q = admittance(N, eta, pol)
        M = M @ layer_matrix(phase_thickness(wl, eta, d), q)

    Ns = internal_index(nk_sub)
    qs = admittance(Ns, normal_eta(Ns, s), pol)

    B = M[..., 0, 0] + M[..., 0, 1] * qs
    C = M[..., 1, 0] + M[..., 1, 1] * qs
    den = q0 * B + C
    if not (np.all(np.isfinite(den)) and np.all(np.abs(den) > 0)):
        raise FloatingPointError("overflow in transfer matrix (layer too absorbing/thick)")
    r = (q0 * B - C) / den
    t = 2.0 * q0 / den
    R = np.abs(r) ** 2
    T = (qs.real / q0.real) * np.abs(t) ** 2
    A = 1.0 - R - T
    return TMMResult(wl, R, T, A, r, t, float(theta_deg), pol)


def solve_structure(wavelength_nm, incident, layers, substrate, theta_deg=0.0, pol="TE"):
    """layers: list of (Material, thickness_nm)."""
    wl = _wavelengths(wavelength_nm)
    return solve(wl, incident.nk(wl), [m.nk(wl) for m, _ in layers], [d for _, d in layers],
                 substrate.nk(wl), theta_deg, pol)


def energy_residual(res: TMMResult):
    """R + T - 1 (zero for lossless structures; equals -A in general)."""
    return res.R + res.T - 1.0
