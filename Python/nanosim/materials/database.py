"""Material models for NanoFilterSim.

Complex index is *stored/returned* as n + ik with k >= 0 (database convention).
The solver converts to its internal e^{+i w t} form N = n - ik (see core/tmm.py).

Wavelengths are in nm everywhere. Extrapolation outside a material's validity range
is refused (ValueError) unless ``allow_extrapolation=True`` is passed explicitly.
"""
from __future__ import annotations
import numpy as np

REQUIRED_TABLE_KEYS = ("material", "wavelength_unit", "source", "reference", "limitations")


class Material:
    """Base class. Subclasses implement ``_nk`` and set ``range_nm``."""
    name = "material"
    meta: dict = {}
    range_nm = (0.0, np.inf)

    def __init__(self, allow_extrapolation: bool = False):
        self.allow_extrapolation = allow_extrapolation

    def nk(self, wavelength_nm):
        wl = np.atleast_1d(np.asarray(wavelength_nm, dtype=float))
        if not np.all(np.isfinite(wl)) or np.any(wl <= 0):
            raise ValueError("wavelength must be finite and > 0 nm")
        lo, hi = self.range_nm
        if (wl.min() < lo or wl.max() > hi) and not self.allow_extrapolation:
            raise ValueError(
                f"{self.name}: wavelength outside validity range [{lo:g}, {hi:g}] nm "
                "(set allow_extrapolation=True to override)")
        return self._nk(wl)


class ConstantMaterial(Material):
    def __init__(self, n, k=0.0, name=None, **kw):
        super().__init__(**kw)
        if not (np.isfinite(n) and n > 0 and np.isfinite(k) and k >= 0):
            raise ValueError("constant material needs n > 0 and k >= 0")
        self.n, self.k = float(n), float(k)
        self.name = name or f"const({n:g},{k:g})"
        self.meta = {"source": "user-specified constant", "limitations": "no dispersion"}

    def _nk(self, wl):
        return np.full(wl.shape, self.n + 1j * self.k)


class DispersionMaterial(Material):
    """n^2 = A + sum_i B_i l^2/(l^2 - C_i) + sum_j D_j/(l^2 - E_j),  l in micrometres.
    k = 0 (dataset gives no extinction coefficient)."""

    def __init__(self, name, A, lam_terms, plain_terms, range_nm, meta, **kw):
        super().__init__(**kw)
        self.name, self.A = name, A
        self.lam_terms, self.plain_terms = lam_terms, plain_terms
        self.range_nm, self.meta = range_nm, meta

    def _nk(self, wl):
        l2 = (wl * 1e-3) ** 2
        n2 = np.full(wl.shape, float(self.A))
        for B, C in self.lam_terms:
            n2 = n2 + B * l2 / (l2 - C)
        for D, E in self.plain_terms:
            n2 = n2 + D / (l2 - E)
        return np.sqrt(n2).astype(complex)


class TabulatedMaterial(Material):
    """Linear interpolation of tabulated n, k. Outside the table: refused, or
    (explicit opt-in) nearest-edge value held constant."""

    def __init__(self, wl_nm, n, k, name, meta, **kw):
        super().__init__(**kw)
        wl_nm, n, k = (np.asarray(a, float) for a in (wl_nm, n, k))
        if wl_nm.size < 2 or np.any(np.diff(wl_nm) <= 0):
            raise ValueError("table wavelengths must be strictly increasing (>= 2 points)")
        if np.any(n <= 0) or np.any(k < 0):
            raise ValueError("table needs n > 0 and k >= 0")
        self.wl, self.n, self.k, self.name, self.meta = wl_nm, n, k, name, meta
        self.range_nm = (wl_nm[0], wl_nm[-1])

    def _nk(self, wl):
        return np.interp(wl, self.wl, self.n) + 1j * np.interp(wl, self.wl, self.k)


def load_table(path, allow_extrapolation=False):
    """CSV with '# key: value' metadata lines, then columns wavelength,n[,k]."""
    meta, rows = {}, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                if ":" in s:
                    key, val = s[1:].split(":", 1)
                    meta[key.strip().lower()] = val.strip()
                continue
            try:
                rows.append([float(x) for x in s.split(",")])
            except ValueError:
                continue  # header line
    missing = [k for k in REQUIRED_TABLE_KEYS if not meta.get(k)]
    if missing:
        raise ValueError(f"{path}: missing required metadata: {missing}")
    unit = meta["wavelength_unit"].lower()
    scale = {"nm": 1.0, "um": 1e3, "µm": 1e3}.get(unit)
    if scale is None:
        raise ValueError(f"unsupported wavelength_unit '{unit}' (use nm or um)")
    a = np.array(rows)
    if a.ndim != 2 or a.shape[1] < 2:
        raise ValueError(f"{path}: need columns wavelength,n[,k]")
    k = a[:, 2] if a.shape[1] > 2 else np.zeros(len(a))
    return TabulatedMaterial(a[:, 0] * scale, a[:, 1], k, meta["material"], meta,
                             allow_extrapolation=allow_extrapolation)


# --- built-in dispersion models -------------------------------------------------
# VERIFY BEFORE PUBLICATION: coefficients were entered from the literature values
# (see 'source'/'reference'), without network access in the build environment.
def _sio2(**kw):
    return DispersionMaterial(
        "SiO2 (fused silica)", 1.0,
        [(0.6961663, 0.0684043 ** 2), (0.4079426, 0.1162414 ** 2), (0.8974794, 9.896161 ** 2)],
        [], (210.0, 3710.0),
        {"form": "bulk fused silica", "wavelength_unit": "nm",
         "source": "Sellmeier fit (refractiveindex.info entry 'Malitson 1965', 20 C)",
         "reference": "I. H. Malitson, J. Opt. Soc. Am. 55, 1205 (1965)",
         "limitations": "bulk glass, not a sputtered/evaporated SiO2 film; k assumed 0; "
                        "valid 0.21-3.71 um; coefficients unverified against the database"}, **kw)


def _tio2(**kw):
    return DispersionMaterial(
        "TiO2 (rutile, ordinary ray)", 5.913, [], [(0.2441, 0.0803)], (430.0, 1530.0),
        {"form": "bulk rutile single crystal, ordinary ray", "wavelength_unit": "nm",
         "source": "refractiveindex.info entry 'DeVore 1951 (n_o)'",
         "reference": "J. R. DeVore, J. Opt. Soc. Am. 41, 416 (1951)",
         "limitations": "bulk crystal index (~2.65 at 550 nm) is HIGHER than typical TiO2 thin films "
                        "(anatase/amorphous, ~2.3-2.5); k assumed 0 although real TiO2 absorbs near "
                        "the UV edge; valid 0.43-1.53 um only; coefficients unverified against database"}, **kw)


def get_material(spec: str, allow_extrapolation: bool = False) -> Material:
    """spec: 'air' | 'const:n' | 'const:n,k' | 'SiO2' | 'fused_silica' | 'TiO2' | 'table:<path>'."""
    s = spec.strip()
    low = s.lower()
    kw = {"allow_extrapolation": allow_extrapolation}
    if low == "air":
        return ConstantMaterial(1.0, 0.0, name="air", **kw)
    if low.startswith("const:"):
        parts = [float(x) for x in s[6:].split(",")]
        if len(parts) not in (1, 2):
            raise ValueError(f"bad constant material spec '{spec}'")
        return ConstantMaterial(parts[0], parts[1] if len(parts) == 2 else 0.0, **kw)
    if low in ("sio2", "fused_silica"):
        return _sio2(**kw)
    if low == "tio2":
        return _tio2(**kw)
    if low.startswith("table:"):
        return load_table(s[6:], allow_extrapolation)
    raise ValueError(f"unknown material '{spec}'")
