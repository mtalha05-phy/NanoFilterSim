"""Plain-text case files shared by the Python, MATLAB and C++ implementations.

    name <str>
    incident <material spec>          (lossless)
    substrate <material spec>
    polarization TE|TM
    theta <deg>
    wavelength <min_nm> <max_nm> <npoints>
    layer <material spec> <thickness_nm>     (repeat; incident side first)
Lines beginning with '#' are comments.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from ..materials.database import get_material
from .tmm import solve_structure


@dataclass
class Case:
    name: str = "case"
    incident: str = "air"
    substrate: str = "SiO2"
    polarization: str = "TE"
    theta: float = 0.0
    wl: tuple = (400.0, 800.0, 401)
    layers: list = field(default_factory=list)


def read_case(path) -> Case:
    c = Case()
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.split("#", 1)[0].strip() if raw.lstrip().startswith("#") else raw.strip()
            if not line:
                continue
            key, _, rest = line.partition(" ")
            rest = rest.strip()
            if key == "name": c.name = rest
            elif key == "incident": c.incident = rest
            elif key == "substrate": c.substrate = rest
            elif key == "polarization": c.polarization = rest.upper()
            elif key == "theta": c.theta = float(rest)
            elif key == "wavelength":
                a, b, n = rest.split()
                c.wl = (float(a), float(b), int(n))
            elif key == "layer":
                spec, d = rest.rsplit(None, 1)
                c.layers.append((spec, float(d)))
            else:
                raise ValueError(f"unknown key '{key}' in {path}")
    return c


def wavelength_grid(wl):
    a, b, n = wl
    if n < 1:
        raise ValueError("npoints must be >= 1")
    if n == 1:
        return np.array([a])
    return a + np.arange(n) * ((b - a) / (n - 1))


def run_case(c: Case, allow_extrapolation=False):
    wl = wavelength_grid(c.wl)
    inc = get_material(c.incident, allow_extrapolation)
    sub = get_material(c.substrate, allow_extrapolation)
    layers = [(get_material(s, allow_extrapolation), d) for s, d in c.layers]
    return solve_structure(wl, inc, layers, sub, c.theta, c.polarization)


def write_csv(path, res):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("wavelength_nm,R,T,A\n")
        for w, R, T, A in zip(res.wavelength_nm, res.R, res.T, res.A):
            fh.write(f"{w:.17g},{R:.17g},{T:.17g},{A:.17g}\n")
