"""Phase-2 validation tests for the Python TMM core.
Runs under pytest if available, or via tests/run_tests.py (no external dependency)."""
import os, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from nanosim.core import tmm
from nanosim.core.tmm import solve, solve_structure, normal_eta, admittance, phase_thickness, layer_matrix
from nanosim.core.case import read_case, run_case, write_csv
from nanosim.materials.database import get_material, load_table, ConstantMaterial

WL = np.linspace(450, 800, 141)
air, sio2, tio2 = get_material("air"), get_material("SiO2"), get_material("TiO2")


def raises(exc, fn, *a, **k):
    try:
        fn(*a, **k)
    except exc:
        return True
    raise AssertionError(f"{exc.__name__} not raised")


def fresnel(n0, ns, th0, pol):
    s = n0 * np.sin(np.deg2rad(th0))
    c0 = np.sqrt(1 - (s / n0) ** 2 + 0j)
    cs = np.sqrt(1 - (s / ns) ** 2 + 0j)
    if pol == "TE":
        return (n0 * c0 - ns * cs) / (n0 * c0 + ns * cs)
    return (ns * c0 - n0 * cs) / (ns * c0 + n0 * cs)


# ---------------- mathematical ----------------
def test_phase_thickness():
    d = phase_thickness(550.0, 1.5 + 0j, 110.0)
    assert abs(d - 2 * np.pi * 1.5 * 110 / 550) < 1e-14


def test_admittance_normal_and_oblique():
    N = np.array([2.0 + 0j]); s = 0.0
    assert abs(admittance(N, normal_eta(N, s), "TE")[0] - 2.0) < 1e-14
    assert abs(admittance(N, normal_eta(N, s), "TM")[0] - 2.0) < 1e-14
    s = 1.0 * np.sin(np.deg2rad(30)); eta = normal_eta(N, s)   # eta = N cos(th)
    cos_t = np.sqrt(1 - (s / 2.0) ** 2)
    assert abs(admittance(N, eta, "TE")[0] - 2 * cos_t) < 1e-14
    assert abs(admittance(N, eta, "TM")[0] - 2 / cos_t) < 1e-14


def test_layer_matrix_unimodular_and_quarter_wave():
    q = 1.7 + 0j
    M = layer_matrix(np.array([0.7]), q)[0]
    assert abs(np.linalg.det(M) - 1) < 1e-14                      # lossless: det = 1
    Mq = layer_matrix(np.array([np.pi / 2]), q)[0]
    assert np.allclose(Mq, [[0, 1j / q], [1j * q, 0]], atol=1e-14)


def test_matrix_order_matters_and_is_incident_first():
    # two different layers: reversing the stack must change the (asymmetric) result
    kw = dict(theta_deg=0.0, pol="TE")
    wl = np.array([600.0])
    n = lambda v: np.array([v + 0j])
    a = solve(wl, n(1.0), [n(2.4), n(1.46)], [60, 120], n(1.7), **kw)
    b = solve(wl, n(1.0), [n(1.46), n(2.4)], [120, 60], n(1.7), **kw)
    assert abs(a.R[0] - b.R[0]) > 1e-3
    # and the manual product M1 M2 agrees with solve()
    q = lambda v: v + 0j
    M = layer_matrix(phase_thickness(600.0, q(2.4), 60), q(2.4)) @ \
        layer_matrix(phase_thickness(600.0, q(1.46), 120), q(1.46))
    B = M[0, 0] + M[0, 1] * 1.7; C = M[1, 0] + M[1, 1] * 1.7
    r = (B - C) / (B + C)
    assert abs(abs(r) ** 2 - a.R[0]) < 1e-14


# ---------------- physical ----------------
def test_single_interface_fresnel_all_cases():
    for pol in ("TE", "TM"):
        for th in (0, 30, 60, 85):
            res = solve_structure(WL, air, [], sio2, th, pol)
            ns = sio2.nk(WL).real
            ref = np.array([abs(fresnel(1.0, n, th, pol)) ** 2 for n in ns])
            assert np.max(np.abs(res.R - ref)) < 1e-13
            assert np.max(np.abs(res.R + res.T - 1)) < 1e-13


def test_empty_layer_sequence_is_bare_interface():
    res = solve_structure(WL, air, [], ConstantMaterial(1.5), 0, "TE")
    assert np.allclose(res.R, 0.04, atol=1e-14)


def test_quarter_wave_antireflection_single_layer():
    ns = 1.52; nl = np.sqrt(ns); lam0 = 600.0
    wl = np.array([lam0])
    res = solve(wl, np.array([1 + 0j]), [np.array([nl + 0j])], [lam0 / (4 * nl)], np.array([ns + 0j]))
    assert res.R[0] < 1e-28


def test_half_wave_layer_is_absent():
    lam0 = 600.0; wl = np.array([lam0]); nl = 2.3
    res = solve(wl, np.array([1 + 0j]), [np.array([nl + 0j])], [lam0 / (2 * nl)], np.array([1.5 + 0j]))
    assert abs(res.R[0] - 0.04) < 1e-14


def test_energy_conservation_lossless_dbr_all_angles():
    L0 = 550.0
    dH, dL = L0 / (4 * tio2.nk(L0)[0].real), L0 / (4 * sio2.nk(L0)[0].real)
    layers = [(tio2, dH), (sio2, dL)] * 12
    for pol in ("TE", "TM"):
        for th in (0, 20, 45, 70):
            r = solve_structure(WL, air, layers, sio2, th, pol)
            assert np.max(np.abs(r.R + r.T - 1)) < 1e-12, (pol, th)


def test_absorptance_nonnegative_and_positive_for_k_gt_0():
    lay = [(ConstantMaterial(2.4, 0.05), 100.0)]
    for pol in ("TE", "TM"):
        r = solve_structure(WL, air, lay, sio2, 40, pol)
        assert r.A.min() > 1e-4          # convention check: loss, not gain
        assert np.all(r.R + r.T < 1)


def test_thick_absorber_limit_equals_interface_fresnel():
    N = 2.0 + 0.5j
    r = solve_structure(np.array([600.0]), air, [(ConstantMaterial(2.0, 0.5), 5000.0)], sio2, 0, "TE")
    ref = abs((1 - N) / (1 + N)) ** 2
    assert abs(r.R[0] - ref) < 1e-12 and r.T[0] < 1e-12


def test_total_internal_reflection_branch():
    prism, vac = ConstantMaterial(1.5), ConstantMaterial(1.0)
    for pol in ("TE", "TM"):
        r = solve_structure(WL, prism, [], vac, 60, pol)
        assert np.max(np.abs(r.R - 1)) < 1e-13 and np.max(np.abs(r.T)) < 1e-13


def test_frustrated_tir_lossless_and_decays_with_gap():
    prism = ConstantMaterial(1.5)
    for pol in ("TE", "TM"):
        prev = None
        for gap in (100.0, 400.0, 2000.0):
            r = solve_structure(WL, prism, [(air, gap)], prism, 50, pol)
            assert np.max(np.abs(r.R + r.T - 1)) < 1e-12
            assert np.all(r.T >= 0) and np.all(r.T < 1)
            if prev is not None:
                assert np.all(r.T < prev)          # evanescent coupling decreases with gap
            prev = r.T
        assert prev.max() < 1e-4                   # 2000 nm gap: essentially fully reflected


def test_quarter_wave_dbr_matches_analytic_peak():
    L0 = 550.0
    nH, nL, ns = tio2.nk(L0)[0].real, sio2.nk(L0)[0].real, sio2.nk(L0)[0].real
    dH, dL = L0 / (4 * nH), L0 / (4 * nL)
    for N in (1, 4, 8, 15):
        r = solve_structure(np.array([L0]), air, [(tio2, dH), (sio2, dL)] * N, sio2, 0, "TE")
        Y = (nH / nL) ** (2 * N) * ns            # exact input admittance at design wavelength
        Rref = ((1 - Y) / (1 + Y)) ** 2
        assert abs(r.R[0] - Rref) < 1e-12, N


def test_stopband_width_approximation_is_only_approximate():
    # informational: finite-N numerical width vs idealised (4/pi) asin formula, loose agreement
    L0 = 550.0
    nH, nL = tio2.nk(L0)[0].real, sio2.nk(L0)[0].real
    dH, dL = L0 / (4 * nH), L0 / (4 * nL)
    wl = np.linspace(450, 800, 3501)
    r = solve_structure(wl, air, [(tio2, dH), (sio2, dL)] * 15, sio2, 0, "TE")
    band = wl[r.R > 0.99]
    width_frac = (band.max() - band.min()) / L0
    approx = 4 / np.pi * np.arcsin((nH - nL) / (nH + nL))
    assert 0.5 * approx < width_frac < 1.5 * approx, (width_frac, approx)


def _fp_layers(L0=550.0, m=2, N=5):
    nH, nL = tio2.nk(L0)[0].real, sio2.nk(L0)[0].real
    dH, dL = L0 / (4 * nH), L0 / (4 * nL)
    mirror = [(tio2, dH), (sio2, dL)] * N
    return mirror + [(sio2, m * L0 / (2 * nL))] + mirror[::-1]


def test_cavity_resonance_symmetric_structure_reaches_unity():
    # incident medium = substrate = SiO2 -> structure is mirror-symmetric -> T_max = 1 at resonance
    wl = np.linspace(540, 560, 4001)
    r = solve_structure(wl, sio2, _fp_layers(), sio2, 0, "TE")
    assert r.T.max() > 0.9999 and abs(wl[np.argmax(r.T)] - 550.0) < 0.05


def test_cavity_resonance_asymmetric_air_to_substrate():
    # air-side and substrate-side mirrors differ -> resonance still at 550 nm but T_max < 1
    wl = np.linspace(540, 560, 4001)
    r = solve_structure(wl, air, _fp_layers(), sio2, 0, "TE")
    assert abs(wl[np.argmax(r.T)] - 550.0) < 0.05
    assert 0.5 < r.T.max() < 0.999
    off = solve_structure(np.array([500.0]), air, _fp_layers(), sio2, 0, "TE")
    assert off.T[0] < 0.5


# ---------------- software ----------------
def test_invalid_inputs():
    one = np.array([1 + 0j])
    raises(ValueError, solve, np.array([0.0]), one, [], [], one)
    raises(ValueError, solve, np.array([-5.0]), one, [], [], one)
    raises(ValueError, solve, np.array([np.nan]), one, [], [], one)
    raises(ValueError, solve, np.array([]), one, [], [], one)
    raises(ValueError, solve, np.array([500.0]), one, [one], [-10.0], one)
    raises(ValueError, solve, np.array([500.0]), one, [one], [np.nan], one)
    raises(ValueError, solve, np.array([500.0]), one, [one, one], [10.0], one)
    raises(ValueError, solve, np.array([500.0]), one, [], [], one, theta_deg=90.0)
    raises(ValueError, solve, np.array([500.0]), one, [], [], one, pol="XX")
    raises(ValueError, solve, np.array([500.0]), np.array([1 + 0.1j]), [], [], one)  # lossy incident


def test_material_range_and_definition_checks():
    raises(ValueError, tio2.nk, np.array([400.0]))               # below TiO2 validity range
    assert get_material("TiO2", allow_extrapolation=True).nk(np.array([400.0]))[0].real > 2.5
    raises(ValueError, get_material, "unobtainium")
    raises(ValueError, ConstantMaterial, -1.0)
    raises(ValueError, ConstantMaterial, 1.5, -0.1)
    raises(ValueError, get_material, "const:1,2,3")


def test_sellmeier_reference_values():
    assert abs(sio2.nk(587.5618)[0].real - 1.45846) < 2e-4       # fused silica n_d ~ 1.4585
    assert abs(sio2.nk(550.0)[0].real - 1.4599) < 2e-4


def test_tabulated_material_interpolation_and_metadata():
    with tempfile.TemporaryDirectory() as d:
        good = os.path.join(d, "good.csv")
        with open(good, "w") as f:
            f.write("# material: test\n# form: synthetic\n# wavelength_unit: um\n# source: unit test\n"
                    "# reference: none\n# limitations: synthetic data\nwavelength,n,k\n0.4,2.0,0.1\n0.8,3.0,0.3\n")
        m = load_table(good)
        v = m.nk(np.array([600.0]))[0]
        assert abs(v - (2.5 + 0.2j)) < 1e-12
        raises(ValueError, m.nk, np.array([900.0]))
        assert abs(load_table(good, allow_extrapolation=True).nk(np.array([900.0]))[0] - (3.0 + 0.3j)) < 1e-12
        bad = os.path.join(d, "bad.csv")
        with open(bad, "w") as f:
            f.write("# material: test\nwavelength,n\n400,2.0\n800,3.0\n")
        raises(ValueError, load_table, bad)


def test_csv_export_roundtrip():
    with tempfile.TemporaryDirectory() as d:
        res = solve_structure(WL, air, [], sio2, 0, "TE")
        p = os.path.join(d, "o.csv"); write_csv(p, res)
        a = np.loadtxt(p, delimiter=",", skiprows=1)
        assert a.shape == (len(WL), 4) and np.array_equal(a[:, 1], res.R)


def test_case_file_roundtrip():
    here = os.path.join(os.path.dirname(__file__), "..", "..", "simulations", "validation_cases")
    c = read_case(os.path.join(here, "dbr_N8_TE_0.txt"))
    assert len(c.layers) == 16 and c.polarization == "TE"
    res = run_case(c)
    assert np.max(np.abs(res.R + res.T - 1)) < 1e-12
