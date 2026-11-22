"""Generate the shared validation cases (simulations/validation_cases/*.txt).
Quarter-wave thicknesses are computed once here and written as numbers so that all three
languages use identical inputs."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "Python"))
from nanosim.materials.database import get_material

OUT = os.path.join(os.path.dirname(__file__), "..", "simulations", "validation_cases")
os.makedirs(OUT, exist_ok=True)
L0 = 550.0
nH = get_material("TiO2").nk(L0)[0].real
nL = get_material("SiO2").nk(L0)[0].real
dH, dL = L0 / (4 * nH), L0 / (4 * nL)


def write(name, inc, sub, pol, th, wl, layers):
    with open(os.path.join(OUT, name + ".txt"), "w") as f:
        f.write(f"name {name}\nincident {inc}\nsubstrate {sub}\npolarization {pol}\ntheta {th}\n")
        f.write(f"wavelength {wl[0]} {wl[1]} {wl[2]}\n")
        for s, d in layers:
            f.write(f"layer {s} {d:.15g}\n")


W = (450, 800, 351)
write("fresnel_normal", "air", "const:1.5", "TE", 0, W, [])
write("fresnel_45_TE", "air", "SiO2", "TE", 45, W, [])
write("fresnel_45_TM", "air", "SiO2", "TM", 45, W, [])
write("single_layer_TE_30", "air", "SiO2", "TE", 30, W, [("const:2.1", 80)])
write("single_layer_TM_30", "air", "SiO2", "TM", 30, W, [("const:2.1", 80)])
dbr = [("TiO2", dH), ("SiO2", dL)] * 8
write("dbr_N8_TE_0", "air", "SiO2", "TE", 0, W, dbr)
write("dbr_N8_TE_45", "air", "SiO2", "TE", 45, W, dbr)
write("dbr_N8_TM_45", "air", "SiO2", "TM", 45, W, dbr)
m = 2
dC = m * L0 / (2 * nL)
mirror = [("TiO2", dH), ("SiO2", dL)] * 5
write("fp_HL5_C_LH5_TE_0", "air", "SiO2", "TE", 0, (540, 560, 2001),
      mirror + [("SiO2", dC)] + mirror[::-1])
write("absorbing_layer_TE_0", "air", "SiO2", "TE", 0, W, [("const:2.4,0.05", 100)])
write("absorbing_layer_TM_40", "air", "SiO2", "TM", 40, W, [("const:2.4,0.05", 100)])
write("thick_absorber_TE_0", "air", "SiO2", "TE", 0, W, [("const:2.0,0.5", 5000)])
write("tir_prism_TE_60", "const:1.5", "air", "TE", 60, W, [])
write("tir_prism_TM_60", "const:1.5", "air", "TM", 60, W, [])
write("frustrated_tir_TE_50", "const:1.5", "const:1.5", "TE", 50, W, [("air", 100)])
write("frustrated_tir_TM_50", "const:1.5", "const:1.5", "TM", 50, W, [("air", 100)])
print("wrote", len(os.listdir(OUT)), "cases; nH(550)=%.6f nL(550)=%.6f dH=%.4f dL=%.4f nm" % (nH, nL, dH, dL))
