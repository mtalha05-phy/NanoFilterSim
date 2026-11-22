// Transfer Matrix Method core. Conventions identical to Python/nanosim/core/tmm.py:
// e^{+i w t}; solver index N = n - ik; eta = sqrt(N^2 - s^2) with Im(eta) <= 0;
// q_TE = eta, q_TM = N^2/eta; delta = 2 pi eta d / lambda; M = M_1 ... M_N (incident side first).
#pragma once
#include <array>
#include <complex>
#include <string>
#include <utility>
#include <vector>
#include "nanofilter/material.hpp"

namespace nfs {
using cd = std::complex<double>;
enum class Pol { TE, TM };

struct Structure {
    Material incident = Material::constant(1.0, 0.0, "air");
    Material substrate = Material::constant(1.5, 0.0, "substrate");
    std::vector<std::pair<Material, double>> layers;  // (material, thickness nm)
};

struct Result {
    std::vector<double> wl, R, T, A;
    std::vector<cd> r, t;
};

cd normalEta(cd N, double s);                       // N cos(theta), passive branch
cd admittance(cd N, cd eta, Pol pol);
cd phaseThickness(double wl_nm, cd eta, double d_nm);
std::array<cd, 4> layerMatrix(cd delta, cd q);      // row-major {m11,m12,m21,m22}
std::array<cd, 4> matmul(const std::array<cd, 4>& a, const std::array<cd, 4>& b);

Pol parsePol(const std::string& s);
Result solve(const Structure& st, const std::vector<double>& wl_nm, double theta_deg, Pol pol);

}  // namespace nfs
