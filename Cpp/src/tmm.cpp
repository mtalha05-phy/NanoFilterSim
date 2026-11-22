#include "nanofilter/tmm.hpp"
#include <cmath>
#include <stdexcept>

namespace nfs {
static const double PI = 3.14159265358979323846;

cd normalEta(cd N, double s) {
    cd eta = std::sqrt(N * N - s * s);
    return eta.imag() > 0 ? -eta : eta;
}

cd admittance(cd N, cd eta, Pol pol) {
    if (std::abs(eta) < 1e-14)
        throw std::invalid_argument("zero normal wavevector component (grazing/critical)");
    return pol == Pol::TE ? eta : N * N / eta;
}

cd phaseThickness(double wl, cd eta, double d) { return 2.0 * PI * eta * d / wl; }

std::array<cd, 4> layerMatrix(cd delta, cd q) {
    cd c = std::cos(delta), s = std::sin(delta);
    return {c, cd(0, 1) * s / q, cd(0, 1) * q * s, c};
}

std::array<cd, 4> matmul(const std::array<cd, 4>& a, const std::array<cd, 4>& b) {
    return {a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3],
            a[2] * b[0] + a[3] * b[2], a[2] * b[1] + a[3] * b[3]};
}

Pol parsePol(const std::string& s) {
    if (s == "TE" || s == "te") return Pol::TE;
    if (s == "TM" || s == "tm") return Pol::TM;
    throw std::invalid_argument("polarization must be TE or TM");
}

Result solve(const Structure& st, const std::vector<double>& wl, double theta_deg, Pol pol) {
    if (wl.empty()) throw std::invalid_argument("wavelength array is empty");
    if (!(std::isfinite(theta_deg) && theta_deg >= 0.0 && theta_deg < 90.0))
        throw std::invalid_argument("incidence angle must satisfy 0 <= theta < 90 deg");
    for (auto& l : st.layers)
        if (!std::isfinite(l.second) || l.second < 0) throw std::invalid_argument("layer thickness must be finite and >= 0");
    for (double w : wl)
        if (!std::isfinite(w) || w <= 0) throw std::invalid_argument("wavelengths must be finite and > 0");

    Result out;
    const double sinT = std::sin(theta_deg * PI / 180.0);
    for (double lam : wl) {
        cd n0c = st.incident.nk(lam);
        if (std::abs(n0c.imag()) > 1e-12) throw std::invalid_argument("incident medium must be lossless");
        double n0 = n0c.real(), s = n0 * sinT;
        cd q0 = admittance(cd(n0, 0), normalEta(cd(n0, 0), s), pol);

        std::array<cd, 4> M = {cd(1), cd(0), cd(0), cd(1)};
        for (auto& l : st.layers) {
            cd N = std::conj(l.first.nk(lam));
            cd eta = normalEta(N, s);
            M = matmul(M, layerMatrix(phaseThickness(lam, eta, l.second), admittance(N, eta, pol)));
        }
        cd Ns = std::conj(st.substrate.nk(lam));
        cd qs = admittance(Ns, normalEta(Ns, s), pol);

        cd B = M[0] + M[1] * qs, C = M[2] + M[3] * qs, den = q0 * B + C;
        if (!std::isfinite(den.real()) || !std::isfinite(den.imag()) || std::abs(den) == 0.0)
            throw std::runtime_error("overflow in transfer matrix");
        cd r = (q0 * B - C) / den, t = 2.0 * q0 / den;
        double R = std::norm(r), T = (qs.real() / q0.real()) * std::norm(t);
        out.wl.push_back(lam); out.R.push_back(R); out.T.push_back(T); out.A.push_back(1.0 - R - T);
        out.r.push_back(r); out.t.push_back(t);
    }
    return out;
}

}  // namespace nfs
