#include "nanofilter/material.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace nfs {

static std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
}
static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
    return a == std::string::npos ? "" : s.substr(a, b - a + 1);
}

Material Material::constant(double n, double k, std::string name) {
    if (!(std::isfinite(n) && n > 0 && std::isfinite(k) && k >= 0))
        throw std::invalid_argument("constant material needs n > 0 and k >= 0");
    Material m; m.kind_ = Kind::Constant; m.n_ = n; m.k_ = k; m.name = std::move(name);
    return m;
}

Material Material::dispersion(std::string name, double A, std::vector<std::pair<double, double>> lamTerms,
                              std::vector<std::pair<double, double>> plainTerms, double lo, double hi) {
    Material m; m.kind_ = Kind::Dispersion; m.name = std::move(name); m.A_ = A;
    m.lamTerms_ = std::move(lamTerms); m.plainTerms_ = std::move(plainTerms); m.lo_ = lo; m.hi_ = hi;
    return m;
}

Material Material::table(const std::string& path, bool allowExtrap) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open material table: " + path);
    Material m; m.kind_ = Kind::Table; m.allowExtrap_ = allowExtrap;
    std::vector<std::vector<double>> rows;
    std::string line;
    while (std::getline(f, line)) {
        std::string s = trim(line);
        if (s.empty()) continue;
        if (s[0] == '#') {
            size_t c = s.find(':');
            if (c != std::string::npos) m.meta[lower(trim(s.substr(1, c - 1)))] = trim(s.substr(c + 1));
            continue;
        }
        std::vector<double> row; std::stringstream ss(s); std::string tok; bool ok = true;
        while (std::getline(ss, tok, ',')) {
            try { size_t pos; double v = std::stod(tok, &pos); if (pos != trim(tok).size()) ok = false; row.push_back(v); }
            catch (...) { ok = false; break; }
        }
        if (ok && row.size() >= 2) rows.push_back(row);
    }
    for (const char* key : {"material", "wavelength_unit", "source", "reference", "limitations"})
        if (m.meta[key].empty()) throw std::invalid_argument(std::string("table missing metadata: ") + key);
    std::string unit = lower(m.meta["wavelength_unit"]);
    double scale = unit == "nm" ? 1.0 : (unit == "um" ? 1e3 : 0.0);
    if (scale == 0.0) throw std::invalid_argument("unsupported wavelength_unit (use nm or um)");
    if (rows.size() < 2) throw std::invalid_argument("table needs >= 2 rows");
    for (auto& r : rows) {
        m.tw_.push_back(r[0] * scale); m.tn_.push_back(r[1]); m.tk_.push_back(r.size() > 2 ? r[2] : 0.0);
        if (m.tn_.back() <= 0 || m.tk_.back() < 0) throw std::invalid_argument("table needs n > 0, k >= 0");
        if (m.tw_.size() > 1 && m.tw_.back() <= m.tw_[m.tw_.size() - 2])
            throw std::invalid_argument("table wavelengths must be strictly increasing");
    }
    m.name = m.meta["material"]; m.lo_ = m.tw_.front(); m.hi_ = m.tw_.back();
    return m;
}

Material Material::fromSpec(const std::string& spec, bool allowExtrap) {
    std::string s = trim(spec), low = lower(s);
    Material m;
    if (low == "air") m = constant(1.0, 0.0, "air");
    else if (low.rfind("const:", 0) == 0) {
        std::vector<double> v; std::stringstream ss(s.substr(6)); std::string tok;
        while (std::getline(ss, tok, ',')) v.push_back(std::stod(tok));
        if (v.size() != 1 && v.size() != 2) throw std::invalid_argument("bad constant material spec: " + spec);
        m = constant(v[0], v.size() == 2 ? v[1] : 0.0);
    } else if (low == "sio2" || low == "fused_silica") {
        const double c1 = 0.0684043, c2 = 0.1162414, c3 = 9.896161;
        m = dispersion("SiO2 (fused silica)", 1.0,
                       {{0.6961663, c1 * c1}, {0.4079426, c2 * c2}, {0.8974794, c3 * c3}}, {}, 210.0, 3710.0);
    } else if (low == "tio2") {
        m = dispersion("TiO2 (rutile, ordinary ray)", 5.913, {}, {{0.2441, 0.0803}}, 430.0, 1530.0);
    } else if (low.rfind("table:", 0) == 0) m = table(s.substr(6), allowExtrap);
    else throw std::invalid_argument("unknown material: " + spec);
    m.allowExtrap_ = allowExtrap;
    return m;
}

std::complex<double> Material::nk(double wl) const {
    if (!std::isfinite(wl) || wl <= 0) throw std::invalid_argument("wavelength must be finite and > 0 nm");
    if (kind_ != Kind::Constant && (wl < lo_ || wl > hi_) && !allowExtrap_)
        throw std::out_of_range(name + ": wavelength outside validity range");
    switch (kind_) {
        case Kind::Constant: return {n_, k_};
        case Kind::Dispersion: {
            double l2 = (wl * 1e-3) * (wl * 1e-3), n2 = A_;
            for (auto& t : lamTerms_) n2 += t.first * l2 / (l2 - t.second);
            for (auto& t : plainTerms_) n2 += t.first / (l2 - t.second);
            return {std::sqrt(n2), 0.0};
        }
        case Kind::Table: {
            size_t i = std::upper_bound(tw_.begin(), tw_.end(), wl) - tw_.begin();
            if (i == 0) return {tn_.front(), tk_.front()};
            if (i >= tw_.size()) return {tn_.back(), tk_.back()};
            double f = (wl - tw_[i - 1]) / (tw_[i] - tw_[i - 1]);
            return {tn_[i - 1] + f * (tn_[i] - tn_[i - 1]), tk_[i - 1] + f * (tk_[i] - tk_[i - 1])};
        }
    }
    return {1.0, 0.0};
}

}  // namespace nfs
