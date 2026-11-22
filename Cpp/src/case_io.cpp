#include "nanofilter/case_io.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace nfs {
static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
    return a == std::string::npos ? "" : s.substr(a, b - a + 1);
}

Case readCase(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open case file: " + path);
    Case c; std::string raw;
    while (std::getline(f, raw)) {
        std::string line = trim(raw);
        if (line.empty() || line[0] == '#') continue;
        size_t sp = line.find(' ');
        std::string key = line.substr(0, sp), rest = sp == std::string::npos ? "" : trim(line.substr(sp + 1));
        if (key == "name") c.name = rest;
        else if (key == "incident") c.incident = rest;
        else if (key == "substrate") c.substrate = rest;
        else if (key == "polarization") c.polarization = rest;
        else if (key == "theta") c.theta = std::stod(rest);
        else if (key == "wavelength") { std::stringstream ss(rest); ss >> c.wlMin >> c.wlMax >> c.wlN; }
        else if (key == "layer") {
            size_t p = rest.find_last_of(" \t");
            c.layers.emplace_back(trim(rest.substr(0, p)), std::stod(rest.substr(p + 1)));
        } else throw std::runtime_error("unknown key '" + key + "'");
    }
    return c;
}

std::vector<double> wavelengthGrid(const Case& c) {
    if (c.wlN < 1) throw std::invalid_argument("npoints must be >= 1");
    std::vector<double> w(c.wlN);
    if (c.wlN == 1) { w[0] = c.wlMin; return w; }
    double step = (c.wlMax - c.wlMin) / (c.wlN - 1);
    for (int i = 0; i < c.wlN; ++i) w[i] = c.wlMin + i * step;
    return w;
}

Result runCase(const Case& c, bool allowExtrap) {
    Structure st;
    st.incident = Material::fromSpec(c.incident, allowExtrap);
    st.substrate = Material::fromSpec(c.substrate, allowExtrap);
    for (auto& l : c.layers) st.layers.emplace_back(Material::fromSpec(l.first, allowExtrap), l.second);
    return solve(st, wavelengthGrid(c), c.theta, parsePol(c.polarization));
}

void writeCsv(const std::string& path, const Result& r) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f << "wavelength_nm,R,T,A\n";
    char buf[256];
    for (size_t i = 0; i < r.wl.size(); ++i) {
        std::snprintf(buf, sizeof buf, "%.17g,%.17g,%.17g,%.17g\n", r.wl[i], r.R[i], r.T[i], r.A[i]);
        f << buf;
    }
}
}  // namespace nfs
