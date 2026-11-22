// NanoFilterSim C++ materials. Complex index returned as n + ik (k >= 0). Wavelengths in nm.
#pragma once
#include <complex>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace nfs {

class Material {
public:
    enum class Kind { Constant, Dispersion, Table };
    std::string name;
    std::map<std::string, std::string> meta;

    static Material constant(double n, double k = 0.0, std::string name = "const");
    static Material dispersion(std::string name, double A,
                               std::vector<std::pair<double, double>> lamTerms,
                               std::vector<std::pair<double, double>> plainTerms,
                               double lo_nm, double hi_nm);
    static Material table(const std::string& path, bool allowExtrapolation = false);
    // 'air' | 'const:n[,k]' | 'SiO2' | 'fused_silica' | 'TiO2' | 'table:<path>'
    static Material fromSpec(const std::string& spec, bool allowExtrapolation = false);

    // throws std::invalid_argument for bad wavelengths, std::out_of_range outside validity range
    std::complex<double> nk(double lambda_nm) const;
    void setAllowExtrapolation(bool v) { allowExtrap_ = v; }

private:
    Kind kind_ = Kind::Constant;
    double n_ = 1.0, k_ = 0.0, A_ = 1.0, lo_ = 0.0, hi_ = 1e300;
    std::vector<std::pair<double, double>> lamTerms_, plainTerms_;
    std::vector<double> tw_, tn_, tk_;
    bool allowExtrap_ = false;
};

}  // namespace nfs
