#pragma once
#include <string>
#include <vector>
#include "nanofilter/tmm.hpp"

namespace nfs {
struct Case {
    std::string name = "case", incident = "air", substrate = "SiO2", polarization = "TE";
    double theta = 0.0, wlMin = 400, wlMax = 800;
    int wlN = 401;
    std::vector<std::pair<std::string, double>> layers;
};
Case readCase(const std::string& path);
std::vector<double> wavelengthGrid(const Case& c);
Result runCase(const Case& c, bool allowExtrapolation = false);
void writeCsv(const std::string& path, const Result& r);
}  // namespace nfs
