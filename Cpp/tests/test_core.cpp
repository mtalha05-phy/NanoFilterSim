// Self-contained C++ validation tests (no framework). Exit code 0 only if every check passes.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <vector>
#include "nanofilter/case_io.hpp"
#include "nanofilter/tmm.hpp"
using namespace nfs;

static int failures = 0, total = 0;
#define CHECK(cond) do { ++total; if (!(cond)) { ++failures; std::printf("  CHECK FAILED %s:%d  %s\n", __FILE__, __LINE__, #cond); } } while (0)
static bool throws(const std::function<void()>& f) { try { f(); } catch (...) { return true; } return false; }
static void run(const char* name, const std::function<void()>& f) {
    int before = failures; f(); std::printf("%s  %s\n", failures == before ? "PASS" : "FAIL", name);
}
static std::vector<double> grid(double a, double b, int n) {
    std::vector<double> w(n); for (int i = 0; i < n; ++i) w[i] = a + i * (b - a) / (n - 1); return w;
}
static Structure make(const std::string& inc, const std::string& sub) {
    Structure s; s.incident = Material::fromSpec(inc); s.substrate = Material::fromSpec(sub); return s;
}
static cd fresnel(double n0, double ns, double th, Pol p) {
    double s = n0 * std::sin(th * M_PI / 180.0);
    cd c0 = std::sqrt(cd(1 - (s / n0) * (s / n0))), cs = std::sqrt(cd(1 - (s / ns) * (s / ns)));
    return p == Pol::TE ? (n0 * c0 - ns * cs) / (n0 * c0 + ns * cs) : (ns * c0 - n0 * cs) / (ns * c0 + n0 * cs);
}
static double maxabs(const std::vector<double>& v, double ref) { double m = 0; for (double x : v) m = std::max(m, std::fabs(x - ref)); return m; }

int main() {
    auto wl = grid(450, 800, 141);
    const double L0 = 550.0;
    double nH = Material::fromSpec("TiO2").nk(L0).real(), nL = Material::fromSpec("SiO2").nk(L0).real();
    double dH = L0 / (4 * nH), dL = L0 / (4 * nL);

    run("phase_thickness", [&] { CHECK(std::abs(phaseThickness(550, cd(1.5), 110) - 2 * M_PI * 1.5 * 110 / 550.0) < 1e-14); });

    run("admittance", [&] {
        cd N(2.0, 0);
        CHECK(std::abs(admittance(N, normalEta(N, 0), Pol::TM) - 2.0) < 1e-14);
        double s = std::sin(30 * M_PI / 180.0), c = std::sqrt(1 - (s / 2) * (s / 2));
        CHECK(std::abs(admittance(N, normalEta(N, s), Pol::TE) - 2 * c) < 1e-14);
        CHECK(std::abs(admittance(N, normalEta(N, s), Pol::TM) - 2 / c) < 1e-14);
    });

    run("layer_matrix_unimodular_quarter_wave", [&] {
        auto M = layerMatrix(0.7, cd(1.7));
        CHECK(std::abs(M[0] * M[3] - M[1] * M[2] - 1.0) < 1e-14);
        auto Mq = layerMatrix(M_PI / 2, cd(1.7));
        CHECK(std::abs(Mq[1] - cd(0, 1) / 1.7) < 1e-14 && std::abs(Mq[2] - cd(0, 1) * 1.7) < 1e-14);
    });

    run("matrix_order_incident_first", [&] {
        Structure a = make("air", "const:1.7"), b = a;
        a.layers = {{Material::constant(2.4), 60}, {Material::constant(1.46), 120}};
        b.layers = {{Material::constant(1.46), 120}, {Material::constant(2.4), 60}};
        CHECK(std::fabs(solve(a, {600}, 0, Pol::TE).R[0] - solve(b, {600}, 0, Pol::TE).R[0]) > 1e-3);
    });

    run("fresnel_single_interface", [&] {
        for (Pol p : {Pol::TE, Pol::TM}) for (double th : {0.0, 30.0, 60.0, 85.0}) {
            auto r = solve(make("air", "SiO2"), wl, th, p);
            for (size_t i = 0; i < wl.size(); ++i) {
                double ns = Material::fromSpec("SiO2").nk(wl[i]).real();
                CHECK(std::fabs(r.R[i] - std::norm(fresnel(1.0, ns, th, p))) < 1e-13);
                CHECK(std::fabs(r.R[i] + r.T[i] - 1) < 1e-13);
            }
        }
    });

    run("empty_layers_bare_interface", [&] { CHECK(maxabs(solve(make("air", "const:1.5"), wl, 0, Pol::TE).R, 0.04) < 1e-14); });

    run("quarter_wave_AR_and_half_wave_absent", [&] {
        double ns = 1.52, nl = std::sqrt(ns);
        Structure s = make("air", "const:1.52"); s.layers = {{Material::constant(nl), 600 / (4 * nl)}};
        CHECK(solve(s, {600}, 0, Pol::TE).R[0] < 1e-28);
        Structure h = make("air", "const:1.5"); h.layers = {{Material::constant(2.3), 600 / (2 * 2.3)}};
        CHECK(std::fabs(solve(h, {600}, 0, Pol::TE).R[0] - 0.04) < 1e-14);
    });

    run("energy_conservation_lossless_dbr", [&] {
        Structure s = make("air", "SiO2");
        for (int i = 0; i < 12; ++i) { s.layers.emplace_back(Material::fromSpec("TiO2"), dH); s.layers.emplace_back(Material::fromSpec("SiO2"), dL); }
        for (Pol p : {Pol::TE, Pol::TM}) for (double th : {0.0, 20.0, 45.0, 70.0}) {
            auto r = solve(s, wl, th, p);
            for (size_t i = 0; i < wl.size(); ++i) CHECK(std::fabs(r.R[i] + r.T[i] - 1) < 1e-12);
        }
    });

    run("absorption_positive_and_thick_absorber_limit", [&] {
        Structure s = make("air", "SiO2"); s.layers = {{Material::constant(2.4, 0.05), 100}};
        for (Pol p : {Pol::TE, Pol::TM}) { auto r = solve(s, wl, 40, p); for (double a : r.A) CHECK(a > 1e-4); }
        Structure t = make("air", "SiO2"); t.layers = {{Material::constant(2.0, 0.5), 5000}};
        auto r = solve(t, {600}, 0, Pol::TE);
        cd N(2.0, 0.5);
        CHECK(std::fabs(r.R[0] - std::norm((1.0 - N) / (1.0 + N))) < 1e-12 && r.T[0] < 1e-12);
    });

    run("total_internal_reflection_and_frustrated", [&] {
        for (Pol p : {Pol::TE, Pol::TM}) {
            auto r = solve(make("const:1.5", "air"), wl, 60, p);
            CHECK(maxabs(r.R, 1.0) < 1e-13 && maxabs(r.T, 0.0) < 1e-13);
            Structure s = make("const:1.5", "const:1.5"); s.layers = {{Material::constant(1.0), 2000}};
            auto f = solve(s, wl, 50, p);
            for (size_t i = 0; i < wl.size(); ++i) CHECK(std::fabs(f.R[i] + f.T[i] - 1) < 1e-12 && f.T[i] >= 0 && f.T[i] < 1e-4);
        }
    });

    run("quarter_wave_dbr_matches_analytic_peak", [&] {
        for (int N : {1, 4, 8, 15}) {
            Structure s = make("air", "SiO2");
            for (int i = 0; i < N; ++i) { s.layers.emplace_back(Material::fromSpec("TiO2"), dH); s.layers.emplace_back(Material::fromSpec("SiO2"), dL); }
            double Y = std::pow(nH / nL, 2 * N) * nL, Rref = std::pow((1 - Y) / (1 + Y), 2);
            CHECK(std::fabs(solve(s, {L0}, 0, Pol::TE).R[0] - Rref) < 1e-12);
        }
    });

    run("fabry_perot_symmetric_resonance", [&] {
        Structure s = make("SiO2", "SiO2");
        std::vector<std::pair<Material, double>> mir;
        for (int i = 0; i < 5; ++i) { mir.emplace_back(Material::fromSpec("TiO2"), dH); mir.emplace_back(Material::fromSpec("SiO2"), dL); }
        s.layers = mir; s.layers.emplace_back(Material::fromSpec("SiO2"), 2 * L0 / (2 * nL));
        for (auto it = mir.rbegin(); it != mir.rend(); ++it) s.layers.push_back(*it);
        auto w = grid(540, 560, 4001); auto r = solve(s, w, 0, Pol::TE);
        size_t im = 0; for (size_t i = 0; i < w.size(); ++i) if (r.T[i] > r.T[im]) im = i;
        CHECK(r.T[im] > 0.9999 && std::fabs(w[im] - 550.0) < 0.05);
    });

    run("invalid_inputs", [&] {
        Structure s = make("air", "const:1.5");
        CHECK(throws([&] { solve(s, {0.0}, 0, Pol::TE); }));
        CHECK(throws([&] { solve(s, {-5.0}, 0, Pol::TE); }));
        CHECK(throws([&] { solve(s, {NAN}, 0, Pol::TE); }));
        CHECK(throws([&] { solve(s, {}, 0, Pol::TE); }));
        CHECK(throws([&] { solve(s, {500}, 90.0, Pol::TE); }));
        Structure neg = s; neg.layers = {{Material::constant(2.0), -10}};
        CHECK(throws([&] { solve(neg, {500}, 0, Pol::TE); }));
        Structure nan = s; nan.layers = {{Material::constant(2.0), NAN}};
        CHECK(throws([&] { solve(nan, {500}, 0, Pol::TE); }));
        Structure lossy = make("const:1.0,0.1", "const:1.5");
        CHECK(throws([&] { solve(lossy, {500}, 0, Pol::TE); }));
        CHECK(throws([&] { parsePol("XX"); }));
    });

    run("material_range_and_definition", [&] {
        CHECK(throws([&] { Material::fromSpec("TiO2").nk(400.0); }));
        CHECK(Material::fromSpec("TiO2", true).nk(400.0).real() > 2.5);
        CHECK(throws([&] { Material::fromSpec("unobtainium"); }));
        CHECK(throws([&] { Material::constant(-1.0); }));
        CHECK(throws([&] { Material::constant(1.5, -0.1); }));
        CHECK(throws([&] { Material::fromSpec("const:1,2,3"); }));
        CHECK(std::fabs(Material::fromSpec("SiO2").nk(587.5618).real() - 1.45846) < 2e-4);
    });

    run("table_material", [&] {
        const char* p = "/tmp/nfs_table_test.csv";
        { std::ofstream f(p); f << "# material: test\n# form: synthetic\n# wavelength_unit: um\n# source: unit test\n# reference: none\n# limitations: synthetic\nwavelength,n,k\n0.4,2.0,0.1\n0.8,3.0,0.3\n"; }
        auto m = Material::table(p);
        CHECK(std::abs(m.nk(600) - cd(2.5, 0.2)) < 1e-12);
        CHECK(throws([&] { m.nk(900); }));
        CHECK(std::abs(Material::table(p, true).nk(900) - cd(3.0, 0.3)) < 1e-12);
        { std::ofstream f("/tmp/nfs_table_bad.csv"); f << "# material: x\nwavelength,n\n400,2\n800,3\n"; }
        CHECK(throws([&] { Material::table("/tmp/nfs_table_bad.csv"); }));
    });

    run("csv_export", [&] {
        auto r = solve(make("air", "const:1.5"), wl, 0, Pol::TE);
        writeCsv("/tmp/nfs_out.csv", r);
        std::ifstream f("/tmp/nfs_out.csv"); std::string h; std::getline(f, h);
        CHECK(h == "wavelength_nm,R,T,A");
        double w, R, T, A; char c; f >> w >> c >> R >> c >> T >> c >> A;
        CHECK(R == r.R[0] && T == r.T[0]);
    });

    std::printf("\n%d/%d checks passed\n", total - failures, total);
    return failures ? 1 : 0;
}
