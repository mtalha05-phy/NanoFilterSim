# NanoFilterSim

Design, numerical modeling, simulation and optimization of multilayer nanophotonic structures
(DBR, Fabry-Perot narrowband, long/short-pass filters) with the Transfer Matrix Method,
implemented independently in MATLAB, Python and C++17.

## Status (honest)
**Phase 1 (foundation) and Phase 2 (core TMM solver) are implemented.** Filter modules, optimization, sensitivity
analysis, graphs, LaTeX report and presentation are **not yet implemented**; folders exist as placeholders.

| Component | Implemented | Executed & passing |
|---|---|---|
| Python TMM + materials | yes | yes - see `results/validation` and test output |
| C++17 TMM + materials | yes | yes (built with g++, -Wall -Wextra clean) |
| MATLAB TMM + materials | yes | **not run** (no MATLAB/Octave in the build environment) |
| CMake build | written | **not run** (cmake unavailable; direct g++ build used) |
| Cross-language comparison | Python vs C++ | yes; MATLAB column reads "NOT RUN" until MATLAB CSVs exist |

## Run
```bash
bash tools/run_all_tests.sh            # Python tests, C++ build/tests, cross-language report
# MATLAB (once available):  cd MATLAB; main('all');  then  python tools/cross_validate.py
```
Python only: `python Python/tests/run_tests.py` (written in pytest style, but pytest could not be installed in the build environment, so pytest itself has not been run). Single case:
`python -m nanosim.run_case simulations/validation_cases/dbr_N8_TE_0.txt out.csv` (run from `Python/`) or `build/nfs_run case.txt out.csv`.

## Layout
`Python/nanosim/{core,materials,...}`, `MATLAB/{core,materials,tests,...}`, `Cpp/{include,src,tests}`,
`simulations/validation_cases` (shared case files read by all three languages), `results/validation`, `docs/theory`,
`tools/`. Conventions and code map: `docs/theory/transfer_matrix_method.md`. Materials: `materials/README.md`.

## Known limitations
* Built-in TiO2 is bulk rutile (n(550) ~ 2.65, valid 430-1530 nm), not a thin film; coefficients not yet verified
  against the database. Do not quote results as film performance.
* Extrapolation outside material ranges is refused by default.
* Planar, isotropic, nonmagnetic, specular layers only.

## References (theory)
H. A. Macleod, *Thin-Film Optical Filters*; I. H. Malitson, JOSA 55, 1205 (1965); J. R. DeVore, JOSA 41, 416 (1951).
Full bibliography will be built in `docs/research_report/references.bib` from sources actually consulted.
