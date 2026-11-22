# Material database

Complex index convention in files and APIs: **n + ik, k >= 0**, wavelengths in **nm**.

## Tabulated datasets (`table:<path>` spec)
CSV, metadata as `# key: value` lines, then `wavelength,n[,k]`. **Required keys:** `material`, `wavelength_unit`
(`nm` or `um`), `source`, `reference`, `limitations`; `form` recommended. Missing keys => load error.
Out-of-range wavelengths raise an error unless extrapolation is explicitly enabled (table edge value is held).

## Built-in dispersion models (Phase 2)
| name | model | validity | notes |
|---|---|---|---|
| `SiO2` / `fused_silica` | Sellmeier, Malitson, J. Opt. Soc. Am. 55, 1205 (1965) | 210-3710 nm | bulk glass, k = 0 assumed |
| `TiO2` | DeVore, J. Opt. Soc. Am. 41, 416 (1951), rutile ordinary ray: n^2 = 5.913 + 0.2441/(l^2 - 0.0803) | **430-1530 nm** | bulk crystal (n(550) ~ 2.65), not a thin film; k = 0 assumed |

**Status: coefficients were entered from the published values without network access and have NOT yet been
checked against refractiveindex.info.** Verify (and ideally replace TiO2 by a thin-film dataset) before any
publication. The default 400-800 nm window is not fully covered by the TiO2 model (430 nm lower bound); the library
refuses to extrapolate by default, so validation cases use 450-800 nm.

Other specs: `air`, `const:n`, `const:n,k` (for analytic tests). Empty folders `TiO2/ SiO2/ Si3N4/ MgF2/` are reserved for
tabulated datasets with their source files.
