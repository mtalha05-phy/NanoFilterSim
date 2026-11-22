# Transfer Matrix Method — conventions and code map

| Quantity | Definition | Python (`Python/nanosim/core/tmm.py`) | MATLAB (`MATLAB/core/`) | C++ (`Cpp/src/tmm.cpp`) |
|---|---|---|---|---|
| time convention | e^{+i w t} (Macleod); solver index N = n - ik (database n + ik) | `internal_index` | `conj` in `nfs_tmm_solve` | `std::conj` in `solve` |
| normal component | eta = N cos(th) = sqrt(N^2 - s^2), s = n0 sin(th0), Im(eta) <= 0 | `normal_eta` | `nfs_eta` | `normalEta` |
| admittance | q_TE = eta, q_TM = N^2/eta | `admittance` | `nfs_admittance` | `admittance` |
| phase thickness | delta = 2 pi eta d / lambda | `phase_thickness` | `nfs_phase_thickness` | `phaseThickness` |
| layer matrix | [[cos d, i sin d/q],[i q sin d, cos d]] | `layer_matrix` | `nfs_layer_matrix` | `layerMatrix` |
| stack | M = M_1 M_2 ... M_N, layer 1 on the incident side | loop in `solve` | loop in `nfs_tmm_solve` | loop in `solve` |
| coefficients | B = M11 + M12 qs, C = M21 + M22 qs; r = (q0 B - C)/(q0 B + C); t = 2 q0/(q0 B + C) | `solve` | same | same |
| R, T, A | R = abs(r)^2; T = Re(qs)/Re(q0) abs(t)^2; A = 1 - R - T (never clipped) | | | |

Branch choice: for a passive medium (k >= 0) the wave must decay or be outgoing. With N = n - ik this selects
Im(eta) <= 0; a lossless evanescent wave (total internal reflection, frustrated TIR) therefore gets eta = -i|eta|.
Verified by the TIR (R = 1) and frustrated-TIR tests (R + T = 1, T decays with gap) and by A > 0 for k > 0.

Assumptions/limits: incident medium lossless; 0 <= theta < 90 deg; no optical activity, anisotropy, roughness or
scattering; very thick strongly absorbing layers can overflow cos/sin (raised as an error, not hidden).
Analytic checks used: Fresnel coefficients, quarter-wave AR layer, half-wave "absent" layer, exact quarter-wave stack
admittance Y = (nH/nL)^(2N) ns, thick-absorber limit. The (4/pi) asin stopband formula is an idealized approximation
and is tested only as a loose consistency check.
