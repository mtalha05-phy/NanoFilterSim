function res = nfs_tmm_solve(lam, nkInc, nkLayers, d, nkSub, thetaDeg, pol)
%NFS_TMM_SOLVE  Transfer-matrix solver. Conventions identical to Python/C++ (see Python core/tmm.py).
%  lam [nm]; nkInc, nkSub: complex n+ik arrays on lam; nkLayers: cell of such arrays (incident side
%  first); d: thicknesses [nm]. M = M_1 M_2 ... M_N. R, T, A are never clipped.
if nargin < 6, thetaDeg = 0; end
if nargin < 7, pol = 'TE'; end
lam = lam(:).'; pol = upper(pol);
if isempty(lam) || any(~isfinite(lam)) || any(lam <= 0), error('nfs:badWavelength','wavelengths must be finite and > 0'); end
if ~any(strcmp(pol,{'TE','TM'})), error('nfs:badPol','polarization must be TE or TM'); end
if ~(isfinite(thetaDeg) && thetaDeg >= 0 && thetaDeg < 90), error('nfs:badAngle','need 0 <= theta < 90 deg'); end
if numel(d) ~= numel(nkLayers), error('nfs:badLayers','thickness/layer count mismatch'); end
if any(~isfinite(d)) || any(d < 0), error('nfs:badThickness','thickness must be finite and >= 0'); end
nkInc = nkInc(:).'; nkSub = nkSub(:).';
if any(abs(imag(nkInc)) > 1e-12), error('nfs:lossyIncident','incident medium must be lossless'); end
n0 = real(nkInc); s = n0*sin(thetaDeg*pi/180);
q0 = nfs_admittance(complex(n0), nfs_eta(complex(n0), s), pol);
Ns = conj(nkSub); qs = nfs_admittance(Ns, nfs_eta(Ns, s), pol);
L = numel(lam); R = zeros(1,L); T = R; r = complex(R); t = complex(R);
for i = 1:L
    M = eye(2);
    for j = 1:numel(d)
        N = conj(nkLayers{j}(:).'); N = N(i);
        eta = nfs_eta(N, s(i));
        M = M * nfs_layer_matrix(nfs_phase_thickness(lam(i), eta, d(j)), nfs_admittance(N, eta, pol));
    end
    B = M(1,1) + M(1,2)*qs(i); C = M(2,1) + M(2,2)*qs(i); den = q0(i)*B + C;
    if ~isfinite(den) || den == 0, error('nfs:overflow','overflow in transfer matrix'); end
    r(i) = (q0(i)*B - C)/den; t(i) = 2*q0(i)/den;
    R(i) = abs(r(i))^2; T(i) = real(qs(i))/real(q0(i))*abs(t(i))^2;
end
res = struct('wl',lam,'R',R,'T',T,'A',1-R-T,'r',r,'t',t,'theta',thetaDeg,'pol',pol);
end
