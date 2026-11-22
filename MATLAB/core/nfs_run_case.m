function res = nfs_run_case(c, allowExtrap)
if nargin < 2, allowExtrap = false; end %#ok<NASGU>
a = c.wl(1); b = c.wl(2); n = c.wl(3);
if n < 1, error('nfs:case','npoints must be >= 1'); end
if n == 1, lam = a; else, lam = a + (0:n-1)*((b-a)/(n-1)); end
res = nfs_solve_structure(lam, c.incident, c.layerSpecs, c.d, c.substrate, c.theta, c.polarization);
end
