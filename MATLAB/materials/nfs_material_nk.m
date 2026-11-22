function nk = nfs_material_nk(spec, lam_nm, allowExtrap)
%NFS_MATERIAL_NK  Complex index n+ik (k>=0) for a material spec, same size as lam_nm.
%  spec: 'air' | 'const:n' | 'const:n,k' | 'SiO2' | 'fused_silica' | 'TiO2' | 'table:<path>'
%  Out-of-range wavelengths error unless allowExtrap is true (formula evaluated / table edge held).
%  Built-in coefficients: see Python/nanosim/materials/database.py (unverified vs. database).
if nargin < 3, allowExtrap = false; end
spec = strtrim(spec); low = lower(spec);
if any(~isfinite(lam_nm(:))) || any(lam_nm(:) <= 0)
    error('nfs:badWavelength','wavelength must be finite and > 0 nm');
end
if strcmp(low,'air')
    nk = complex(ones(size(lam_nm)),0); return
elseif strncmp(low,'const:',6)
    v = sscanf(strrep(spec(7:end),',',' '),'%f').';
    if ~(numel(v)==1 || numel(v)==2), error('nfs:badMaterial','bad constant spec %s',spec); end
    if numel(v)==1, v(2) = 0; end
    if ~(v(1) > 0 && v(2) >= 0), error('nfs:badMaterial','need n>0 and k>=0'); end
    nk = complex(v(1)*ones(size(lam_nm)), v(2)); return
elseif any(strcmp(low,{'sio2','fused_silica'}))
    lamT = [0.6961663 0.0684043^2; 0.4079426 0.1162414^2; 0.8974794 9.896161^2];
    nk = disp_model(lam_nm, 1.0, lamT, zeros(0,2), [210 3710], allowExtrap, low);
elseif strcmp(low,'tio2')
    nk = disp_model(lam_nm, 5.913, zeros(0,2), [0.2441 0.0803], [430 1530], allowExtrap, low);
elseif strncmp(low,'table:',6)
    T = nfs_load_table(spec(7:end));
    if ~allowExtrap && (min(lam_nm(:)) < T.wl(1) || max(lam_nm(:)) > T.wl(end))
        error('nfs:range','%s: wavelength outside table range',T.meta.material);
    end
    x = min(max(lam_nm, T.wl(1)), T.wl(end));
    nk = complex(interp1(T.wl,T.n,x,'linear'), interp1(T.wl,T.k,x,'linear'));
else
    error('nfs:badMaterial','unknown material %s',spec);
end
end

function nk = disp_model(lam_nm, A, lamT, plainT, rng, allowExtrap, name)
if ~allowExtrap && (min(lam_nm(:)) < rng(1) || max(lam_nm(:)) > rng(2))
    error('nfs:range','%s: wavelength outside validity range [%g, %g] nm',name,rng(1),rng(2));
end
l2 = (lam_nm*1e-3).^2; n2 = A*ones(size(lam_nm));
for i = 1:size(lamT,1),   n2 = n2 + lamT(i,1)*l2./(l2 - lamT(i,2)); end
for i = 1:size(plainT,1), n2 = n2 + plainT(i,1)./(l2 - plainT(i,2)); end
nk = complex(sqrt(n2), 0);
end
