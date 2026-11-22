function test_core()
%TEST_CORE  Phase-2 validation tests (mirror of Python/tests/test_core.py).
%  NOT YET EXECUTED in the build environment (no MATLAB/Octave available there).
here = fileparts(mfilename('fullpath')); addpath(genpath(fileparts(here)));
global NFS_FAIL NFS_TOTAL; NFS_FAIL = 0; NFS_TOTAL = 0;
lam = linspace(450,800,141); L0 = 550;
nH = real(nfs_material_nk('TiO2',L0)); nL = real(nfs_material_nk('SiO2',L0));
dH = L0/(4*nH); dL = L0/(4*nL);
mk = @(v) complex(v*ones(size(lam)));

check(abs(nfs_phase_thickness(550,1.5+0i,110) - 2*pi*1.5*110/550) < 1e-14, 'phase thickness');

N = 2+0i; s30 = sin(30*pi/180); c = sqrt(1-(s30/2)^2);
check(abs(nfs_admittance(N,nfs_eta(N,0),'TM') - 2) < 1e-14, 'admittance normal');
check(abs(nfs_admittance(N,nfs_eta(N,s30),'TE') - 2*c) < 1e-14 && ...
      abs(nfs_admittance(N,nfs_eta(N,s30),'TM') - 2/c) < 1e-14, 'admittance oblique');

M = nfs_layer_matrix(0.7,1.7); check(abs(det(M)-1) < 1e-14, 'layer matrix det=1');
Mq = nfs_layer_matrix(pi/2,1.7); check(abs(Mq(1,2)-1i/1.7) < 1e-14 && abs(Mq(2,1)-1i*1.7) < 1e-14, 'quarter-wave matrix');

a = nfs_tmm_solve(600,1+0i,{2.4+0i,1.46+0i},[60 120],1.7+0i,0,'TE');
b = nfs_tmm_solve(600,1+0i,{1.46+0i,2.4+0i},[120 60],1.7+0i,0,'TE');
check(abs(a.R-b.R) > 1e-3, 'layer order matters');

ok = true;
for pol = {'TE','TM'}
    for th = [0 30 60 85]
        r = nfs_solve_structure(lam,'air',{},[],'SiO2',th,pol{1});
        ns = real(nfs_material_nk('SiO2',lam));
        ref = arrayfun(@(n) abs(fresnel(1,n,th,pol{1}))^2, ns);
        ok = ok && max(abs(r.R-ref)) < 1e-13 && max(abs(r.R+r.T-1)) < 1e-13;
    end
end
check(ok, 'Fresnel single interface');

r = nfs_solve_structure(lam,'air',{},[],'const:1.5',0,'TE'); check(max(abs(r.R-0.04)) < 1e-14, 'empty layer list');

ns = 1.52; nl = sqrt(ns);
r = nfs_tmm_solve(600,1+0i,{nl+0i},600/(4*nl),ns+0i,0,'TE'); check(r.R < 1e-28, 'quarter-wave AR');
r = nfs_tmm_solve(600,1+0i,{2.3+0i},600/(2*2.3),1.5+0i,0,'TE'); check(abs(r.R-0.04) < 1e-14, 'half-wave absent');

specs = repmat({'TiO2','SiO2'},1,12); dd = repmat([dH dL],1,12); ok = true;
for pol = {'TE','TM'}
    for th = [0 20 45 70]
        r = nfs_solve_structure(lam,'air',specs,dd,'SiO2',th,pol{1}); ok = ok && max(abs(r.R+r.T-1)) < 1e-12;
    end
end
check(ok, 'energy conservation lossless DBR');

ok = true;
for pol = {'TE','TM'}
    r = nfs_solve_structure(lam,'air',{'const:2.4,0.05'},100,'SiO2',40,pol{1}); ok = ok && min(r.A) > 1e-4;
end
check(ok, 'absorptance positive for k>0');
r = nfs_solve_structure(600,'air',{'const:2.0,0.5'},5000,'SiO2',0,'TE');
check(abs(r.R - abs((1-(2+0.5i))/(1+(2+0.5i)))^2) < 1e-12 && r.T < 1e-12, 'thick absorber limit');

ok = true;
for pol = {'TE','TM'}
    r = nfs_solve_structure(lam,'const:1.5',{},[],'air',60,pol{1}); ok = ok && max(abs(r.R-1)) < 1e-13 && max(abs(r.T)) < 1e-13;
    f = nfs_solve_structure(lam,'const:1.5',{'air'},2000,'const:1.5',50,pol{1});
    ok = ok && max(abs(f.R+f.T-1)) < 1e-12 && all(f.T >= 0) && max(f.T) < 1e-4;
end
check(ok, 'TIR and frustrated TIR');

ok = true;
for Nn = [1 4 8 15]
    r = nfs_solve_structure(L0,'air',repmat({'TiO2','SiO2'},1,Nn),repmat([dH dL],1,Nn),'SiO2',0,'TE');
    Y = (nH/nL)^(2*Nn)*nL; ok = ok && abs(r.R - ((1-Y)/(1+Y))^2) < 1e-12;
end
check(ok, 'quarter-wave DBR analytic peak');

mir = repmat({'TiO2','SiO2'},1,5); dm = repmat([dH dL],1,5);
w = linspace(540,560,4001);
r = nfs_solve_structure(w,'SiO2',[mir {'SiO2'} fliplr(mir)],[dm 2*L0/(2*nL) fliplr(dm)],'SiO2',0,'TE');
[Tm,im] = max(r.T); check(Tm > 0.9999 && abs(w(im)-550) < 0.05, 'FP symmetric resonance');

check(throws(@() nfs_tmm_solve(0,1+0i,{},[],1+0i)), 'invalid wavelength 0');
check(throws(@() nfs_tmm_solve(-5,1+0i,{},[],1+0i)), 'invalid wavelength <0');
check(throws(@() nfs_tmm_solve(NaN,1+0i,{},[],1+0i)), 'invalid wavelength NaN');
check(throws(@() nfs_tmm_solve([],1+0i,{},[],1+0i)), 'empty wavelengths');
check(throws(@() nfs_tmm_solve(500,1+0i,{1+0i},-10,1+0i)), 'negative thickness');
check(throws(@() nfs_tmm_solve(500,1+0i,{1+0i},NaN,1+0i)), 'NaN thickness');
check(throws(@() nfs_tmm_solve(500,1+0i,{},[],1+0i,90,'TE')), 'angle 90');
check(throws(@() nfs_tmm_solve(500,1+0i,{},[],1+0i,0,'XX')), 'bad polarization');
check(throws(@() nfs_tmm_solve(500,1+0.1i,{},[],1+0i)), 'lossy incident');
check(throws(@() nfs_material_nk('TiO2',400)), 'TiO2 below range');
check(real(nfs_material_nk('TiO2',400,true)) > 2.5, 'TiO2 explicit extrapolation');
check(throws(@() nfs_material_nk('unobtainium',500)), 'unknown material');
check(throws(@() nfs_material_nk('const:-1',500)), 'negative n');
check(throws(@() nfs_material_nk('const:1,2,3',500)), 'bad const spec');
check(abs(real(nfs_material_nk('SiO2',587.5618)) - 1.45846) < 2e-4, 'SiO2 reference index');

tf = [tempname '.csv']; fid = fopen(tf,'w');
fprintf(fid,'# material: test\n# form: synthetic\n# wavelength_unit: um\n# source: unit test\n# reference: none\n# limitations: synthetic\nwavelength,n,k\n0.4,2.0,0.1\n0.8,3.0,0.3\n'); fclose(fid);
check(abs(nfs_material_nk(['table:' tf],600) - (2.5+0.2i)) < 1e-12, 'table interpolation');
check(throws(@() nfs_material_nk(['table:' tf],900)), 'table out of range');
check(abs(nfs_material_nk(['table:' tf],900,true) - (3.0+0.3i)) < 1e-12, 'table explicit extrapolation');
delete(tf);

tc = [tempname '.csv']; r = nfs_solve_structure(lam,'air',{},[],'const:1.5',0,'TE'); nfs_write_csv(tc,r);
d = readmatrix(tc); check(isequal(size(d),[numel(lam) 4]) && max(abs(d(:,2).'-r.R)) == 0, 'CSV export'); delete(tc);

fprintf('\n%d/%d checks passed\n', NFS_TOTAL-NFS_FAIL, NFS_TOTAL);
if NFS_FAIL > 0, error('test_core: %d checks failed', NFS_FAIL); end
end

function check(cond, name)
global NFS_FAIL NFS_TOTAL; NFS_TOTAL = NFS_TOTAL + 1;
if cond, fprintf('PASS  %s\n',name); else, NFS_FAIL = NFS_FAIL + 1; fprintf('FAIL  %s\n',name); end
end

function r = throws(f)
r = false; try, f(); catch, r = true; end
end

function r = fresnel(n0, ns, th, pol)
s = n0*sin(th*pi/180); c0 = sqrt(complex(1-(s/n0)^2)); cs = sqrt(complex(1-(s/ns)^2));
if strcmp(pol,'TE'), r = (n0*c0-ns*cs)/(n0*c0+ns*cs); else, r = (ns*c0-n0*cs)/(ns*c0+n0*cs); end
end
