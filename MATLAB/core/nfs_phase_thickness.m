function delta = nfs_phase_thickness(lam_nm, eta, d_nm)
%NFS_PHASE_THICKNESS  delta = 2*pi*eta*d/lambda  (eta = N cos(theta)).
delta = 2*pi*eta*d_nm./lam_nm;
end
