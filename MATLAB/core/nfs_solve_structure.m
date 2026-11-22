function res = nfs_solve_structure(lam, incSpec, layerSpecs, d, subSpec, thetaDeg, pol)
%NFS_SOLVE_STRUCTURE  Convenience wrapper: material specs -> index arrays -> TMM.
nkL = cell(1,numel(layerSpecs));
for j = 1:numel(layerSpecs), nkL{j} = nfs_material_nk(layerSpecs{j}, lam); end
res = nfs_tmm_solve(lam, nfs_material_nk(incSpec,lam), nkL, d, nfs_material_nk(subSpec,lam), thetaDeg, pol);
end
