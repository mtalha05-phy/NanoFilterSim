function M = nfs_layer_matrix(delta, q)
%NFS_LAYER_MATRIX  [cos d, i sin d/q; i q sin d, cos d]  (scalar delta, q).
M = [cos(delta), 1i*sin(delta)/q; 1i*q*sin(delta), cos(delta)];
end
