function eta = nfs_eta(N, s)
%NFS_ETA  eta = N cos(theta) = sqrt(N^2 - s^2), branch with Im(eta) <= 0 (passive medium).
%  N = n - ik (solver convention, e^{+i w t});  s = n0 sin(theta0).
eta = sqrt(complex(N.*N - s.*s));
flip = imag(eta) > 0;
eta(flip) = -eta(flip);
end
