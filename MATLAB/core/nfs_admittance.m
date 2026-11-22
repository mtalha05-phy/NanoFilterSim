function q = nfs_admittance(N, eta, pol)
%NFS_ADMITTANCE  q_TE = eta,  q_TM = N^2/eta.
if any(abs(eta(:)) < 1e-14), error('nfs:grazing','zero normal wavevector component'); end
switch upper(pol)
    case 'TE', q = eta;
    case 'TM', q = N.*N./eta;
    otherwise, error('nfs:badPol','polarization must be TE or TM');
end
end
