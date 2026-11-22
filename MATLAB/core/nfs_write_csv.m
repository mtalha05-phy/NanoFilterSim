function nfs_write_csv(path, res)
fid = fopen(path,'w'); if fid < 0, error('nfs:io','cannot write %s',path); end
cleaner = onCleanup(@() fclose(fid));
fprintf(fid,'wavelength_nm,R,T,A\n');
for i = 1:numel(res.wl)
    fprintf(fid,'%.17g,%.17g,%.17g,%.17g\n',res.wl(i),res.R(i),res.T(i),res.A(i));
end
end
