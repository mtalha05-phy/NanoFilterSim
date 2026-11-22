function c = nfs_read_case(path)
%NFS_READ_CASE  Parse the shared plain-text case format (see Python/nanosim/core/case.py).
c = struct('name','case','incident','air','substrate','SiO2','polarization','TE','theta',0, ...
           'wl',[400 800 401],'layerSpecs',{{}},'d',[]);
fid = fopen(path,'r'); if fid < 0, error('nfs:io','cannot open %s',path); end
cleaner = onCleanup(@() fclose(fid));
while true
    ln = fgetl(fid); if ~ischar(ln), break; end
    s = strtrim(ln); if isempty(s) || s(1) == '#', continue; end
    sp = find(isspace(s),1); if isempty(sp), key = s; rest = ''; else, key = s(1:sp-1); rest = strtrim(s(sp+1:end)); end
    switch key
        case 'name', c.name = rest;
        case 'incident', c.incident = rest;
        case 'substrate', c.substrate = rest;
        case 'polarization', c.polarization = upper(rest);
        case 'theta', c.theta = str2double(rest);
        case 'wavelength', c.wl = sscanf(rest,'%f').';
        case 'layer'
            p = find(isspace(rest),1,'last');
            c.layerSpecs{end+1} = strtrim(rest(1:p-1)); c.d(end+1) = str2double(rest(p+1:end));
        otherwise, error('nfs:case','unknown key %s',key);
    end
end
end
