function T = nfs_load_table(path)
%NFS_LOAD_TABLE  CSV with '# key: value' metadata then columns wavelength,n[,k].
fid = fopen(path,'r'); if fid < 0, error('nfs:io','cannot open %s',path); end
cleaner = onCleanup(@() fclose(fid));
meta = struct(); rows = [];
while true
    ln = fgetl(fid); if ~ischar(ln), break; end
    s = strtrim(ln); if isempty(s), continue; end
    if s(1) == '#'
        c = find(s == ':', 1);
        if ~isempty(c)
            key = lower(strtrim(s(2:c-1))); key = regexprep(key,'[^a-z0-9_]','_');
            meta.(key) = strtrim(s(c+1:end));
        end
        continue
    end
    v = str2double(strsplit(s, ','));
    if numel(v) >= 2 && ~any(isnan(v)), rows(end+1,1:numel(v)) = v; end %#ok<AGROW>
end
req = {'material','wavelength_unit','source','reference','limitations'};
for i = 1:numel(req)
    if ~isfield(meta,req{i}) || isempty(meta.(req{i})), error('nfs:meta','missing metadata: %s',req{i}); end
end
switch lower(meta.wavelength_unit)
    case 'nm', sc = 1; case {'um','µm'}, sc = 1e3;
    otherwise, error('nfs:meta','unsupported wavelength_unit (use nm or um)');
end
if size(rows,1) < 2, error('nfs:table','table needs >= 2 rows'); end
T.wl = rows(:,1).*sc; T.n = rows(:,2);
if size(rows,2) > 2, T.k = rows(:,3); else, T.k = zeros(size(T.n)); end
if any(diff(T.wl) <= 0), error('nfs:table','wavelengths must be strictly increasing'); end
if any(T.n <= 0) || any(T.k < 0), error('nfs:table','need n > 0 and k >= 0'); end
T.meta = meta;
end
