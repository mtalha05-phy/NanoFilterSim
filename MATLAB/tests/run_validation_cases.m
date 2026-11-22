function run_validation_cases()
%RUN_VALIDATION_CASES  Run every shared case and write results/validation/matlab/<case>.csv.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(root,'MATLAB')));
files = dir(fullfile(root,'simulations','validation_cases','*.txt'));
out = fullfile(root,'results','validation','matlab'); if ~exist(out,'dir'), mkdir(out); end
for i = 1:numel(files)
    c = nfs_read_case(fullfile(files(i).folder, files(i).name));
    nfs_write_csv(fullfile(out, [erase(files(i).name,'.txt') '.csv']), nfs_run_case(c));
end
fprintf('wrote %d MATLAB case outputs to %s\n', numel(files), out);
end
