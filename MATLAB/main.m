function main(action)
%MAIN  NanoFilterSim MATLAB entry point.  main('test') | main('validate') | main('all')
%  Phase 2 provides the TMM core only; filter modules arrive in later phases.
if nargin < 1, action = 'all'; end
here = fileparts(mfilename('fullpath')); addpath(genpath(here));
if any(strcmp(action,{'test','all'})), test_core(); end
if any(strcmp(action,{'validate','all'})), run_validation_cases(); end
end
