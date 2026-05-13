% setup.m - Add all function directories to the MATLAB path
%
% Run this script once from the repository root directory before
% executing any figure scripts.
%
% Usage:
%   >> cd /path/to/path-independence-ebm
%   >> setup
%   >> cd paper1_origin_of_path_independence/figures
%   >> PaperFig1

reporoot = fileparts(mfilename('fullpath'));

% Shared model functions (two-box EBM, forcing, CO2 concentration)
addpath(fullfile(reporoot, 'shared'));

% Paper 1 helper functions (CO2 emissions scenario variants)
addpath(fullfile(reporoot, 'paper1_origin_of_path_independence', 'helpers'));

% Paper 2 helper functions (single-forcer EBM, concentration model)
addpath(fullfile(reporoot, 'paper2_cumulative_emissions_accounting', 'helpers'));

fprintf('Paths added. cd to a figures/ directory and run any script.\n');
